# ============================================================
# MODEL 5: SITE DIFFERENCES WITHIN EACH DIEL PERIOD
#
# Points = estimated marginal means
# Error bars = 95% confidence intervals
#
# Lowercase letters compare SITES WITHIN THE SAME PERIOD.
#
# Sites sharing at least one lowercase letter within the same
# period and acoustic index are not significantly different
# based on Holm-adjusted comparisons at alpha = 0.05.
#
# Letters must not be compared across different periods.
#
# Sites are horizontally offset in this order:
# AA1, AA2, KK1, KK2, KK3
#
# No connecting lines are plotted.
# ============================================================


# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------

library(tidyverse)
library(multcompView)


# ------------------------------------------------------------
# 2. FILE PATHS
# ------------------------------------------------------------

results_dir <- paste0(
  "D:/Flor/SaudiArabia/Data/Results/",
  "SiteAndReserveComparison/"
)

# Estimated marginal means for every Site x Period combination
emm_file <- paste0(
  results_dir,
  "Model5_SiteByPeriod_EstimatedMeans_DawnMorningDusk.csv"
)

# Pairwise site comparisons within each period
site_contrast_file <- paste0(
  results_dir,
  "Model5_SiteDifferencesWithinPeriod_DawnMorningDusk_Holm.csv"
)

output_png <- paste0(
  results_dir,
  "Figure4_SiteDifferencesWithinPeriod_Model5_SymbolPlot.png"
)

output_pdf <- paste0(
  results_dir,
  "Figure4_SiteDifferencesWithinPeriod_Model5_SymbolPlot.pdf"
)


# ------------------------------------------------------------
# 3. CHECK THAT INPUT FILES EXIST
# ------------------------------------------------------------

if(!file.exists(emm_file)){
  
  stop(
    paste0(
      "The estimated marginal means file was not found:\n",
      emm_file
    )
  )
}

if(!file.exists(site_contrast_file)){
  
  stop(
    paste0(
      "The site contrast file was not found:\n",
      site_contrast_file
    )
  )
}


# ------------------------------------------------------------
# 4. READ INPUT FILES
# ------------------------------------------------------------

emm <- read.csv(
  emm_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

site_contrasts <- read.csv(
  site_contrast_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 5. CHECK REQUIRED COLUMNS
# ------------------------------------------------------------

required_emm_columns <- c(
  "Period",
  "SiteID",
  "emmean",
  "lower.CL",
  "upper.CL",
  "Index"
)

required_contrast_columns <- c(
  "contrast",
  "Period",
  "p.value",
  "Index"
)

missing_emm_columns <- setdiff(
  required_emm_columns,
  names(emm)
)

missing_contrast_columns <- setdiff(
  required_contrast_columns,
  names(site_contrasts)
)

if(length(missing_emm_columns) > 0){
  
  stop(
    paste0(
      "The EMM file is missing these columns:\n",
      paste(
        missing_emm_columns,
        collapse = ", "
      )
    )
  )
}

if(length(missing_contrast_columns) > 0){
  
  stop(
    paste0(
      "The site contrast file is missing these columns:\n",
      paste(
        missing_contrast_columns,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 6. CLEAN ESTIMATED-MEAN DATA
# ------------------------------------------------------------

emm <- emm %>%
  mutate(
    SiteID = str_replace_all(
      SiteID,
      "\\\\_",
      "_"
    ),
    
    Index = str_replace_all(
      Index,
      "\\\\_",
      "_"
    ),
    
    Period = str_trim(
      Period
    ),
    
    emmean = as.numeric(
      emmean
    ),
    
    lower.CL = as.numeric(
      lower.CL
    ),
    
    upper.CL = as.numeric(
      upper.CL
    )
  )


# ------------------------------------------------------------
# 7. CLEAN SITE CONTRAST DATA
# ------------------------------------------------------------

site_contrasts <- site_contrasts %>%
  mutate(
    contrast = str_replace_all(
      contrast,
      "\\\\_",
      "_"
    ),
    
    contrast = str_squish(
      contrast
    ),
    
    Index = str_replace_all(
      Index,
      "\\\\_",
      "_"
    ),
    
    Period = str_trim(
      Period
    ),
    
    p.value = as.numeric(
      p.value
    )
  )


# ------------------------------------------------------------
# 8. DEFINE PERIOD ORDER
# ------------------------------------------------------------

period_levels <- c(
  "Dawn",
  "Morning",
  "Dusk"
)

emm <- emm %>%
  mutate(
    Period = factor(
      Period,
      levels = period_levels
    )
  )

site_contrasts <- site_contrasts %>%
  mutate(
    Period = factor(
      Period,
      levels = period_levels
    )
  )


# ------------------------------------------------------------
# 9. SPLIT SITE CONTRASTS
#
# Expected format:
# IARR_A_NorthRK - IARR_A_SouthRK
# ------------------------------------------------------------

site_contrasts2 <- site_contrasts %>%
  separate(
    col = contrast,
    
    into = c(
      "Site1",
      "Site2"
    ),
    
    sep = "\\s+-\\s+",
    
    remove = FALSE,
    
    extra = "merge",
    
    fill = "right"
  ) %>%
  
  mutate(
    Site1 = str_trim(
      Site1
    ),
    
    Site2 = str_trim(
      Site2
    ),
    
    Comparison = paste0(
      Site1,
      "-",
      Site2
    )
  )


# ------------------------------------------------------------
# 10. CHECK THAT SITE CONTRASTS WERE SPLIT CORRECTLY
# ------------------------------------------------------------

bad_contrasts <- site_contrasts2 %>%
  filter(
    is.na(Site1) |
      is.na(Site2) |
      Site1 == "" |
      Site2 == ""
  )

if(nrow(bad_contrasts) > 0){
  
  print(
    bad_contrasts
  )
  
  stop(
    paste0(
      "One or more site contrasts could not be separated. ",
      "The expected format is 'Site1 - Site2'."
    )
  )
}


# ------------------------------------------------------------
# 11. REMOVE ACCIDENTAL DUPLICATE COMPARISONS
#
# Each comparison must be unique within:
# Index x Period
# ------------------------------------------------------------

duplicate_contrasts <- site_contrasts2 %>%
  count(
    Index,
    Period,
    Comparison
  ) %>%
  filter(
    n > 1
  )

if(nrow(duplicate_contrasts) > 0){
  
  print(
    duplicate_contrasts
  )
  
  warning(
    paste0(
      "Duplicate site comparisons were found. ",
      "Only the first occurrence will be used."
    )
  )
}

site_contrasts2 <- site_contrasts2 %>%
  distinct(
    Index,
    Period,
    Comparison,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# 12. FUNCTION TO CALCULATE SITE LETTERS
#
# Letters are calculated separately for every:
# Index x Period
#
# Therefore, letters compare sites only within one period.
# ------------------------------------------------------------

get_site_letters <- function(dat){
  
  dat <- dat %>%
    filter(
      !is.na(p.value),
      is.finite(p.value),
      !is.na(Comparison)
    )
  
  if(nrow(dat) == 0){
    
    return(
      tibble(
        SiteID = character(0),
        SiteLetter = character(0)
      )
    )
  }
  
  p_values <- dat$p.value
  
  names(p_values) <- dat$Comparison
  
  letter_result <- multcompView::multcompLetters(
    p_values,
    threshold = 0.05
  )$Letters
  
  tibble(
    SiteID = names(
      letter_result
    ),
    
    SiteLetter = unname(
      letter_result
    )
  )
}


# ------------------------------------------------------------
# 13. CALCULATE SITE LETTERS FOR EACH INDEX AND PERIOD
# ------------------------------------------------------------

site_letters <- site_contrasts2 %>%
  
  group_by(
    Index,
    Period
  ) %>%
  
  group_modify(
    ~ get_site_letters(.x)
  ) %>%
  
  ungroup() %>%
  
  mutate(
    Period = factor(
      Period,
      levels = period_levels
    )
  )


# ------------------------------------------------------------
# 14. CHECK SITE LETTER OUTPUT
# ------------------------------------------------------------

cat(
  "\nNumber of site-letter records:",
  nrow(site_letters),
  "\n"
)

print(
  site_letters %>%
    arrange(
      Index,
      Period,
      SiteID
    )
)


# ------------------------------------------------------------
# 15. JOIN SITE LETTERS TO ESTIMATED MARGINAL MEANS
# ------------------------------------------------------------

plot_dat <- emm %>%
  
  left_join(
    site_letters,
    
    by = c(
      "Index",
      "Period",
      "SiteID"
    )
  )


# ------------------------------------------------------------
# 16. CHECK LETTER JOIN
#
# Expected:
# Number without letters = 0
# ------------------------------------------------------------

cat(
  "\nNumber of EMM observations:",
  nrow(plot_dat),
  "\n"
)

cat(
  "Number with letters:",
  sum(
    !is.na(plot_dat$SiteLetter)
  ),
  "\n"
)

cat(
  "Number without letters:",
  sum(
    is.na(plot_dat$SiteLetter)
  ),
  "\n"
)

if(any(is.na(plot_dat$SiteLetter))){
  
  cat(
    "\nRows without matched site letters:\n"
  )
  
  print(
    plot_dat %>%
      filter(
        is.na(SiteLetter)
      ) %>%
      select(
        Index,
        Period,
        SiteID,
        emmean
      )
  )
}


# ------------------------------------------------------------
# 17. EXCLUDE SOUND-LEVEL DIAGNOSTIC VARIABLES
# ------------------------------------------------------------

excluded_indices <- c(
  "mean_NoiseFloor_dB",
  "mean_SignalLevel_dB",
  "mean_SNR_dB"
)

plot_dat <- plot_dat %>%
  filter(
    !Index %in% excluded_indices
  )


# ------------------------------------------------------------
# 18. DEFINE SITE LABELS AND ORDER
#
# Change AA1 and AA2 to IA1 and IA2 here if needed.
# ------------------------------------------------------------

site_labels <- c(
  "IARR_A_NorthRK"   = "AA1",
  "IARR_A_SouthRK"   = "AA2",
  "KKRR_B_Falaq"     = "KK1",
  "KKRR_A_Catchment" = "KK2",
  "KKRR_B_Wadi"      = "KK3"
)

site_order <- c(
  "AA1",
  "AA2",
  "KK1",
  "KK2",
  "KK3"
)

plot_dat <- plot_dat %>%
  mutate(
    Site = recode(
      SiteID,
      !!!site_labels,
      .default = SiteID
    ),
    
    Site = factor(
      Site,
      levels = site_order
    )
  )


# ------------------------------------------------------------
# 19. DEFINE INDEX LABELS AND PANEL ORDER
# ------------------------------------------------------------

index_labels <- c(
  "mean_ACI300"  = "ACI 300-6000 Hz",
  "mean_ACI2000" = "ACI 2000-6000 Hz",
  "mean_ADI300"  = "ADI 300-6000 Hz",
  "mean_ADI2000" = "ADI 2000-6000 Hz",
  "mean_BIO300"  = "BIO 300-6000 Hz",
  "mean_BIO2000" = "BIO 2000-6000 Hz",
  "mean_NP300"   = "NP 300-6000 Hz",
  "mean_NP2000"  = "NP 2000-6000 Hz",
  "mean_NDSI"    = "NDSI",
  "mean_A"       = "Anthrophony",
  "mean_B"       = "Biophony"
)

index_order <- c(
  "ACI 300-6000 Hz",
  "ACI 2000-6000 Hz",
  "ADI 300-6000 Hz",
  "ADI 2000-6000 Hz",
  "BIO 300-6000 Hz",
  "BIO 2000-6000 Hz",
  "NP 300-6000 Hz",
  "NP 2000-6000 Hz",
  "NDSI",
  "Anthrophony",
  "Biophony"
)

plot_dat <- plot_dat %>%
  mutate(
    IndexLabel = recode(
      Index,
      !!!index_labels,
      .default = Index
    ),
    
    IndexLabel = factor(
      IndexLabel,
      levels = index_order
    ),
    
    Period = factor(
      Period,
      levels = period_levels
    )
  )


# ------------------------------------------------------------
# 20. REMOVE INVALID OR UNRECOGNIZED ROWS
# ------------------------------------------------------------

plot_dat <- plot_dat %>%
  filter(
    !is.na(Site),
    !is.na(IndexLabel),
    !is.na(Period),
    
    !is.na(emmean),
    is.finite(emmean),
    
    !is.na(lower.CL),
    is.finite(lower.CL),
    
    !is.na(upper.CL),
    is.finite(upper.CL)
  ) %>%
  
  droplevels()


# ------------------------------------------------------------
# 21. VERIFY SITE, PERIOD, AND INDEX ORDER
# ------------------------------------------------------------

cat(
  "\nObservations by site:\n"
)

print(
  table(
    plot_dat$Site,
    useNA = "ifany"
  )
)

cat(
  "\nObservations by period:\n"
)

print(
  table(
    plot_dat$Period,
    useNA = "ifany"
  )
)

cat(
  "\nObservations by index:\n"
)

print(
  table(
    plot_dat$IndexLabel,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# CALCULATE EXACT HORIZONTAL POSITIONS
#
# Site order within every period:
# AA1, AA2, KK1, KK2, KK3
# ------------------------------------------------------------

site_offsets <- c(
  "AA1" = -0.36,
  "AA2" = -0.18,
  "KK1" =  0.00,
  "KK2" =  0.18,
  "KK3" =  0.36
)

plot_dat <- plot_dat %>%
  mutate(
    PeriodNumber = as.numeric(
      Period
    ),
    
    SiteOffset = unname(
      site_offsets[
        as.character(Site)
      ]
    ),
    
    PlotX = PeriodNumber +
      SiteOffset
  )


# ------------------------------------------------------------
# POSITION LETTERS ABOVE CONFIDENCE INTERVALS
#
# letter_offset controls the vertical distance between the
# upper confidence interval and the letter.
# ------------------------------------------------------------

letter_offset <- 0.12

plot_dat <- plot_dat %>%
  group_by(
    IndexLabel
  ) %>%
  mutate(
    panel_min = min(
      lower.CL,
      na.rm = TRUE
    ),
    
    panel_max = max(
      upper.CL,
      na.rm = TRUE
    ),
    
    panel_range = panel_max -
      panel_min,
    
    panel_range = if_else(
      is.finite(panel_range) &
        panel_range > 0,
      panel_range,
      1
    ),
    
    letter_y = upper.CL +
      letter_offset * panel_range
  ) %>%
  ungroup()


# ------------------------------------------------------------
# SITE FILL COLORS
# ------------------------------------------------------------

site_fill_colours <- c(
  "AA1" = "#2166AC",
  "AA2" = "#67A9CF",
  "KK1" = "#B35806",
  "KK2" = "#E08214",
  "KK3" = "#FDB863"
)


# ------------------------------------------------------------
# CREATE SIDE-BY-SIDE BAR GRAPH
#
# No position_dodge() is needed because PlotX contains the
# precise horizontal location of each site.
# ------------------------------------------------------------

bar_width <- 0.15

p_site_period <- ggplot(
  data = plot_dat,
  mapping = aes(
    x = PlotX,
    y = emmean,
    fill = Site
  )
) +
  
  # Estimated marginal means
  geom_col(
    width = bar_width,
    colour = "black",
    linewidth = 0.45
  ) +
  
  # 95% confidence intervals
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.045,
    linewidth = 0.70,
    colour = "black"
  ) +
  
  # Site-comparison letters
  geom_text(
    aes(
      y = letter_y,
      label = SiteLetter
    ),
    size = 4.5,
    fontface = "plain",
    colour = "black",
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    facets = vars(
      IndexLabel
    ),
    scales = "free_y",
    ncol = 3,
    drop = TRUE
  ) +
  
  # Place period labels at the center of each group
  scale_x_continuous(
    breaks = c(
      1,
      2,
      3
    ),
    labels = c(
      "Dawn",
      "Morning",
      "Dusk"
    ),
    limits = c(
      0.55,
      3.45
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  scale_fill_manual(
    values = site_fill_colours,
    breaks = site_order,
    limits = site_order,
    name = "Site",
    drop = FALSE
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.02,
        0.25
      )
    )
  ) +
  
  labs(
    x = "Diel period",
    y = "Estimated marginal mean"
  ) +
  
  guides(
    fill = guide_legend(
      order = 1,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        colour = "black",
        linewidth = 0.4
      )
    )
  ) +
  
  theme_bw(
    base_size = 15
  ) +
  
  theme(
    strip.text = element_text(
      size = 16,
      face = "bold"
    ),
    
    strip.background = element_rect(
      fill = "grey95",
      colour = "black",
      linewidth = 0.6
    ),
    
    axis.title = element_text(
      size = 16
    ),
    
    axis.text.x = element_text(
      size = 15,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 14,
      colour = "black"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 15
    ),
    
    legend.key.width = grid::unit(
      1.2,
      "cm"
    ),
    
    panel.grid.major.x = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    panel.spacing = grid::unit(
      1,
      "lines"
    ),
    
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )


# ------------------------------------------------------------
# DISPLAY FIGURE
# ------------------------------------------------------------

print(
  p_site_period
)


# ------------------------------------------------------------
# SAVE FIGURE
# ------------------------------------------------------------

ggsave(
  filename = output_png,
  plot = p_site_period,
  width = 14,
  height = 14,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = output_pdf,
  plot = p_site_period,
  width = 14,
  height = 14,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)