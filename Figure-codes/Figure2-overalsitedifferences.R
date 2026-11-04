# ============================================================
# MODEL 4: SPATIAL SITE HETEROGENEITY
#
# GAMM estimated marginal means by site
# Points = estimated marginal means
# Error bars = 95% confidence intervals
# Letters = Holm-adjusted site comparisons within each index
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

emm_file <- paste0(
  results_dir,
  "Model4_OverallSiteMeans_DawnMorningDusk.csv"
)

# Change this filename if your contrast output has another name
contrast_file <- paste0(
  results_dir,
  "Model4_OverallSiteComparisons_DawnMorningDusk_Holm.csv"
)

output_png <- paste0(
  results_dir,
  "Figure2-overallsite differences -Model4-95conf.png"
)

output_pdf <- paste0(
  results_dir,
  "Figure2-overallsite differences -Model4-95conf.pdf"
)


# ------------------------------------------------------------
# 3. READ DATA
# ------------------------------------------------------------

emm <- read.csv(
  emm_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

site_contrasts <- read.csv(
  contrast_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. CHECK REQUIRED COLUMNS
# ------------------------------------------------------------

required_emm_columns <- c(
  "SiteID",
  "emmean",
  "lower.CL",
  "upper.CL",
  "Index"
)

required_contrast_columns <- c(
  "contrast",
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

if (length(missing_emm_columns) > 0) {
  
  stop(
    paste(
      "The EMM file is missing:",
      paste(
        missing_emm_columns,
        collapse = ", "
      )
    )
  )
}

if (length(missing_contrast_columns) > 0) {
  
  stop(
    paste(
      "The contrast file is missing:",
      paste(
        missing_contrast_columns,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 5. CLEAN TEXT
#
# This removes accidental backslashes before underscores.
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
    )
  )

site_contrasts <- site_contrasts %>%
  mutate(
    contrast = str_replace_all(
      contrast,
      "\\\\_",
      "_"
    ),
    
    contrast = str_squish(contrast),
    
    Index = str_replace_all(
      Index,
      "\\\\_",
      "_"
    )
  )


# ------------------------------------------------------------
# 6. SPLIT EACH CONTRAST INTO ITS TWO SITES
#
# Expected input:
# IARR_A_NorthRK - KKRR_B_Wadi
# ------------------------------------------------------------

site_contrasts <- site_contrasts %>%
  
  separate(
    contrast,
    into = c("Site1", "Site2"),
    sep = "\\s+-\\s+",
    remove = FALSE,
    extra = "merge",
    fill = "right"
  ) %>%
  
  mutate(
    Site1 = str_trim(Site1),
    Site2 = str_trim(Site2),
    
    # multcompLetters expects names such as Site1-Site2
    Comparison = paste0(
      Site1,
      "-",
      Site2
    )
  )


# ------------------------------------------------------------
# 7. CHECK THAT CONTRASTS WERE SPLIT CORRECTLY
# ------------------------------------------------------------

bad_contrasts <- site_contrasts %>%
  filter(
    is.na(Site1) |
      is.na(Site2) |
      Site1 == "" |
      Site2 == ""
  )

if (nrow(bad_contrasts) > 0) {
  
  print(bad_contrasts)
  
  stop(
    paste0(
      "One or more contrasts could not be separated. ",
      "The expected format is 'Site1 - Site2'."
    )
  )
}


# ------------------------------------------------------------
# 8. FUNCTION TO CALCULATE COMPACT-LETTER DISPLAYS
#
# Letters are calculated independently for every index.
# ------------------------------------------------------------

make_site_letters <- function(dat) {
  
  dat <- dat %>%
    filter(
      !is.na(p.value),
      !is.na(Comparison)
    )
  
  if (nrow(dat) == 0) {
    
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
    SiteID = names(letter_result),
    SiteLetter = unname(letter_result)
  )
}


# ------------------------------------------------------------
# 9. CALCULATE LETTERS FOR EACH INDEX
# ------------------------------------------------------------

site_letters <- site_contrasts %>%
  
  group_by(Index) %>%
  
  group_modify(
    ~ make_site_letters(.x)
  ) %>%
  
  ungroup()


# ------------------------------------------------------------
# 10. CHECK LETTER OUTPUT
# ------------------------------------------------------------

cat(
  "\nNumber of letter records:",
  nrow(site_letters),
  "\n"
)

print(
  site_letters %>%
    arrange(Index, SiteID)
)


# ------------------------------------------------------------
# 11. JOIN LETTERS TO ESTIMATED MARGINAL MEANS
# ------------------------------------------------------------

plot_dat <- emm %>%
  
  left_join(
    site_letters,
    by = c(
      "Index",
      "SiteID"
    )
  )


# ------------------------------------------------------------
# 12. CRITICAL JOIN DIAGNOSTIC
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
  sum(!is.na(plot_dat$SiteLetter)),
  "\n"
)

cat(
  "Number without letters:",
  sum(is.na(plot_dat$SiteLetter)),
  "\n"
)

if (any(is.na(plot_dat$SiteLetter))) {
  
  cat(
    "\nRows without matched letters:\n"
  )
  
  print(
    plot_dat %>%
      filter(is.na(SiteLetter)) %>%
      select(
        Index,
        SiteID,
        emmean
      )
  )
}


# ------------------------------------------------------------
# 13. EXCLUDE SOUND-LEVEL DIAGNOSTIC VARIABLES
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
# 14. CLEAN SITE LABELS
# ------------------------------------------------------------

site_labels <- c(
  "IARR_A_NorthRK"   = "IA1",
  "IARR_A_SouthRK"   = "IA2",
  "KKRR_B_Falaq"     = "KK1",
  "KKRR_A_Catchment" = "KK2",
  "KKRR_B_Wadi"      = "KK3"
)

site_order <- c(
  "IA1",
  "IA2",
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
# 15. CLEAN INDEX LABELS
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
  "mean_A"       = "A",
  "mean_B"       = "B"
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
  "A",
  "B"
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
    )
  )


# ------------------------------------------------------------
# 16. CREATE RESERVE VARIABLE
# ------------------------------------------------------------

plot_dat <- plot_dat %>%
  mutate(
    Reserve = case_when(
      str_starts(SiteID, "IARR") ~ "IARR",
      str_starts(SiteID, "KKRR") ~ "KKRR",
      TRUE ~ "Other"
    )
  )


# ------------------------------------------------------------
# 17. POSITION LETTERS ABOVE CONFIDENCE INTERVALS
#
# The offset is calculated independently for each index
# because the response scales differ greatly.
# ------------------------------------------------------------

plot_dat <- plot_dat %>%
  
  group_by(IndexLabel) %>%
  
  mutate(
    panel_min = min(
      lower.CL,
      na.rm = TRUE
    ),
    
    panel_max = max(
      upper.CL,
      na.rm = TRUE
    ),
    
    panel_range = panel_max - panel_min,
    
    panel_range = if_else(
      is.finite(panel_range) &
        panel_range > 0,
      panel_range,
      1
    ),
    
    letter_y = upper.CL +
      0.1 * panel_range
  ) %>%
  
  ungroup()


# ------------------------------------------------------------
# 18. COLOURS AND SYMBOLS
#
# Colour identifies reserve.
# Shape identifies the individual site.
# ------------------------------------------------------------

reserve_colours <- c(
  "IARR" = "#2166AC",
  "KKRR" = "#D95F02",
  "Other" = "#585858"
)

site_shapes <- c(
  "IA1"     = 16,
  "IA2"     = 17,
  "KK1"     = 15,
  "KK2"     = 18,
  "KK3"      = 25
)


# ------------------------------------------------------------
# 19. CREATE SPATIAL HETEROGENEITY FIGURE
# ------------------------------------------------------------

p_spatial <- ggplot(
  plot_dat,
  aes(
    x = Site,
    y = emmean,
    colour = Reserve,
    shape = Site
  )
) +
  
  # Confidence intervals
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.12,
    linewidth = 0.7
  ) +
  
  # Estimated marginal means
  geom_point(
    size = 4,
    stroke = 1.1
  ) +
  
  # Compact-letter displays
  # Letters inherit the reserve colour used for their points
  geom_text(
    aes(
      y = letter_y,
      label = SiteLetter
    ),
    size = 6,
    fontface = "bold",
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~ IndexLabel,
    scales = "free_y",
    ncol = 3
  ) +
  
  scale_colour_manual(
    values = reserve_colours,
    name = "Reserve"
  ) +
  
  scale_shape_manual(
    values = site_shapes,
    name = "Site"
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.17)
    )
  ) +
  
  labs(
    x = NULL,
    y = "Estimated marginal mean"
  ) +
  
  guides(
    colour = guide_legend(
      order = 1,
      nrow = 1
    ),
    
    shape = guide_legend(
      order = 2,
      nrow = 1,
      byrow = TRUE
    )
  ) +
  
  theme_bw(
    base_size = 15
  ) +
  
  theme(
    strip.text = element_text(
      size = 15,
      face = "bold"
    ),
    
    strip.background = element_rect(
      fill = "grey95",
      colour = "black",
      linewidth = 0.6
    ),
    
    axis.title.y = element_text(
      size = 18
    ),
    
    axis.text.y = element_text(
      size = 18,
      colour = "black"
    ),
    
    axis.text.x = element_text(
      size = 15,
      colour = "black",
      angle = 90,
      hjust = 1,
      vjust = 1
    ),
    
    legend.position = "bottom",
    
    legend.box = "vertical",
    
    legend.title = element_text(
      size = 15,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 15
    ),
    
    panel.spacing = grid::unit(
      1,
      "lines"
    ),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )


# ------------------------------------------------------------
# 20. DISPLAY FIGURE
# ------------------------------------------------------------

print(p_spatial)


# ------------------------------------------------------------
# 21. SAVE HIGH-RESOLUTION PNG
# ------------------------------------------------------------

ggsave(
  filename = output_png,
  plot = p_spatial,
  width = 14,
  height = 14,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 22. SAVE VECTOR PDF
# ------------------------------------------------------------

ggsave(
  filename = output_pdf,
  plot = p_spatial,
  width = 14,
  height = 14,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)