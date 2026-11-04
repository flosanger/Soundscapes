# ============================================================
# MODEL 5
# SITE x PERIOD ESTIMATED MARGINAL MEANS
#
# Lowercase letters compare PERIODS WITHIN EACH SITE.
# Sites sharing a lowercase letter across periods are not
# significantly different based on Holm-adjusted comparisons.
# ============================================================

library(tidyverse)
library(multcompView)

# ------------------------------------------------------------
# 1. FILE PATHS
# ------------------------------------------------------------

results_dir <- paste0(
  "D:/Flor/SaudiArabia/Data/Results/",
  "SiteAndReserveComparison/"
)

emm_file <- paste0(
  results_dir,
  "Model5_SiteByPeriod_EstimatedMeans_DawnMorningDusk.csv"
)

# This must contain pairwise PERIOD comparisons within each site
period_contrast_file <- paste0(
  results_dir,
  "Model5_PeriodDifferencesWithinSite_DawnMorningDusk_Holm.csv"
)

output_figure <- paste0(
  results_dir,
  "Figure3-Period Differences within site- Model5.png"
)


# ------------------------------------------------------------
# 2. READ FILES
# ------------------------------------------------------------

emm <- read.csv(
  emm_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

period_contrasts <- read.csv(
  period_contrast_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. CHECK REQUIRED COLUMNS
# ------------------------------------------------------------

required_emm_columns <- c(
  "Period", "SiteID", "emmean",
  "lower.CL", "upper.CL", "Index"
)

required_contrast_columns <- c(
  "contrast", "SiteID", "p.value", "Index"
)

missing_emm <- setdiff(required_emm_columns, names(emm))

missing_contrasts <- setdiff(
  required_contrast_columns,
  names(period_contrasts)
)

if (length(missing_emm) > 0) {
  stop(
    paste(
      "The EMM file is missing:",
      paste(missing_emm, collapse = ", ")
    )
  )
}

if (length(missing_contrasts) > 0) {
  stop(
    paste(
      "The period contrast file is missing:",
      paste(missing_contrasts, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 4. STANDARDISE TEXT
# ------------------------------------------------------------

# Remove accidental backslashes before underscores if present
emm <- emm %>%
  mutate(
    SiteID = str_replace_all(SiteID, "\\\\_", "_"),
    Index  = str_replace_all(Index, "\\\\_", "_")
  )

period_contrasts <- period_contrasts %>%
  mutate(
    SiteID  = str_replace_all(SiteID, "\\\\_", "_"),
    Index   = str_replace_all(Index, "\\\\_", "_"),
    contrast = str_replace_all(contrast, "\\\\_", "_"),
    contrast = str_squish(contrast)
  )


# ------------------------------------------------------------
# 5. ORDER PERIODS
# ------------------------------------------------------------

period_levels <- c("Dawn", "Morning", "Dusk")

emm <- emm %>%
  mutate(
    Period = factor(
      Period,
      levels = period_levels
    )
  )


# ------------------------------------------------------------
# 6. CALCULATE COMPACT-LETTER DISPLAYS
#
# Letters are calculated separately for every:
# Index x SiteID combination.
#
# Therefore, the letters compare Dawn, Morning and Dusk
# only within the same site and acoustic index.
# ------------------------------------------------------------

# ------------------------------------------------------------
# CALCULATE LOWERCASE LETTERS FOR PERIOD DIFFERENCES
# WITHIN EACH SITE
# ------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(multcompView)

period_levels <- c("Dawn", "Morning", "Dusk")

# Clean contrast names
period_contrasts2 <- period_contrasts %>%
  mutate(
    SiteID = str_replace_all(SiteID, "\\\\_", "_"),
    Index  = str_replace_all(Index, "\\\\_", "_"),
    contrast = str_squish(contrast)
  ) %>%
  
  # Split "Dawn - Morning" into the two periods
  separate(
    contrast,
    into = c("Period1", "Period2"),
    sep = "\\s+-\\s+",
    remove = FALSE
  ) %>%
  
  mutate(
    Period1 = str_trim(Period1),
    Period2 = str_trim(Period2),
    
    # multcompLetters wants comparison names such as Dawn-Morning
    Comparison = paste0(Period1, "-", Period2)
  )


# ------------------------------------------------------------
# FUNCTION TO CONVERT PAIRWISE P-VALUES TO LETTERS
# ------------------------------------------------------------

get_letters <- function(df) {
  
  pvals <- df$p.value
  names(pvals) <- df$Comparison
  
  letters_out <- multcompView::multcompLetters(
    pvals,
    threshold = 0.05
  )$Letters
  
  tibble(
    Period = names(letters_out),
    PeriodLetter = unname(letters_out)
  )
}


# ------------------------------------------------------------
# CREATE LETTERS SEPARATELY FOR EACH INDEX x SITE
# ------------------------------------------------------------

period_letters <- period_contrasts2 %>%
  group_by(Index, SiteID) %>%
  group_modify(~ get_letters(.x)) %>%
  ungroup()


# Force period names to match EMM file exactly
period_letters <- period_letters %>%
  mutate(
    Period = str_trim(Period),
    Period = factor(
      Period,
      levels = period_levels
    )
  )


# ------------------------------------------------------------
# CHECK LETTERS BEFORE JOINING
# ------------------------------------------------------------

print(period_letters)

print(
  table(
    period_letters$Period,
    useNA = "ifany"
  )
)

print(
  table(
    period_letters$PeriodLetter,
    useNA = "ifany"
  )
)

# ------------------------------------------------------------
# 7. JOIN LETTERS TO EMM DATA
# ------------------------------------------------------------

plot_dat <- emm %>%
  left_join(
    period_letters,
    by = c("Index", "SiteID", "Period")
  )


# ------------------------------------------------------------
# 8. EXCLUDE SOUND-LEVEL DIAGNOSTIC VARIABLES
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
# 9. CLEAN AND ORDER SITE LABELS
# ------------------------------------------------------------

site_labels <- c(
  "IARR_A_NorthRK"   = "IA1",
  "IARR_A_SouthRK"   = "IA2",
  "KKRR_B_Falaq"     = "KK1",
  "KKRR_A_Catchment" = "KK2",
  "KKRR_B_Wadi"      = "KK3"
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
      levels = c(
        "IA1",
        "IA2",
        "KK1",
        "KK2",
        "KK3"
      )
    )
  )



# ------------------------------------------------------------
# 10. CLEAN INDEX LABELS
# ------------------------------------------------------------

index_labels <- c(
  "mean_ACI300"  = "ACI 300 - 6000 Hz",
  "mean_ACI2000" = "ACI 2000 - 6000 Hz",
  "mean_ADI300"  = "ADI 300 - 6000 Hz",
  "mean_ADI2000" = "ADI 2000 - 6000 Hz",
  "mean_BIO300"  = "BIO 300 - 6000 Hz",
  "mean_BIO2000" = "BIO 2000 - 6000 Hz",
  "mean_NP300"   = "NP 300 - 6000 Hz",
  "mean_NP2000"  = "NP 2000 - 6000 Hz",
  "mean_NDSI"    = "NDSI",
  "mean_A"       = "A",
  "mean_B"       = "B"
)

index_order <- c(
  "ACI 300 - 6000 Hz",
  "ACI 2000 - 6000 Hz",
  "ADI 300 - 6000 Hz",
  "ADI 2000 - 6000 Hz",
  "BIO 300 - 6000 Hz",
  "BIO 2000 - 6000 Hz",
  "NP 300 - 6000 Hz",
  "NP 2000 - 6000 Hz",
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
# 11. SET LETTER POSITIONS
#
# Because each index has a different scale, calculate the
# vertical offset separately for each index.
#
# The letters are placed above the upper 95% confidence limit.
# ------------------------------------------------------------

plot_dat <- plot_dat %>%
  group_by(IndexLabel) %>%
  mutate(
    panel_range = max(upper.CL, na.rm = TRUE) -
      min(lower.CL, na.rm = TRUE),
    
    panel_range = if_else(
      is.finite(panel_range) & panel_range > 0,
      panel_range,
      1
    ),
    
    letter_y = upper.CL + 0.055 * panel_range
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 12. RESERVE COLOURS AND SITE SYMBOLS
#
# Colour identifies reserve.
# Shape and line type identify site.
# ------------------------------------------------------------

reserve_colours <- c(
  "IA1" = "#2166AC",
  "IA2" = "#2166AC",
  "KK1" = "#D95F02",
  "KK2" = "#D95F02",
  "KK3" = "#D95F02"
)

site_shapes <- c(
  "IA1" = 16,
  "IA2" = 17,
  "KK1" = 18,
  "KK2" = 15,
  "KK3" = 25
)

site_linetypes <- c(
  "IA1" = "solid",
  "IA2" = "dashed",
  "KK1" = "dashed",
  "KK2" = "solid",
  "KK3" = "dotted"
)


# ------------------------------------------------------------
# 13. CREATE FIGURE
# ------------------------------------------------------------

# ------------------------------------------------------------
# POSITION DODGE
# ------------------------------------------------------------

# Small horizontal separation among sites within each period
pd <- position_dodge(width = 0.5)


# ------------------------------------------------------------
# CREATE FIGURE
# ------------------------------------------------------------

p_emm <- ggplot(
  plot_dat,
  aes(
    x = Period,
    y = emmean,
    group = Site,
    colour = Site,
    shape = Site,
    linetype = Site
  )
) +
  
  # Lines
  geom_line(
    linewidth = 0.9,
    position = pd
  ) +
  
  # 95% confidence intervals
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.05,
    linewidth = 0.55,
    position = pd
  ) +
  
  # EMM points
  geom_point(
    size = 3.1,
    stroke = 1,
    position = pd
  ) +
  
  # Significance letters
  geom_text(
    aes(
      y = letter_y,
      label = PeriodLetter
    ),
    size = 6,
    fontface = "bold",
    position = pd,
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
    name = "Site"
  ) +
  
  scale_shape_manual(
    values = site_shapes,
    name = "Site"
  ) +
  
  scale_linetype_manual(
    values = site_linetypes,
    name = "Site"
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.18)
    )
  ) +
  
  labs(
    x = "Diel period",
    y = "Estimated marginal mean"
  ) +
  
  theme_bw(base_size = 15) +
  
  theme(
    strip.text = element_text(
      size = 16,
      face = "bold"
    ),
    
    strip.background = element_rect(
      fill = "grey95",
      colour = "black"
    ),
    
    axis.title = element_text(
      size = 16
    ),
    
    axis.text = element_text(
      size = 16,
      colour = "black"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 12
    ),
    
    panel.spacing = grid::unit(1, "lines")
  )


print(p_emm)


# ------------------------------------------------------------
# 14. DISPLAY AND SAVE
# ------------------------------------------------------------

print(p_emm)

ggsave(
  filename = output_figure,
  plot = p_emm,
  width = 14,
  height = 14,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = str_replace(
    output_figure,
    "\\.png$",
    ".pdf"
  ),
  plot = p_emm,
  width = 14,
  height = 14,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


