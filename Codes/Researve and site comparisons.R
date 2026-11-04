# ==============================================================================
# Hyper-arid soundscapes: an exploratory study of bioacoustics indices
# in Saudi Arabia
# ==============================================================================
#
# AUTHORS
# Florencia Sangermano¹
# Mariska Weijerman¹*
# Aparajita Tyagi¹
# Rab Nawaz²
# Abdullatif Al-Abdulwahab²
#
# AFFILIATIONS
# ¹ Graduate School of Geography, Clark University,
#   Worcester, MA 01610, USA
#
# ² Imam Abdulaziz bin Mohammed Royal Reserve Development Authority,
#   Riyadh, 13531, Saudi Arabia
#
# CODE AUTHORS
# Florencia Sangermano
# Aparajita Tyagi
#
# AI USE STATEMENT
# Generative artificial intelligence (AI) tools were used as programming aids
# during code development, including assistance with code syntax, debugging,
# and code organization. Code was reviewed by the authors to ensure accuracy.
# All analytical decisions, model specifications,implementation, verification, and interpretation of results were conducted
# by the authors.
#
# ==============================================================================







# ============================================================
# MODELING SECTION
#
# All models use:
# Dawn, Morning, and Dusk only
#
# Questions:
# 1. Do reserves differ overall across the retained periods?
# 2. Do sites differ overall across the retained periods?
# 3. Do periods differ within each site?
# 4. Do sites differ within each period? 
# ============================================================

# ============================================================
# PACKAGES
# ============================================================

library(mgcv)
library(dplyr)
library(emmeans)
library(suncalc)

# ============================================================
# OUTPUT DIRECTORIES
# ============================================================

results_directory <- paste0(
  "D:/Flor/SaudiArabia/Data/Results/",
  "SiteAndReserveComparison_outputs1"
)

dir.create(
  results_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. IMPORT DATA
# ============================================================

df <- read.csv(
  paste0(
    "D:/Flor/SaudiArabia/Data/Results/Indices/",
    "SaudiArabia_Allindices_Oct2026_multiplefreqBins_HOURLY.csv"
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 2. KEEP MARCH AND APRIL
# ============================================================

df <- df %>%
  filter(
    Month %in% c(3, 4)
  )

# ============================================================
# 3. FORMAT VARIABLES
# ============================================================

df <- df %>%
  mutate(
    SiteID = factor(SiteID),
    Reserve = factor(Reserve),
    
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = as.integer(Hour),
    
    Date = as.Date(
      sprintf(
        "%04d-%02d-%02d",
        Year,
        Month,
        Day
      )
    ),
    
    DateTime = as.POSIXct(
      sprintf(
        "%04d-%02d-%02d %02d:00:00",
        Year,
        Month,
        Day,
        Hour
      ),
      tz = "Asia/Riyadh"
    )
  )

df$DateTimeNum <- as.numeric(
  df$DateTime
)

# ============================================================
# 4. CALCULATE ASTRONOMICAL TIMES
# ============================================================

# Riyadh coordinates
lat <- 24.7136
lon <- 46.6753

dates <- data.frame(
  date = unique(df$Date),
  lat = lat,
  lon = lon
)

sun <- getSunlightTimes(
  data = dates,
  keep = c(
    "sunrise",
    "sunset",
    "solarNoon"
  ),
  tz = "Asia/Riyadh"
)

# ============================================================
# 5. CONVERT ASTRONOMICAL TIMES TO DECIMAL HOURS
# ============================================================

to_decimal_hour <- function(x){
  
  as.numeric(format(x, "%H")) +
    as.numeric(format(x, "%M")) / 60 +
    as.numeric(format(x, "%S")) / 3600
}

sun <- sun %>%
  mutate(
    SunriseHour = to_decimal_hour(sunrise),
    SunsetHour = to_decimal_hour(sunset),
    SolarNoonHour = to_decimal_hour(solarNoon)
  ) %>%
  dplyr::select(
    Date = date,
    SunriseHour,
    SunsetHour,
    SolarNoonHour
  )

# ============================================================
# 6. JOIN ASTRONOMICAL TIMES
# ============================================================

df <- left_join(
  df,
  sun,
  by = "Date"
)

if(any(is.na(df$SunriseHour))){
  
  stop(
    paste(
      sum(is.na(df$SunriseHour)),
      "observations did not receive astronomical times."
    )
  )
}

# Convert Date to factor for random-effect smoother
df$Date <- factor(df$Date)

# ============================================================
# 7. DEFINE ALL SIX ASTRONOMICAL PERIODS
# ============================================================

df <- df %>%
  mutate(
    
    Period = case_when(
      
      # Dawn: 2 hours before to 2 hours after sunrise
      Hour >= SunriseHour - 2 &
        Hour <= SunriseHour + 2 ~
        "Dawn",
      
      # Morning: after dawn to 2 hours before solar noon
      Hour > SunriseHour + 2 &
        Hour < SolarNoonHour - 2 ~
        "Morning",
      
      # Midday: 2 hours before to 2 hours after solar noon
      Hour >= SolarNoonHour - 2 &
        Hour <= SolarNoonHour + 2 ~
        "Midday",
      
      # Afternoon: after midday to 2 hours before sunset
      Hour > SolarNoonHour + 2 &
        Hour < SunsetHour - 2 ~
        "Afternoon",
      
      # Dusk: 2 hours before to 2 hours after sunset
      Hour >= SunsetHour - 2 &
        Hour <= SunsetHour + 2 ~
        "Dusk",
      
      # All remaining hours
      TRUE ~
        "Night"
    ),
    
    Period = factor(
      Period,
      levels = c(
        "Night",
        "Dawn",
        "Morning",
        "Midday",
        "Afternoon",
        "Dusk"
      )
    )
  )

# Check assignment
print(
  table(
    df$Reserve,
    df$Period
  )
)

# ============================================================
# 8. SORT FULL DATASET
# ============================================================

df <- df %>%
  arrange(
    SiteID,
    DateTime
  )


# ============================================================
# 9. CREATE RESTRICTED DIEL DATASET
# Use Dawn, Morning, and Dusk only
# ============================================================

df_diel <- df %>%
  filter(
    Period %in% c(
      "Dawn",
      "Morning",
      "Dusk"
    )
  ) %>%
  droplevels() %>%
  arrange(
    SiteID,
    DateTime
  ) %>%
  group_by(SiteID) %>%
  mutate(
    AR.start = c(
      TRUE,
      diff(DateTimeNum) != 3600
    )
  ) %>%
  ungroup()

# Set desired display and comparison order
df_diel$Period <- factor(
  df_diel$Period,
  levels = c(
    "Dawn",
    "Morning",
    "Dusk"
  )
)

cat(
  "\nObservations retained for Models 2 and 3:",
  nrow(df_diel),
  "\n"
)

print(
  table(
    df_diel$Reserve,
    df_diel$Period
  )
)

# ============================================================
# 11. INDICES
# ============================================================

indices <- c(
  "mean_ACI300",
  "mean_ACI2000",
  "mean_ADI300",
  "mean_ADI2000",
  "mean_BIO300",
  "mean_BIO2000",
  "mean_NP300",
  "mean_NP2000",
  "mean_NDSI",
  "mean_B",
  "mean_A"
)

# Confirm that all response variables exist
missing_indices <- setdiff(
  indices,
  names(df)
)

if(length(missing_indices) > 0){
  
  stop(
    paste0(
      "These index columns are missing from df:\n",
      paste(
        missing_indices,
        collapse = ", "
      )
    )
  )
}

# ============================================================
# 12. FUNCTION TO ESTIMATE INDEX-SPECIFIC RHO
# note: code includes an" overall model_type not used in paper"
# ============================================================

estimate_rho <- function(
    response_var,
    analysis_data,
    model_type = c(
      "overall",
      "diel"
    )
){
  
  model_type <- match.arg(
    model_type
  )
  
  if(model_type == "overall"){
    
    formula_text <- paste0(
      response_var,
      " ~ Reserve +
         s(Hour, k = 24) +
         s(Date, bs = 're') +
         s(SiteID, bs = 're')"
    )
    
  } else {
    
    formula_text <- paste0(
      response_var,
      " ~ Reserve * Period +
         s(Date, bs = 're') +
         s(SiteID, bs = 're')"
    )
  }
  
  temp_model <- bam(
    formula = as.formula(
      formula_text
    ),
    data = analysis_data,
    method = "fREML",
    discrete = TRUE
  )
  
  rho_value <- acf(
    residuals(temp_model),
    lag.max = 1,
    plot = FALSE,
    na.action = na.pass
  )$acf[2]
  
  as.numeric(
    rho_value
  )
}


# ============================================================
# 13. ESTIMATE RHO
# Dawn, Morning, and Dusk only
# ============================================================

cat(
  "\nEstimating rho values for restricted diel models\n"
)

RhoDielSummary <- data.frame(
  Index = indices,
  Rho = sapply(
    indices,
    estimate_rho,
    analysis_data = df_diel,
    model_type = "diel"
  ),
  row.names = NULL,
  stringsAsFactors = FALSE
)

print(
  RhoDielSummary
)

write.csv(
  RhoDielSummary,
  file.path(
    results_directory,
    "RhoValues_DawnMorningDusk_Mean.csv"
  ),
  row.names = FALSE
)

# ============================================================
# FUNCTION TO CREATE INDEX-SPECIFIC MODEL DATA
# ============================================================
#
# This function:
# 1. Removes missing or non-finite values for the current index.
# 2. Sorts observations chronologically within each site.
# 3. Recalculates AR.start after removing missing observations.
# 4. Drops unused factor levels.
#
# This prevents:
# "variable lengths differ (found for AR.start)"
# ============================================================

prepare_model_data <- function(
    analysis_data,
    response_var
){

  model_data <- analysis_data %>%

    filter(
      !is.na(.data[[response_var]]),
      is.finite(.data[[response_var]])
    ) %>%

    arrange(
      SiteID,
      DateTime
    ) %>%

    group_by(
      SiteID
    ) %>%

    mutate(

      AR.start = c(
        TRUE,
        diff(DateTimeNum) != 3600
      )
    ) %>%

    ungroup() %>%

    droplevels()

  return(
    model_data
  )
}


# ============================================================
# 14. CREATE EMPTY OUTPUT TABLES
# ============================================================

OverallReserveSummary <- data.frame()

OverallSiteSummary <- data.frame()

OverallSiteEstimatedMeans <- data.frame()

PeriodWithinSiteSummary <- data.frame()

SiteWithinPeriodSummary <- data.frame()

SitePeriodEstimatedMeans <- data.frame()


# ============================================================
#15.  MODEL 1
# OVERALL RESERVE DIFFERENCES
#
# Question 1:
# Do IARR and KKRR differ overall when considering only
# Dawn, Morning, and Dusk and adjusting for diel period?
#
# Structure:
# Reserve + Period
# ============================================================

for(i in indices){
  
  cat(
    "\nRunning Model 1:",
    i,
    "\n"
  )
  
  # ----------------------------------------------------------
  # Retrieve index-specific diel rho
  # ----------------------------------------------------------
  
  rho_index <- RhoDielSummary$Rho[
    RhoDielSummary$Index == i
  ]
  
  if(length(rho_index) != 1 ||
     !is.finite(rho_index)){
    
    stop(
      paste(
        "No valid diel rho value was found for",
        i
      )
    )
  }
  
  # ----------------------------------------------------------
  # Create index-specific data and AR.start
  # ----------------------------------------------------------
  
  model_data <- prepare_model_data(
    analysis_data = df_diel,
    response_var = i
  )
  
  cat(
    "Rows used:",
    nrow(model_data),
    "| rho:",
    round(rho_index, 3),
    "\n"
  )
  
  # ----------------------------------------------------------
  # Fit Model 1
  # ----------------------------------------------------------
  
  model1 <- bam(
    
    as.formula(
      paste0(
        i,
        " ~ Reserve +
           Period +
           s(Date, bs = 're') +
           s(SiteID, bs = 're')"
      )
    ),
    
    data = model_data,
    
    method = "fREML",
    
    discrete = TRUE,
    
    rho = rho_index,
    
    AR.start = model_data$AR.start,
    
    na.action = na.exclude
  )
  
  # ----------------------------------------------------------
  # Overall reserve estimated means
  # ----------------------------------------------------------
  
  emm_reserve_overall <- emmeans(
    model1,
    ~ Reserve,
    weights = "equal"
  )
  
  reserve_means <- as.data.frame(
    emm_reserve_overall
  )
  
  reserve_contrast <- as.data.frame(
    pairs(
      emm_reserve_overall
    )
  )
  
  OverallReserveSummary <- bind_rows(
    
    OverallReserveSummary,
    
    data.frame(
      Index = i,
      
      IARR = reserve_means$emmean[
        reserve_means$Reserve == "IARR"
      ],
      
      KKRR = reserve_means$emmean[
        reserve_means$Reserve == "KKRR"
      ],
      
      Difference =
        reserve_contrast$estimate,
      
      SE =
        reserve_contrast$SE,
      
      df =
        reserve_contrast$df,
      
      t_value =
        reserve_contrast$t.ratio,
      
      p_value =
        reserve_contrast$p.value,
      
      N =
        nrow(model_data),
      
      Rho =
        rho_index,
      
      stringsAsFactors = FALSE
    )
  )
}


# ============================================================
# SAVE MODEL 1
# ============================================================

write.csv(
  OverallReserveSummary,
  file.path(
    results_directory,
    "Model1_OverallReserveDifferences_DawnMorningDusk.csv"
  ),
  row.names = FALSE
)



# ============================================================
# 16. MODEL 2
# OVERALL SITE DIFFERENCES
#
# Question 2:
# Do monitored sites differ overall when considering only
# Dawn, Morning, and Dusk and adjusting for period?
#
# Structure:
# SiteID + Period
# ============================================================

for(i in indices){
  
  cat(
    "\nRunning Model 4:",
    i,
    "\n"
  )
  
  # ----------------------------------------------------------
  # Retrieve index-specific diel rho
  # ----------------------------------------------------------
  
  rho_index <- RhoDielSummary$Rho[
    RhoDielSummary$Index == i
  ]
  
  if(length(rho_index) != 1 ||
     !is.finite(rho_index)){
    
    stop(
      paste(
        "No valid diel rho value was found for",
        i
      )
    )
  }
  
  # ----------------------------------------------------------
  # Create index-specific model data
  # ----------------------------------------------------------
  
  model_data <- prepare_model_data(
    analysis_data = df_diel,
    response_var = i
  )
  
  cat(
    "Rows used:",
    nrow(model_data),
    "| rho:",
    round(rho_index, 3),
    "\n"
  )
  
  # ----------------------------------------------------------
  # Fit overall site model
  # ----------------------------------------------------------
  
  site_model <- bam(
    
    as.formula(
      paste0(
        i,
        " ~ SiteID +
           Period +
           s(Date, bs = 're')"
      )
    ),
    
    data = model_data,
    
    method = "fREML",
    
    discrete = TRUE,
    
    rho = rho_index,
    
    AR.start = model_data$AR.start,
    
    na.action = na.exclude
  )
  
  # ----------------------------------------------------------
  # Estimated means by site, averaged equally over periods
  # ----------------------------------------------------------
  
  emm_site <- emmeans(
    site_model,
    ~ SiteID,
    weights = "equal"
  )
  
  site_means <- as.data.frame(
    emm_site
  )
  
  site_means$Index <- i
  site_means$N <- nrow(model_data)
  site_means$Rho <- rho_index
  
  OverallSiteEstimatedMeans <- bind_rows(
    OverallSiteEstimatedMeans,
    site_means
  )
  
  # ----------------------------------------------------------
  # Pairwise site comparisons
  # ----------------------------------------------------------
  
  site_contrasts <- as.data.frame(
    pairs(
      emm_site,
      adjust = "holm"
    )
  )
  
  site_contrasts$Index <- i
  site_contrasts$N <- nrow(model_data)
  site_contrasts$Rho <- rho_index
  
  OverallSiteSummary <- bind_rows(
    OverallSiteSummary,
    site_contrasts
  )
}


# ============================================================
# ADD SIGNIFICANCE LABELS TO MODEL 2
# ============================================================

OverallSiteSummary <- OverallSiteSummary %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

SignificantOverallSiteComparisons <-
  OverallSiteSummary %>%
  filter(
    p.value < 0.05
  )


# ============================================================
# SAVE MODEL2
# ============================================================

write.csv(
  OverallSiteEstimatedMeans,
  file.path(
    results_directory,
    "Model2_OverallSiteMeans_DawnMorningDusk.csv"
  ),
  row.names = FALSE
)

write.csv(
  OverallSiteSummary,
  file.path(
    results_directory,
    "Model2_OverallSiteComparisons_DawnMorningDusk_Holm.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. MODEL 3
# SITE-BY-PERIOD INTERACTION
#
# Question 3:
# Within each monitored site, do indices differ among
# Dawn, Morning, and Dusk?
#
# Question 4 :
# Which sites differ within Dawn, Morning, and Dusk?
#
# Structure:
# SiteID * Period
# ============================================================

for(i in indices){
  
  cat(
    "\nRunning Model 5:",
    i,
    "\n"
  )
  
  # ----------------------------------------------------------
  # Retrieve index-specific diel rho
  # ----------------------------------------------------------
  
  rho_index <- RhoDielSummary$Rho[
    RhoDielSummary$Index == i
  ]
  
  if(length(rho_index) != 1 ||
     !is.finite(rho_index)){
    
    stop(
      paste(
        "No valid diel rho value was found for",
        i
      )
    )
  }
  
  # ----------------------------------------------------------
  # Create index-specific model data
  # ----------------------------------------------------------
  
  model_data <- prepare_model_data(
    analysis_data = df_diel,
    response_var = i
  )
  
  cat(
    "Rows used:",
    nrow(model_data),
    "| rho:",
    round(rho_index, 3),
    "\n"
  )
  
  # ----------------------------------------------------------
  # Fit Site x Period model
  # ----------------------------------------------------------
  
  site_period_model <- bam(
    
    as.formula(
      paste0(
        i,
        " ~ SiteID * Period +
           s(Date, bs = 're')"
      )
    ),
    
    data = model_data,
    
    method = "fREML",
    
    discrete = TRUE,
    
    rho = rho_index,
    
    AR.start = model_data$AR.start,
    
    na.action = na.exclude
  )
  
  # ==========================================================
  # ESTIMATED MEANS FOR EACH SITE-PERIOD COMBINATION
  # ==========================================================
  
  emm_period_site <- emmeans(
    site_period_model,
    ~ Period | SiteID
  )
  
  site_period_means <- as.data.frame(
    emm_period_site
  )
  
  site_period_means$Index <- i
  site_period_means$N <- nrow(model_data)
  site_period_means$Rho <- rho_index
  
  SitePeriodEstimatedMeans <- bind_rows(
    SitePeriodEstimatedMeans,
    site_period_means
  )
  
  # ==========================================================
  # PERIOD DIFFERENCES WITHIN EACH SITE
  #
  # Directly answers Question 3
  # ==========================================================
  
  period_within_site <- as.data.frame(
    pairs(
      emm_period_site,
      adjust = "holm"
    )
  )
  
  period_within_site$Index <- i
  period_within_site$N <- nrow(model_data)
  period_within_site$Rho <- rho_index
  
  PeriodWithinSiteSummary <- bind_rows(
    PeriodWithinSiteSummary,
    period_within_site
  )
  
  # ==========================================================
  # SITE DIFFERENCES WITHIN EACH PERIOD
  #
  # Question 4:  Site x Period interpretation
  # ==========================================================
  
  emm_site_period <- emmeans(
    site_period_model,
    ~ SiteID | Period
  )
  
  site_within_period <- as.data.frame(
    pairs(
      emm_site_period,
      adjust = "holm"
    )
  )
  
  site_within_period$Index <- i
  site_within_period$N <- nrow(model_data)
  site_within_period$Rho <- rho_index
  
  SiteWithinPeriodSummary <- bind_rows(
    SiteWithinPeriodSummary,
    site_within_period
  )
}


# ============================================================
# ADD SIGNIFICANCE LABELS TO MODEL 3
# ============================================================

PeriodWithinSiteSummary <- PeriodWithinSiteSummary %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

SiteWithinPeriodSummary <- SiteWithinPeriodSummary %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )


# ============================================================
# CREATE SIGNIFICANT MODEL 3 TABLES
# ============================================================

SignificantPeriodsWithinSite <-
  PeriodWithinSiteSummary %>%
  filter(
    p.value < 0.05
  )

SignificantSitesWithinPeriod <-
  SiteWithinPeriodSummary %>%
  filter(
    p.value < 0.05
  )


# ============================================================
# SAVE MODEL 3
# ============================================================

write.csv(
  SitePeriodEstimatedMeans,
  file.path(
    results_directory,
    "Model3_SiteByPeriod_EstimatedMeans_DawnMorningDusk.csv"
  ),
  row.names = FALSE
)

write.csv(
  PeriodWithinSiteSummary,
  file.path(
    results_directory,
    "Model3_PeriodDifferencesWithinSite_DawnMorningDusk_Holm.csv"
  ),
  row.names = FALSE
)



write.csv(
  SiteWithinPeriodSummary,
  file.path(
    results_directory,
    "Model3_SiteDifferencesWithinPeriod_DawnMorningDusk_Holm.csv"
  ),
  row.names = FALSE
)

