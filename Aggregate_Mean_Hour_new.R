library(dplyr)

# ============================================================
# 1. IMPORT MINUTE-LEVEL INDEX DATA
# ============================================================

df <- read.csv(
  "D:/Flor/SaudiArabia/Data/Results/Indices/SaudiArabia_Allindices_Oct2026_multiplefreqBins.csv",
  stringsAsFactors = FALSE
)

# ============================================================
# 2. FORMAT VARIABLES
# ============================================================

df <- df %>%
  mutate(
    RecorderID = as.character(RecorderID),
    SiteID = as.character(SiteID),
    Reserve = as.character(Reserve),
    
    Date = as.character(Date),
    
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = as.integer(Hour)
  )

# ============================================================
# 3. CALCULATE HOURLY MEANS
# ============================================================

df_hourly_means <- df %>%
  
  group_by(
    RecorderID,
    SiteID,
    Reserve,
    Date,
    Year,
    Month,
    Day,
    Hour
  ) %>%
  
  summarise(
    
    mean_ACI300 = mean(
      ACI300,
      na.rm = TRUE
    ),
    
    mean_ACI2000 = mean(
      ACI2000,
      na.rm = TRUE
    ),
    
    mean_ADI300 = mean(
      ADI300,
      na.rm = TRUE
    ),
    
    mean_ADI2000 = mean(
      ADI2000,
      na.rm = TRUE
    ),
    
    mean_BIO300 = mean(
      BIO300,
      na.rm = TRUE
    ),
    
    mean_BIO2000 = mean(
      BIO2000,
      na.rm = TRUE
    ),
    
    mean_NP300 = mean(
      NP300,
      na.rm = TRUE
    ),
    
    mean_NP2000 = mean(
      NP2000,
      na.rm = TRUE
    ),
    
    mean_NDSI = mean(
      NDSI,
      na.rm = TRUE
    ),
    
    mean_B = mean(
      B,
      na.rm = TRUE
    ),
    
    mean_A = mean(
      A,
      na.rm = TRUE
    ),
    
    mean_NoiseFloor_dB = mean(
      NoiseFloor_dB,
      na.rm = TRUE
    ),
    
    mean_SignalLevel_dB = mean(
      SignalLevel_dB,
      na.rm = TRUE
    ),
    
    mean_SNR_dB = mean(
      SNR_dB,
      na.rm = TRUE
    ),
    
    # Number of usable 1-minute recordings
    # contributing to this hourly value
    n_records = n(),
    
    .groups = "drop"
  )

# ============================================================
# 4. CHECK NUMBER OF RECORDINGS PER HOUR
# ============================================================

print(
  table(
    df_hourly_means$n_records
  )
)

# ============================================================
# 5. SAVE HOURLY DATA
# ============================================================

write.csv(
  df_hourly_means,
  "D:/Flor/SaudiArabia/Data/Results/Indices/SaudiArabia_Allindices_Oct2026_multiplefreqBins_HOURLY.csv",
  row.names = FALSE
)

# ============================================================
# 6. DISPLAY
# ============================================================

print(
  head(df_hourly_means)
)

cat(
  "\nHourly rows created:",
  nrow(df_hourly_means),
  "\n"
)