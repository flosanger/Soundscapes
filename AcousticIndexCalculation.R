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


library(tuneR)
library(seewave)
library(soundecology)
library(tidyverse)
library(dplyr)

#=================================================
# Paths
#=================================================

setwd("D:/Flor/SaudiArabia/Data/WavsChannel60")

weather_file <- "D:/Flor/SaudiArabia/Data/Results/WeatherClassification_NewSep2026.csv"

filtered_file <- "D:/Flor/SaudiArabia/Data/Results/Indices/FilteredWeatherTable_Oct2026.csv" #output from filter

outfile <- "D:/Flor/SaudiArabia/Data/Results/Indices/SaudiArabia_Allindices_Oct2026_multiplefreqBins.csv"

checkpoint_file <- "D:/Flor/SaudiArabia/Data/Results/Indices/AcousticIndices_checkpoint_Oct2026_multiplefreqBins.RData"

#=================================================
# Read weather table
#=================================================

WeatherTable <- read.csv(
  weather_file,
  stringsAsFactors = FALSE
)


#=================================================
# 1. PRE-FILTERING
#=================================================


Filtered <- WeatherTable %>%
  filter(WeatherClasses == "NONE" | MaxConfidence <= 0.2203)

cat(nrow(Filtered), "channels kept after weather filtering\n")

#=================================================
# 2. CHANNEL SELECTION PER RECORDING MOMENT
#=================================================
# Group by Site, Date, Hour, Min
# If both channels survive → keep Left (0)
# If only one survives → keep that one

Filtered <- Filtered %>%
  group_by(Site, Date, Hour, Min) %>%
  summarise(
    ChannelUsed = ifelse(
      any(Channel == 0),   # left survives
      0,
      Channel[1]           # only one channel survived
    ),
    File = Filename[Channel == ChannelUsed][1],
    WeatherClasses = WeatherClasses[Channel == ChannelUsed][1],
    MaxConfidence = MaxConfidence[Channel == ChannelUsed][1],
    .groups = "drop"
  )

cat(nrow(Filtered), "final recordings after channel selection\n")

#=================================================
# Save filtered table
#=================================================

write.csv(
  Filtered,
  filtered_file,
  row.names = FALSE
)

cat("Filtered table saved to:", filtered_file, "\n")

#=================================================
# Use Filtered table for index calculation
#=================================================
ProcessTable <- Filtered

#=================================================
# Create output file
#=================================================

  header <- data.frame(
    File = character(),
    ChannelUsed = character(),
    
    RecorderID = character(),
    SiteID = character(),
    Reserve = character(),
    Date = as.Date(character()),
    Year = integer(),
    Month = integer(),
    Day = integer(),
    Hour = integer(),
    Min = integer(),
    Sec = integer(),
    
    
    ACI300 = numeric(),
    ACI2000 = numeric(),
    ADI300 = numeric(),
    ADI2000 = numeric(),
    BIO300 = numeric(),
    BIO2000 = numeric(),
    NP300 = numeric(),
    NP2000 = numeric(),
    NDSI = numeric(),
    B = numeric(),
    A = numeric(),
    NoiseFloor_dB = numeric(),
    SignalLevel_dB = numeric(),
    SNR_dB= numeric(),
    
    stringsAsFactors = FALSE
  )
  
  #=================================================
  # Create output CSV only when absent
  #=================================================
  
  if(!file.exists(outfile)){
    
    write.csv(
      header,
      outfile,
      row.names = FALSE
    )
}

#=================================================
# Completed files
#=================================================

completed <- read.csv(
  outfile,
  stringsAsFactors = FALSE
)$File

cat(
  length(completed),
  "files already processed\n"
)


#=========================================
# Seewave- Indices Functions
#=========================================    



# Number of Peaks (NP)

NP_seewave <- function(
    Recording,
    min_freq = 300,
    max_freq = 6000,
    fft_w = 1024,
    amp_threshold = c(1 / 1000, 1 / 1000)
){
  
  spec <- seewave::meanspec(
    Recording,
    f = Recording@samp.rate,
    wl = fft_w,
    from = min_freq / 1000,
    to = max_freq / 1000,
    plot = FALSE
  )
  
  peaks <- seewave::fpeaks(
    spec,
    amp = amp_threshold,
    plot = FALSE
  )
  
  if(is.null(peaks) || length(peaks) == 0){
    return(0)
  }
  
  return(nrow(peaks))
}



#SNR

estimate_snr <- function(
    Recording,
    min_freq = 300,
    max_freq = 5900,
    fft_w = 1024,
    overlap = 75,
    noise_quantile = 0.10,
    signal_quantile = 0.95
){
  
  # Calculate spectrogram in dB relative to its maximum
  sp <- seewave::spectro(
    wave = Recording,
    f = Recording@samp.rate,
    wl = fft_w,
    ovlp = overlap,
    dB = "max0",
    norm = TRUE,
    plot = FALSE
  )
  
  # Frequency values returned by spectro() are in kHz
  frequency_hz <- sp$freq * 1000
  
  # Retain analysis-frequency range
  keep <- frequency_hz >= min_freq &
    frequency_hz <= max_freq
  
  amplitude_db <- sp$amp[
    keep,
    ,
    drop = FALSE
  ]
  
  # Remove invalid values
  valid_values <- amplitude_db[
    is.finite(amplitude_db)
  ]
  
  if(length(valid_values) == 0){
    return(
      data.frame(
        NoiseFloor_dB = NA_real_,
        SignalLevel_dB = NA_real_,
        SNR_dB = NA_real_
      )
    )
  }
  
  # Estimate background and signal levels
  noise_floor <- as.numeric(
    quantile(
      valid_values,
      probs = noise_quantile,
      na.rm = TRUE
    )
  )
  
  signal_level <- as.numeric(
    quantile(
      valid_values,
      probs = signal_quantile,
      na.rm = TRUE
    )
  )
  
  snr_db <- signal_level - noise_floor
  
  return(
    data.frame(
      NoiseFloor_dB = noise_floor,
      SignalLevel_dB = signal_level,
      SNR_dB = snr_db
    )
  )
}


#=================================================
# Full Table from Name
#=================================================

#=================================================
# Recorder-to-site lookup
#=================================================

RecorderLookup <- data.frame(
  RecorderID = c(
    "S4A21787",
    "S4A23764",
    "S4A23770",
    "S4A23775",
    "S4A23837"
  ),
  
  SiteID = c(
    "IARR_A_SouthRK",
    "IARR_A_NorthRK",
    "KKRR_B_Falaq",
    "KKRR_A_Catchment",
    "KKRR_B_Wadi"
  ),
  
  Reserve = c(
    "IARR",
    "IARR",
    "KKRR",
    "KKRR",
    "KKRR"
  ),
  
  stringsAsFactors = FALSE
)

#=================================================
# Function to extract metadata from filename
#=================================================

parse_filename_metadata <- function(file_path, recorder_lookup){
  
  # Remove directory and file extension
  file_name <- basename(file_path)
  file_stem <- tools::file_path_sans_ext(file_name)
  
  # Split filename at underscores
  parts <- strsplit(
    file_stem,
    split = "_",
    fixed = TRUE
  )[[1]]
  
  # Expected structure:
  # RecorderID_Channel_Date_Time_000
  if(length(parts) < 5){
    
    stop(
      paste(
        "Unexpected filename structure:",
        file_name
      )
    )
  }
  
  recorder_id <- parts[1]
  recording_date <- parts[3]
  recording_time <- parts[4]
  
  # Validate date and time strings
  if(nchar(recording_date) != 8){
    
    stop(
      paste(
        "Date does not contain eight digits:",
        file_name
      )
    )
  }
  
  if(nchar(recording_time) != 6){
    
    stop(
      paste(
        "Time does not contain six digits:",
        file_name
      )
    )
  }
  
  # Extract date components
  year_value <- as.integer(
    substr(recording_date, 1, 4)
  )
  
  month_value <- as.integer(
    substr(recording_date, 5, 6)
  )
  
  day_value <- as.integer(
    substr(recording_date, 7, 8)
  )
  
  # Extract time components
  hour_value <- as.integer(
    substr(recording_time, 1, 2)
  )
  
  minute_value <- as.integer(
    substr(recording_time, 3, 4)
  )
  
  second_value <- as.integer(
    substr(recording_time, 5, 6)
  )
  
  # YYYYMMDD 
  date_value <- recording_date
    format = "%Y%m%d"
  
  
  # Match recorder to site and reserve
  lookup_row <- recorder_lookup[
    recorder_lookup$RecorderID == recorder_id,
    ,
    drop = FALSE
  ]
  
  if(nrow(lookup_row) == 0){
    
    warning(
      paste(
        "RecorderID was not found in RecorderLookup:",
        recorder_id
      )
    )
    
    site_id <- NA_character_
    reserve_value <- NA_character_
    
  } else {
    
    site_id <- lookup_row$SiteID[1]
    reserve_value <- lookup_row$Reserve[1]
  }
  
  data.frame(
    RecorderID = recorder_id,
    SiteID = site_id,
    Reserve = reserve_value,
    Date = date_value,
    Year = year_value,
    Month = month_value,
    Day = day_value,
    Hour = hour_value,
    Min = minute_value,
    Sec = second_value,
    stringsAsFactors = FALSE
  )
}


#=================================================
# Main Loop
#=================================================

for(i in seq_len(nrow(ProcessTable)))
{
  
  current_file <- ProcessTable$File[i]
  
  current_channel <- ProcessTable$ChannelUsed[i]
  
  channel_name <- ifelse(
    current_channel == 0,
    "Left",
    "Right"
  )
  
  if(current_file %in% completed)
  {
    next
  }
  
  cat(
    "\nProcessing",
    i,
    "of",
    nrow(ProcessTable),
    "\n"
  )
  
  cat(
    basename(current_file),
    "\n"
  )
  
  cat(
    "Channel:",
    channel_name,
    "\n"
  )
  
  tryCatch({
    
    #=========================================
    # Read WAV
    #=========================================
    
    Recording <- readWave(
      current_file,
      from = 0,
      to = 60,
      units = "seconds"
    )
    
    #=========================================
    # Sampling-rate standardization
    #=========================================
    
    if(Recording@samp.rate == 24000)
    {
      Recording <- downsample(
        Recording,
        12000
      )
      
      cat(
        "Downsampled to 12000 Hz\n"
      )
    }
    

    #=========================================
    # Extract metadata from filename
    #=========================================
    
    metadata <- parse_filename_metadata(
      file_path = current_file,
      recorder_lookup = RecorderLookup
    )
    
    cat(
      "Recorder:",
      metadata$RecorderID,
      "| Site:",
      metadata$SiteID,
      "| Reserve:",
      metadata$Reserve,
      "| Date:",
      as.character(metadata$Date),
      "| Time:",
      sprintf(
        "%02d:%02d:%02d",
        metadata$Hour,
        metadata$Min,
        metadata$Sec
      ),
      "\n"
    )
 
    
    #=========================================
    # Calculate Indices
    #=========================================
    
    ACI300 <- acoustic_complexity(
      Recording,
      min_freq = 300,
      max_freq = 6000, fft_w = 1024
    )$AciTotAll_left
    
    ACI2000 <- acoustic_complexity(
      Recording,
      min_freq = 2000,
      max_freq = 6000, fft_w = 1024
    )$AciTotAll_left
    
    Recording_ADI300 <- seewave::ffilter(
      wave = Recording,
      f = Recording@samp.rate,
      from = 300,
      to = 5999,
      bandpass = TRUE,
      wl = 1024,
      ovlp = 75,
      output = "Wave"
    )
    
    ADI300 <- acoustic_diversity(
      Recording_ADI300,
      max_freq = 6000,
      db_threshold = -30,
      freq_step = 500
    )$adi_left
    
    Recording_ADI2000 <- seewave::ffilter(
      wave = Recording,
      f = Recording@samp.rate,
      from = 2000,
      to = 5999,
      bandpass = TRUE,
      wl = 1024,
      ovlp = 75,
      output = "Wave"
    )
    
    ADI2000 <- acoustic_diversity(
      Recording_ADI2000,
      max_freq = 6000,
      db_threshold = -30,
      freq_step = 500
    )$adi_left
    

    BIO300 <- bioacoustic_index(
      Recording,
      min_freq = 300,
      max_freq = 6000,
      fft_w = 1024
    )$left_area
    
    BIO2000 <- bioacoustic_index(
      Recording,
      min_freq = 2000,
      max_freq = 6000,
      fft_w = 1024
    )$left_area
    
    NP300 <- NP_seewave(
      Recording,
      min_freq = 300,
      max_freq = 6000,
      fft_w = 1024
    )
    
    NP2000 <- NP_seewave(
      Recording,
      min_freq = 2000,
      max_freq = 6000,
      fft_w = 1024
    )
    
    ndsi <- ndsi(
      Recording,
      fft_w = 1024,
      anthro_min = 300,
      anthro_max = 2000,
      bio_min = 2000,
      bio_max = 6000
    )
    
    snr_result <- estimate_snr(
      Recording = Recording,
      min_freq = 300,
      max_freq = 5900,
      fft_w = 1024,
      overlap = 75,
      noise_quantile = 0.10,
      signal_quantile = 0.95
    )
    
    NoiseFloor_dB <- snr_result$NoiseFloor_dB
    SignalLevel_dB <- snr_result$SignalLevel_dB
    SNR_dB <- snr_result$SNR_dB
    
    #=========================================
    # Results
    #=========================================
    
    result <- data.frame(
      File = current_file,
      ChannelUsed = channel_name,
      
      
      RecorderID = metadata$RecorderID,
      SiteID = metadata$SiteID,
      Reserve = metadata$Reserve,
      Date = as.character(metadata$Date),
      Year = metadata$Year,
      Month = metadata$Month,
      Day = metadata$Day,
      Hour = metadata$Hour,
      Min = metadata$Min,
      Sec = metadata$Sec,
      
      
      ACI300 = ACI300,
      ACI2000 = ACI2000,
      ADI300 = ADI300,
      ADI2000 = ADI2000,
      BIO300 = BIO300,
      BIO2000 = BIO2000,
      NP300 = NP300,
      NP2000 = NP2000,
      NDSI = ndsi$ndsi_left,
      B = ndsi$biophony_left,
      A = ndsi$anthrophony_left,
      NoiseFloor_dB = NoiseFloor_dB,
      SignalLevel_dB = SignalLevel_dB,
      SNR_dB = SNR_dB
    )
    
    
    
    #=========================================
    # Append to CSV immediately
    #=========================================
    
    # check columns
    expected_columns <- names(header)
    
    if(!identical(names(result), expected_columns)){
      
      stop(
        paste0(
          "Result columns do not match the header.\n",
          "Result columns: ",
          paste(names(result), collapse = ", "),
          "\nHeader columns: ",
          paste(expected_columns, collapse = ", ")
        )
      )
    }
    
    # Write table
    
    write.table(
      result,
      file = outfile,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
    
    
    #=========================================
    # Update processed list
    #=========================================
    
    completed <- unique(
      c(
        completed,
        current_file
      )
    )
    
    #=========================================
    # Save checkpoint after every successful file
    #=========================================
    
    save(
      completed,
      i,
      current_file,
      file = checkpoint_file
    )
    
    if(!file.exists(checkpoint_file)){
      warning(
        paste(
          "Checkpoint was not created:",
          checkpoint_file
        )
      )
    } else {
      cat(
        "Checkpoint saved after",
        length(completed),
        "completed files\n"
      )
    }
    
    cat("Finished\n")
    
  },  # closes the successful portion of tryCatch()
  
  error = function(e)
  {
    error_message <- paste0(
      format(Sys.time()),
      " | ERROR: ",
      basename(current_file),
      " | ",
      conditionMessage(e),
      "\n"
    )
    
    cat(error_message)
    
    cat(
      error_message,
      file = file.path(
        dirname(outfile),
        "AcousticIndices_errors.txt"
      ),
      append = TRUE
    )
  }
  
  )  # closes tryCatch()
  
}  # closes the main for loop

#=================================================
# Final checkpoint
#=================================================

save(
  completed,
  file = checkpoint_file
)

if(file.exists(checkpoint_file))
{
  cat(
    "\nFinal checkpoint saved to:\n",
    checkpoint_file,
    "\n"
  )
}
else
{
  warning(
    paste(
      "Final checkpoint was not created:",
      checkpoint_file
    )
  )
}

cat(
  "\nAll processing completed\n"
)