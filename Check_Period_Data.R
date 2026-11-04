library(dplyr)
library(tidyr)
library(suncalc)

# ============================================================
# 1. PATHS
# ============================================================

index_file <- paste0(
  "D:/Flor/SaudiArabia/Data/Results/Indices/",
  "SaudiArabia_Allindices_Oct2026_multiplefreqBins.csv"
)

output_directory <- paste0(
  "D:/Flor/SaudiArabia/Data/Results/",
  "Coverage"
)

if(!dir.exists(output_directory)){
  dir.create(
    output_directory,
    recursive = TRUE
  )
}

# ============================================================
# 2. READ MINUTE-LEVEL ACOUSTIC-INDEX RESULTS
# ============================================================

indices <- read.csv(
  index_file,
  stringsAsFactors = FALSE
)

# Confirm required columns are available
required_columns <- c(
  "File",
  "RecorderID",
  "SiteID",
  "Reserve",
  "Date",
  "Year",
  "Month",
  "Day",
  "Hour",
  "Min"
)

missing_columns <- setdiff(
  required_columns,
  names(indices)
)

if(length(missing_columns) > 0){
  
  stop(
    paste0(
      "The following required columns are missing:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}

cat(
  nrow(indices),
  "usable minute-level recordings read\n"
)

# ============================================================
# 3. FORMAT VARIABLES
# ============================================================

indices <- indices %>%
  mutate(
    RecorderID = as.character(RecorderID),
    SiteID = as.character(SiteID),
    Reserve = as.character(Reserve),
    
    Date = as.character(Date),
    
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = as.integer(Hour),
    Min = as.integer(Min),
    
    DateValue = as.Date(
      Date,
      format = "%Y%m%d"
    )
  )

# Check date conversion
if(any(is.na(indices$DateValue))){
  
  warning(
    paste(
      sum(is.na(indices$DateValue)),
      "rows contained dates that could not be converted."
    )
  )
}

# ============================================================
# 4. KEEP MARCH AND APRIL
# ============================================================

indices <- indices %>%
  filter(
    Month %in% c(3, 4),
    !is.na(DateValue),
    !is.na(Hour)
  )

cat(
  nrow(indices),
  "usable minute-level recordings retained for March-April\n"
)

# ============================================================
# 5. REMOVE ANY ACCIDENTAL DUPLICATES
# ============================================================

# There should be no more than one selected recording for the
# same recorder, date, hour, and minute.

indices <- indices %>%
  distinct(
    RecorderID,
    SiteID,
    Reserve,
    DateValue,
    Hour,
    Min,
    .keep_all = TRUE
  )

# ============================================================
# 6. CALCULATE DAILY ASTRONOMICAL TIMES
# ============================================================

# Coordinates used in the BAM period calculations
latitude <- 24.7136
longitude <- 46.6753

astronomical_dates <- indices %>%
  distinct(DateValue) %>%
  transmute(
    date = DateValue,
    lat = latitude,
    lon = longitude
  )

SunTimes <- getSunlightTimes(
  data = astronomical_dates,
  keep = c(
    "sunrise",
    "sunset",
    "solarNoon"
  ),
  tz = "Asia/Riyadh"
)

# Convert astronomical times to decimal clock hours
to_decimal_hour <- function(x){
  
  as.numeric(format(x, "%H")) +
    as.numeric(format(x, "%M")) / 60 +
    as.numeric(format(x, "%S")) / 3600
}

SunTimes <- SunTimes %>%
  mutate(
    SunriseHour = to_decimal_hour(sunrise),
    SunsetHour = to_decimal_hour(sunset),
    SolarNoonHour = to_decimal_hour(solarNoon)
  ) %>%
  transmute(
    DateValue = date,
    SunriseHour,
    SolarNoonHour,
    SunsetHour
  )

# ============================================================
# 7. JOIN ASTRONOMICAL TIMES TO INDEX RESULTS
# ============================================================

indices <- indices %>%
  left_join(
    SunTimes,
    by = "DateValue"
  )

if(any(is.na(indices$SunriseHour))){
  
  stop(
    paste(
      sum(is.na(indices$SunriseHour)),
      "rows did not receive astronomical times."
    )
  )
}

# ============================================================
# 8. DEFINE ASTRONOMICAL PERIODS
# ============================================================

indices <- indices %>%
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

# Check period assignment
print(
  table(
    indices$Period,
    useNA = "ifany"
  )
)

# ============================================================
# 9. CREATE ONE ROW PER SITE-DATE-HOUR
# ============================================================

# Each row in the original table is one usable minute.
# This table counts the usable minutes remaining in each hour.

HourlyAvailability <- indices %>%
  group_by(
    RecorderID,
    SiteID,
    Reserve,
    DateValue,
    Year,
    Month,
    Day,
    Hour,
    Period
  ) %>%
  summarise(
    UsableMinutes = n_distinct(Min),
    .groups = "drop"
  ) %>%
  mutate(
    AvailableHour = UsableMinutes > 0,
    CompleteHour = UsableMinutes >= 5,
    MinuteCoverage = pmin(UsableMinutes / 5, 1)
  )

# Validate counts
print(
  table(
    HourlyAvailability$UsableMinutes,
    useNA = "ifany"
  )
)

# ============================================================
# 10. CALCULATE EXPECTED SITE-DATE-HOUR COMBINATIONS
# ============================================================

# To measure missingness, we need both observed and unobserved
# combinations. Build a complete grid for every site, date,
# and hour in March-April.

SiteDates <- indices %>%
  distinct(
    RecorderID,
    SiteID,
    Reserve,
    DateValue,
    Year,
    Month,
    Day
  )

ExpectedCoverage <- SiteDates %>%
  tidyr::crossing(
    Hour = 0:23
  ) %>%
  left_join(
    SunTimes,
    by = "DateValue"
  ) %>%
  mutate(
    
    Period = case_when(
      
      Hour >= SunriseHour - 2 &
        Hour <= SunriseHour + 2 ~
        "Dawn",
      
      Hour > SunriseHour + 2 &
        Hour < SolarNoonHour - 2 ~
        "Morning",
      
      Hour >= SolarNoonHour - 2 &
        Hour <= SolarNoonHour + 2 ~
        "Midday",
      
      Hour > SolarNoonHour + 2 &
        Hour < SunsetHour - 2 ~
        "Afternoon",
      
      Hour >= SunsetHour - 2 &
        Hour <= SunsetHour + 2 ~
        "Dusk",
      
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
  ) %>%
  left_join(
    HourlyAvailability %>%
      dplyr::select(
        RecorderID,
        SiteID,
        Reserve,
        DateValue,
        Hour,
        UsableMinutes
      ),
    by = c(
      "RecorderID",
      "SiteID",
      "Reserve",
      "DateValue",
      "Hour"
    )
  ) %>%
  mutate(
    UsableMinutes = replace_na(
      UsableMinutes,
      0L
    ),
    
    AvailableHour = UsableMinutes > 0,
    
    CompleteHour = UsableMinutes >= 5,
    
    MinuteCoverage = pmin(
      UsableMinutes / 5,
      1
    )
  )

# ============================================================
# 11. COVERAGE BY CLOCK HOUR AND SITE
# ============================================================

HourlyCoverageBySite <- ExpectedCoverage %>%
  group_by(
    RecorderID,
    SiteID,
    Reserve,
    Hour
  ) %>%
  summarise(
    PossibleDateHours = n(),
    
    AvailableDateHours = sum(
      AvailableHour,
      na.rm = TRUE
    ),
    
    CompleteDateHours = sum(
      CompleteHour,
      na.rm = TRUE
    ),
    
    UsableMinutes = sum(
      UsableMinutes,
      na.rm = TRUE
    ),
    
    PossibleMinutes = PossibleDateHours * 5,
    
    PercentHoursAvailable =
      100 * AvailableDateHours / PossibleDateHours,
    
    PercentHoursComplete =
      100 * CompleteDateHours / PossibleDateHours,
    
    PercentMinutesAvailable =
      100 * UsableMinutes / PossibleMinutes,
    
    .groups = "drop"
  )

# ============================================================
# 12. COVERAGE BY ASTRONOMICAL PERIOD AND SITE
# ============================================================

PeriodCoverageBySite <- ExpectedCoverage %>%
  group_by(
    RecorderID,
    SiteID,
    Reserve,
    Period
  ) %>%
  summarise(
    PossibleDateHours = n(),
    
    AvailableDateHours = sum(
      AvailableHour,
      na.rm = TRUE
    ),
    
    CompleteDateHours = sum(
      CompleteHour,
      na.rm = TRUE
    ),
    
    UsableMinutes = sum(
      UsableMinutes,
      na.rm = TRUE
    ),
    
    PossibleMinutes = PossibleDateHours * 5,
    
    DatesPossible = n_distinct(
      DateValue
    ),
    
    DatesWithData = n_distinct(
      DateValue[AvailableHour]
    ),
    
    PercentHoursAvailable =
      100 * AvailableDateHours / PossibleDateHours,
    
    PercentHoursComplete =
      100 * CompleteDateHours / PossibleDateHours,
    
    PercentMinutesAvailable =
      100 * UsableMinutes / PossibleMinutes,
    
    PercentDatesRepresented =
      100 * DatesWithData / DatesPossible,
    
    .groups = "drop"
  )

# ============================================================
# 13. COVERAGE BY ASTRONOMICAL PERIOD AND RESERVE
# ============================================================

PeriodCoverageByReserve <- ExpectedCoverage %>%
  group_by(
    Reserve,
    Period
  ) %>%
  summarise(
    PossibleDateHours = n(),
    
    AvailableDateHours = sum(
      AvailableHour,
      na.rm = TRUE
    ),
    
    CompleteDateHours = sum(
      CompleteHour,
      na.rm = TRUE
    ),
    
    UsableMinutes = sum(
      UsableMinutes,
      na.rm = TRUE
    ),
    
    PossibleMinutes = PossibleDateHours * 5,
    
    SitesExpected = n_distinct(
      SiteID
    ),
    
    SitesWithData = n_distinct(
      SiteID[AvailableHour]
    ),
    
    PercentHoursAvailable =
      100 * AvailableDateHours / PossibleDateHours,
    
    PercentHoursComplete =
      100 * CompleteDateHours / PossibleDateHours,
    
    PercentMinutesAvailable =
      100 * UsableMinutes / PossibleMinutes,
    
    .groups = "drop"
  )

# ============================================================
# 14. FLAG SITE-PERIOD COMBINATIONS WITH LOW COVERAGE
# ============================================================

# These are transparent screening criteria, not universal rules.
# Adjust them if needed after inspecting the distributions.

minimum_available_hours <- 20
minimum_percent_hours <- 50
minimum_percent_dates <- 50

PeriodCoverageBySite <- PeriodCoverageBySite %>%
  mutate(
    
    CoverageFlag = case_when(
      
      AvailableDateHours < minimum_available_hours ~
        "Insufficient: fewer than 20 available date-hours",
      
      PercentHoursAvailable < minimum_percent_hours ~
        "Insufficient: less than 50% of date-hours available",
      
      PercentDatesRepresented < minimum_percent_dates ~
        "Insufficient: represented on less than 50% of dates",
      
      TRUE ~
        "Adequate"
    )
  )

# ============================================================
# 15. CREATE RESERVE-PERIOD RECOMMENDATIONS
# ============================================================

# A reserve-period comparison is flagged if one or more of its
# constituent sites has insufficient coverage.

ReservePeriodRecommendation <- PeriodCoverageBySite %>%
  group_by(
    Reserve,
    Period
  ) %>%
  summarise(
    SitesExpected = n_distinct(SiteID),
    
    SitesAdequate = sum(
      CoverageFlag == "Adequate"
    ),
    
    MinimumAvailableHours = min(
      AvailableDateHours,
      na.rm = TRUE
    ),
    
    MinimumSiteHourCoverage = min(
      PercentHoursAvailable,
      na.rm = TRUE
    ),
    
    MinimumSiteDateCoverage = min(
      PercentDatesRepresented,
      na.rm = TRUE
    ),
    
    MeanSiteHourCoverage = mean(
      PercentHoursAvailable,
      na.rm = TRUE
    ),
    
    MeanSiteMinuteCoverage = mean(
      PercentMinutesAvailable,
      na.rm = TRUE
    ),
    
    Recommendation = ifelse(
      SitesAdequate == SitesExpected,
      "Keep",
      "Review or exclude"
    ),
    
    .groups = "drop"
  )

# ============================================================
# 16. DIRECTLY COMPARE IARR AND KKRR COVERAGE
# ============================================================

ReserveComparisonCoverage <- PeriodCoverageByReserve %>%
  dplyr::select(
    Reserve,
    Period,
    AvailableDateHours,
    PercentHoursAvailable,
    PercentMinutesAvailable
  ) %>%
  pivot_wider(
    names_from = Reserve,
    values_from = c(
      AvailableDateHours,
      PercentHoursAvailable,
      PercentMinutesAvailable
    )
  ) %>%
  mutate(
    DifferenceInHourCoverage =
      PercentHoursAvailable_IARR -
      PercentHoursAvailable_KKRR,
    
    DifferenceInMinuteCoverage =
      PercentMinutesAvailable_IARR -
      PercentMinutesAvailable_KKRR
  )

# ============================================================
# 17. SAVE RESULTS
# ============================================================

write.csv(
  HourlyAvailability,
  file.path(
    output_directory,
    "Observed_Hourly_Availability.csv"
  ),
  row.names = FALSE
)

write.csv(
  HourlyCoverageBySite,
  file.path(
    output_directory,
    "Coverage_By_ClockHour_And_Site.csv"
  ),
  row.names = FALSE
)

write.csv(
  PeriodCoverageBySite,
  file.path(
    output_directory,
    "Coverage_By_AstronomicalPeriod_And_Site.csv"
  ),
  row.names = FALSE
)

write.csv(
  PeriodCoverageByReserve,
  file.path(
    output_directory,
    "Coverage_By_AstronomicalPeriod_And_Reserve.csv"
  ),
  row.names = FALSE
)

write.csv(
  ReservePeriodRecommendation,
  file.path(
    output_directory,
    "Reserve_Period_Coverage_Recommendations.csv"
  ),
  row.names = FALSE
)

write.csv(
  ReserveComparisonCoverage,
  file.path(
    output_directory,
    "IARR_KKRR_Coverage_Comparison.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 18. DISPLAY RESULTS
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "COVERAGE BY SITE AND PERIOD\n"
)

cat(
  "========================================\n"
)

print(
  PeriodCoverageBySite %>%
    arrange(
      Reserve,
      SiteID,
      Period
    )
)

cat(
  "\n========================================\n"
)

cat(
  "RESERVE-PERIOD RECOMMENDATIONS\n"
)

cat(
  "========================================\n"
)

print(
  ReservePeriodRecommendation %>%
    arrange(
      Reserve,
      Period
    )
)

cat(
  "\n========================================\n"
)

cat(
  "IARR-KKRR COVERAGE COMPARISON\n"
)

cat(
  "========================================\n"
)

print(
  ReserveComparisonCoverage
)