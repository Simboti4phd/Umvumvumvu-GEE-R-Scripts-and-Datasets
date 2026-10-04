# ============================================================
# SPI-3 CALCULATION USING GAMMA DISTRIBUTION
# UMVUMVUMVU WATERSHED, ZIMBABWE
#
# Study Period:  January 1989 - December 2018
# SPI Timescale: 3 months (SPI-3)
# Distribution:  Gamma
# Zero handling: Center of Mass (Stagge et al., 2015)
#
# ============================================================


# ============================================================
# 1. INSTALL REQUIRED PACKAGES (run only if not installed)
# ============================================================

# install.packages("readr")
# install.packages("dplyr")
# install.packages("lubridate")
# install.packages("ggplot2")
# install.packages("zoo")


# ============================================================
# 2. LOAD PACKAGES
# ============================================================

library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(zoo)      # for rollapply()

# ============================================================
# 3. IMPORT MONTHLY PRECIPITATION DATA
# ============================================================

data <- read_csv(
  "Monthly_PrecipitationCM1.csv",
  show_col_types = FALSE
)

# ============================================================
# 4. INSPECT THE DATA
# ============================================================

print(head(data))
print(names(data))
str(data)


# ============================================================
# 5. STANDARDISE PRECIPITATION COLUMN
# ============================================================

data <- data %>%
  rename(
    Precipitation = Precipitation_mm
  )


# ============================================================
# 6. CONVERT VARIABLES TO APPROPRIATE DATA TYPES
# ============================================================

data <- data %>%
  mutate(
    Year          = as.numeric(Year),
    Month         = as.numeric(Month),
    Precipitation = as.numeric(Precipitation)
  )


# ============================================================
# 7. CREATE DATE VARIABLE
# ============================================================

data <- data %>%
  mutate(
    Date = make_date(
      year  = Year,
      month = Month,
      day   = 1
    )
  )


# ============================================================
# 8. SORT DATA CHRONOLOGICALLY
# ============================================================

data <- data %>%
  arrange(Date)


# ============================================================
# 9. CHECK STUDY PERIOD
# ============================================================

cat("\n============================================\n")
cat("STUDY PERIOD CHECK\n")
cat("============================================\n")

cat("Start date:",
    format(min(data$Date, na.rm = TRUE), "%Y-%m-%d"), "\n")

cat("End date:",
    format(max(data$Date, na.rm = TRUE), "%Y-%m-%d"), "\n")


# ============================================================
# 10. CHECK EXPECTED NUMBER OF MONTHS
# ============================================================

expected_months <- seq(
  from = as.Date("1989-01-01"),
  to   = as.Date("2018-12-01"),
  by   = "month"
)

cat("\nExpected number of months:",
    length(expected_months), "\n")

cat("Number of rows in CSV:",
    nrow(data), "\n")


# ============================================================
# 11. CHECK FOR MISSING MONTHS
# ============================================================

missing_dates <- expected_months[
  !expected_months %in% data$Date
]

if (length(missing_dates) == 0) {
  cat("\nNo months are missing from the dataset.\n")
} else {
  cat("\nWARNING:", length(missing_dates), "months are missing.\n")
  print(missing_dates)
}


# ============================================================
# 12. CHECK FOR DUPLICATE MONTHS
# ============================================================

duplicate_dates <- data$Date[duplicated(data$Date)]

if (length(duplicate_dates) == 0) {
  cat("\nNo duplicate months detected.\n")
} else {
  cat("\nWARNING:", length(duplicate_dates), "duplicate months detected.\n")
  print(duplicate_dates)
}


# ============================================================
# 13. CHECK FOR MISSING PRECIPITATION VALUES
# ============================================================

missing_precip <- sum(is.na(data$Precipitation))

cat("\nNumber of missing precipitation values:",
    missing_precip, "\n")


# ============================================================
# 14. CHECK FOR NEGATIVE PRECIPITATION VALUES
# ============================================================

negative_precip <- sum(data$Precipitation < 0, na.rm = TRUE)

cat("Number of negative precipitation values:",
    negative_precip, "\n")


# ============================================================
# 15. STOP IF DATA QUALITY PROBLEMS EXIST
# ============================================================

if (length(missing_dates) > 0) {
  stop("\nThe dataset contains missing months. ",
       "Correct the monthly time series before calculating SPI.")
}

if (length(duplicate_dates) > 0) {
  stop("\nThe dataset contains duplicate months. ",
       "Correct the duplicate records before calculating SPI.")
}

if (missing_precip > 0) {
  stop("\nThe dataset contains missing precipitation values. ",
       "The SPI calculation requires complete data.")
}

if (negative_precip > 0) {
  stop("\nNegative precipitation values were detected. ",
       "Check the precipitation data before calculating SPI.")
}


# ============================================================
# 16. CREATE 3-MONTH ROLLING PRECIPITATION TOTAL
# ============================================================
#
# SPI-3 is based on the 3-month accumulated precipitation.
# e.g., the value for March 1989 = Jan + Feb + Mar 1989.
#
# The first 2 months of the series will be NA
# because there aren't yet 3 months to sum.

data <- data %>%
  mutate(
    Precip_3M = rollapply(
      Precipitation,
      width     = 3,
      FUN       = sum,
      align     = "right",
      fill      = NA
    )
  )


# ============================================================
# 17. CREATE 3-MONTH PRECIPITATION TIME SERIES
# ============================================================
#
# Note: We drop the first 2 NAs so the ts object starts
# at March 1989 (the first month with a full 3-month total).

precip_3m_values <- data$Precip_3M
valid_idx        <- which(!is.na(precip_3m_values))

precip_ts_3m <- ts(
  precip_3m_values[valid_idx],
  start     = c(
    year(data$Date[valid_idx[1]]),
    month(data$Date[valid_idx[1]])
  ),
  frequency = 12
)


# ============================================================
# 18. VERIFY THE TIME-SERIES OBJECT
# ============================================================

cat("\n============================================\n")
cat("3-MONTH TIME SERIES CHECK\n")
cat("============================================\n")

cat("Start:",
    paste(start(precip_ts_3m), collapse = "-"), "\n")

cat("End:",
    paste(end(precip_ts_3m), collapse = "-"), "\n")

cat("Frequency:",
    frequency(precip_ts_3m), "\n")

cat("Number of observations:",
    length(precip_ts_3m), "\n")


# ============================================================
# 19. SPI-3 USING GAMMA + CENTER OF MASS ZERO HANDLING
# ============================================================
#
# Method (Stagge et al., 2015) — same as SPI-1 but on
# 3-month accumulated precipitation.
#
# For each calendar month (Jan, Feb, ..., Dec):
#
#   1. q = probability of zero 3-month precipitation.
#   2. Fit gamma to the POSITIVE 3-month totals only.
#   3. Assign cumulative probabilities:
#
#        H(x) = q / 2                  if x == 0
#        H(x) = q + (1 - q) * G(x)     if x  > 0
#
#   4. Clamp H(x) to [0.0001, 0.9999] to avoid ±Inf.
#   5. Transform with qnorm() to get SPI.
#
# ============================================================

calculate_spi_center_of_mass <- function(precip_ts) {
  
  x      <- as.numeric(precip_ts)
  n      <- length(x)
  months <- cycle(precip_ts)
  
  spi_values <- rep(NA_real_, n)
  
  for (m in 1:12) {
    
    idx <- which(months == m)
    xm  <- x[idx]
    
    if (length(xm) < 3) next
    
    q <- mean(xm == 0)
    
    positive <- xm[xm > 0]
    
    # All-zero case
    if (length(positive) == 0) {
      H <- rep(q / 2, length(xm))
      H <- pmin(pmax(H, 1e-4), 1 - 1e-4)
      spi_values[idx] <- qnorm(H)
      next
    }
    
    mu     <- mean(positive)
    sigma2 <- var(positive)
    
    if (is.na(mu) || is.na(sigma2) || sigma2 <= 0 || mu <= 0) {
      next
    }
    
    alpha <- mu^2 / sigma2
    beta  <- sigma2 / mu
    
    H <- numeric(length(xm))
    
    for (i in seq_along(xm)) {
      
      if (xm[i] == 0) {
        H[i] <- q / 2
      } else {
        H[i] <- q + (1 - q) * pgamma(xm[i], shape = alpha, scale = beta)
      }
      
    }
    
    H <- pmin(pmax(H, 1e-4), 1 - 1e-4)
    
    spi_values[idx] <- qnorm(H)
    
  }
  
  return(spi_values)
}


# ============================================================
# 20. APPLY THE FUNCTION AND EXTRACT SPI-3
# ============================================================

spi3_values <- calculate_spi_center_of_mass(precip_ts_3m)

# Align back to the full data frame (first 2 months will be NA)
data$SPI_3 <- NA_real_
data$SPI_3[valid_idx] <- spi3_values

# Round for Excel compatibility
data <- data %>%
  mutate(
    SPI_3 = round(SPI_3, 3)
  )


# ============================================================
# 21. CLASSIFY SPI-3
# ============================================================

data <- data %>%
  mutate(
    SPI_3_Classification = case_when(
      SPI_3 >=  2.00 ~ "Extremely Wet",
      SPI_3 >=  1.50 ~ "Very Wet",
      SPI_3 >=  1.00 ~ "Moderately Wet",
      SPI_3 >  -1.00 ~ "Near Normal",
      SPI_3 >  -1.50 ~ "Moderately Dry",
      SPI_3 >  -2.00 ~ "Severely Dry",
      SPI_3 <= -2.00 ~ "Extremely Dry",
      TRUE           ~ NA_character_
    )
  )


# ============================================================
# 22. CREATE FINAL SPI-3 DATASET
# ============================================================

SPI3_data <- data %>%
  dplyr::select(
    Date,
    Year,
    Month,
    Precipitation,
    Precip_3M,
    SPI_3,
    SPI_3_Classification
  )


# ============================================================
# 23. DISPLAY FIRST 12 SPI-3 VALUES
# ============================================================

cat("\n============================================\n")
cat("FIRST 12 SPI-3 VALUES\n")
cat("============================================\n")

print(head(SPI3_data, 12))


# ============================================================
# 24. DISPLAY LAST 12 SPI-3 VALUES
# ============================================================

cat("\n============================================\n")
cat("LAST 12 SPI-3 VALUES\n")
cat("============================================\n")

print(tail(SPI3_data, 12))


# ============================================================
# 25. SPI-3 SUMMARY STATISTICS
# ============================================================

cat("\n============================================\n")
cat("SPI-3 SUMMARY STATISTICS\n")
cat("============================================\n")

print(summary(SPI3_data$SPI_3))

cat("\nSPI-3 Classification Counts:\n")
print(table(SPI3_data$SPI_3_Classification, useNA = "ifany"))


# ============================================================
# 26. PLOT SPI-3 TIME SERIES
# ============================================================

par(mar = c(4, 4, 3, 1))

plot(
  x    = data$Date,
  y    = data$SPI_3,
  type = "h",
  col  = ifelse(
    data$SPI_3 < 0,
    "firebrick",
    "steelblue"
  ),
  lwd  = 2,
  xlab = "Date",
  ylab = "SPI-3",
  main = "Gamma-Distribution SPI-3 (Center of Mass)\nUmvumvumvu Watershed (1989-2018)"
)

abline(h =  0,    lty = 1, col = "black")
abline(h = -1,    lty = 2, col = "orange")
abline(h = -1.5,  lty = 2, col = "red")
abline(h = -2,    lty = 2, col = "darkred")
abline(h =  1,    lty = 2, col = "orange")
abline(h =  1.5,  lty = 2, col = "blue")
abline(h =  2,    lty = 2, col = "darkblue")

legend(
  "topright",
  legend = c("Wet (SPI > 0)", "Dry (SPI < 0)", "+/-1", "+/-1.5", "+/-2"),
  col  = c("steelblue", "firebrick", "orange", "red", "darkred"),
  lty  = c(NA, NA, 2, 2, 2),
  pch  = c(15, 15, NA, NA, NA),
  bty  = "n",
  cex  = 0.8
)


# ============================================================
# 27. SAVE THE PLOT AS A PNG
# ============================================================

png(
  filename = "C:/Users/Chidzaa/Documents/Mutambara/Objective 2/SPI3_Umvumvumvu_1989_2018.png",
  width    = 1600,
  height   = 900,
  res      = 150
)

plot(
  x    = data$Date,
  y    = data$SPI_3,
  type = "h",
  col  = ifelse(
    data$SPI_3 < 0,
    "firebrick",
    "steelblue"
  ),
  lwd  = 2,
  xlab = "Date",
  ylab = "SPI-3",
  main = "Gamma-Distribution SPI-3 (Center of Mass)\nUmvumvumvu Watershed (1989-2018)"
)

abline(h =  0,    lty = 1)
abline(h = -1,    lty = 2, col = "orange")
abline(h = -1.5,  lty = 2, col = "red")
abline(h = -2,    lty = 2, col = "darkred")
abline(h =  1,    lty = 2, col = "orange")
abline(h =  1.5,  lty = 2, col = "blue")
abline(h =  2,    lty = 2, col = "darkblue")

dev.off()

cat("\nPlot saved as PNG.\n")


# ============================================================
# 28. EXPORT THE FINAL SPI-3 DATASET
# ============================================================

SPI3_output <- data %>%
  dplyr::select(
    Date,
    Year,
    Month,
    Precipitation,
    Precip_3M,
    SPI_3,
    SPI_3_Classification
  )

write_csv(
  SPI3_output,
  "C:/Users/Chidzaa/Documents/Mutambara/Objective 2/SPI_3_Umvumvumvu_1989_2018.csv"
)

cat("\nSPI-3 data exported successfully.\n")


# ============================================================
# 29. EXPORT SUMMARY / SANITY CHECK
# ============================================================

cat("\n--- EXPORT SUMMARY ---\n")
cat("Rows exported:  ", nrow(SPI3_output), "\n")
cat("Columns exported:", ncol(SPI3_output), "\n")

cat("\nFirst 6 rows:\n")
print(head(SPI3_output))

cat("\nSPI-3 summary:\n")
print(summary(SPI3_output$SPI_3))

cat("\nClassification counts:\n")
print(table(SPI3_output$SPI_3_Classification, useNA = "ifany"))

cat("\nMonths with 3-month precipitation = 0:\n")
print(sum(SPI3_output$Precip_3M == 0, na.rm = TRUE))

cat("\nAny remaining NA in SPI_3? (expected: 2, at the start)\n")
print(sum(is.na(SPI3_output$SPI_3)))

# ============================================================
# END OF SCRIPT
# ============================================================