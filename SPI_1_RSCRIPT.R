# ============================================================
# SPI-1 CALCULATION USING GAMMA DISTRIBUTION
# UMVUMVUMVU WATERSHED, ZIMBABWE
#
# Study Period:  January 1989 - December 2018
# SPI Timescale: 1 month (SPI-1)
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

# ============================================================
# 2. LOAD PACKAGES
# ============================================================

library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)


# ============================================================
# 3. IMPORT MONTHLY PRECIPITATION DATA
# ============================================================

data <- read_csv(
  "Monthly_PrecipitationCM.csv",
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
# Assumed CSV structure:
#   Year | Month | Precipitation_mm
#
# We rename Precipitation_mm -> Precipitation for clarity.

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

cat(
  "Start date:",
  format(min(data$Date, na.rm = TRUE), "%Y-%m-%d"),
  "\n"
)

cat(
  "End date:",
  format(max(data$Date, na.rm = TRUE), "%Y-%m-%d"),
  "\n"
)


# ============================================================
# 10. CHECK EXPECTED NUMBER OF MONTHS
# ============================================================

expected_months <- seq(
  from = as.Date("1989-01-01"),
  to   = as.Date("2018-12-31"),
  by   = "month"
)

cat(
  "\nExpected number of months:",
  length(expected_months),
  "\n"
)

cat(
  "Number of rows in CSV:",
  nrow(data),
  "\n"
)


# ============================================================
# 11. CHECK FOR MISSING MONTHS
# ============================================================

missing_dates <- expected_months[
  !expected_months %in% data$Date
]

if (length(missing_dates) == 0) {
  
  cat("\nNo months are missing from the dataset.\n")
  
} else {
  
  cat(
    "\nWARNING:",
    length(missing_dates),
    "months are missing.\n"
  )
  
  print(missing_dates)
  
}


# ============================================================
# 12. CHECK FOR DUPLICATE MONTHS
# ============================================================

duplicate_dates <- data$Date[
  duplicated(data$Date)
]

if (length(duplicate_dates) == 0) {
  
  cat("\nNo duplicate months detected.\n")
  
} else {
  
  cat(
    "\nWARNING:",
    length(duplicate_dates),
    "duplicate months detected.\n"
  )
  
  print(duplicate_dates)
  
}


# ============================================================
# 13. CHECK FOR MISSING PRECIPITATION VALUES
# ============================================================

missing_precip <- sum(
  is.na(data$Precipitation)
)

cat(
  "\nNumber of missing precipitation values:",
  missing_precip,
  "\n"
)


# ============================================================
# 14. CHECK FOR NEGATIVE PRECIPITATION VALUES
# ============================================================

negative_precip <- sum(
  data$Precipitation < 0,
  na.rm = TRUE
)

cat(
  "Number of negative precipitation values:",
  negative_precip,
  "\n"
)


# ============================================================
# 15. STOP IF DATA QUALITY PROBLEMS EXIST
# ============================================================

if (length(missing_dates) > 0) {
  stop(
    "\nThe dataset contains missing months. ",
    "Correct the monthly time series before calculating SPI."
  )
}

if (length(duplicate_dates) > 0) {
  stop(
    "\nThe dataset contains duplicate months. ",
    "Correct the duplicate records before calculating SPI."
  )
}

if (missing_precip > 0) {
  stop(
    "\nThe dataset contains missing precipitation values. ",
    "The SPI calculation requires complete data."
  )
}

if (negative_precip > 0) {
  stop(
    "\nNegative precipitation values were detected. ",
    "Check the precipitation data before calculating SPI."
  )
}


# ============================================================
# 16. CREATE MONTHLY TIME-SERIES OBJECT
# ============================================================

precip_ts <- ts(
  data$Precipitation,
  start     = c(1989, 1),
  frequency = 12
)


# ============================================================
# 17. VERIFY THE TIME-SERIES OBJECT
# ============================================================

cat("\n============================================\n")
cat("TIME SERIES CHECK\n")
cat("============================================\n")

cat("Start:",
    paste(start(precip_ts), collapse = "-"), "\n")

cat("End:",
    paste(end(precip_ts), collapse = "-"), "\n")

cat("Frequency:",
    frequency(precip_ts), "\n")

cat("Number of observations:",
    length(precip_ts), "\n")


# ============================================================
# 18. SPI-1 USING GAMMA + CENTER OF MASS ZERO HANDLING
# ============================================================
#
# Method (Stagge et al., 2015):
#
#   For each calendar month (Jan, Feb, ..., Dec):
#
#   1. Calculate q = probability of zero precipitation.
#
#   2. Fit a gamma distribution to the POSITIVE values only.
#
#   3. Assign cumulative probabilities:
#
#        H(x) = q / 2                  if x == 0   (center of mass)
#        H(x) = q + (1 - q) * G(x)     if x  > 0
#
#   4. Clamp H(x) to [0.0001, 0.9999] to avoid ±Inf.
#
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
    
    # Probability of zero precipitation in this calendar month
    q <- mean(xm == 0)
    
    # Positive values for gamma fit
    positive <- xm[xm > 0]
    
    # Handle the all-zeros case
    if (length(positive) == 0) {
      H <- rep(q / 2, length(xm))
      H <- pmin(pmax(H, 1e-4), 1 - 1e-4)
      spi_values[idx] <- qnorm(H)
      next
    }
    
    # Method of moments for gamma parameters
    mu     <- mean(positive)
    sigma2 <- var(positive)
    
    if (is.na(mu) || is.na(sigma2) || sigma2 <= 0 || mu <= 0) {
      next
    }
    
    alpha <- mu^2 / sigma2
    beta  <- sigma2 / mu
    
    # Cumulative probability (center of mass)
    H <- numeric(length(xm))
    
    for (i in seq_along(xm)) {
      
      if (xm[i] == 0) {
        H[i] <- q / 2
      } else {
        H[i] <- q + (1 - q) * pgamma(xm[i], shape = alpha, scale = beta)
      }
      
    }
    
    # Clamp to avoid ±Inf
    H <- pmin(pmax(H, 1e-4), 1 - 1e-4)
    
    spi_values[idx] <- qnorm(H)
    
  }
  
  return(spi_values)
}


# ============================================================
# 19. APPLY THE FUNCTION AND EXTRACT SPI-1
# ============================================================

data$SPI_1 <- calculate_spi_center_of_mass(precip_ts)

# Round for Excel compatibility (avoids #NAME?)
data <- data %>%
  mutate(
    SPI_1 = round(SPI_1, 3)
  )


# ============================================================
# 20. CLASSIFY SPI-1
# ============================================================

data <- data %>%
  mutate(
    SPI_Classification = case_when(
      SPI_1 >=  2.00 ~ "Extremely Wet",
      SPI_1 >=  1.50 ~ "Very Wet",
      SPI_1 >=  1.00 ~ "Moderately Wet",
      SPI_1 >  -1.00 ~ "Near Normal",
      SPI_1 >  -1.50 ~ "Moderately Dry",
      SPI_1 >  -2.00 ~ "Severely Dry",
      SPI_1 <= -2.00 ~ "Extremely Dry",
      TRUE           ~ NA_character_
    )
  )


# ============================================================
# 21. CREATE FINAL SPI DATASET
# ============================================================

SPI_data <- data %>%
  dplyr::select(
    Date,
    Year,
    Month,
    Precipitation,
    SPI_1,
    SPI_Classification
  )


# ============================================================
# 22. DISPLAY FIRST 12 SPI-1 VALUES
# ============================================================

cat("\n============================================\n")
cat("FIRST 12 SPI-1 VALUES\n")
cat("============================================\n")

print(head(SPI_data, 12))


# ============================================================
# 23. DISPLAY LAST 12 SPI-1 VALUES
# ============================================================

cat("\n============================================\n")
cat("LAST 12 SPI-1 VALUES\n")
cat("============================================\n")

print(tail(SPI_data, 12))


# ============================================================
# 24. SPI-1 SUMMARY STATISTICS
# ============================================================

cat("\n============================================\n")
cat("SPI-1 SUMMARY STATISTICS\n")
cat("============================================\n")

print(summary(SPI_data$SPI_1))

cat("\nSPI-1 Classification Counts:\n")
print(table(SPI_data$SPI_Classification, useNA = "ifany"))


# ============================================================
# 25. PLOT SPI-1 TIME SERIES
# ============================================================

par(mar = c(4, 4, 3, 1))

plot(
  x    = data$Date,
  y    = data$SPI_1,
  type = "h",
  col  = ifelse(
    data$SPI_1 < 0,
    "firebrick",
    "steelblue"
  ),
  lwd  = 2,
  xlab = "Date",
  ylab = "SPI-1",
  main = "Gamma-Distribution SPI-1 (Center of Mass)\nUmvumvumvu Watershed (1989-2018)"
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
  legend = c(
    "Wet (SPI > 0)",
    "Dry (SPI < 0)",
    "+/-1", "+/-1.5", "+/-2"
  ),
  col  = c("steelblue", "firebrick", "orange", "red", "darkred"),
  lty  = c(NA, NA, 2, 2, 2),
  pch  = c(15, 15, NA, NA, NA),
  bty  = "n",
  cex  = 0.8
)


# ============================================================
# 26. SAVE THE PLOT AS A PNG
# ============================================================

png(
  filename = "C:/Users/Chidzaa/Documents/Mutambara/Objective 2/SPI1_Umvumvumvu_1989_2018.png",
  width    = 1600,
  height   = 900,
  res      = 150
)

plot(
  x    = data$Date,
  y    = data$SPI_1,
  type = "h",
  col  = ifelse(
    data$SPI_1 < 0,
    "firebrick",
    "steelblue"
  ),
  lwd  = 2,
  xlab = "Date",
  ylab = "SPI-1",
  main = "Gamma-Distribution SPI-1 (Center of Mass)\nUmvumvumvu Watershed (1989-2018)"
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
# 27. EXPORT THE FINAL SPI DATASET
# ============================================================

SPI_output <- data %>%
  dplyr::select(
    Date,
    Year,
    Month,
    Precipitation,
    SPI_1,
    SPI_Classification
  )

write_csv(
  SPI_output,
  "C:/Users/Chidzaa/Documents/Mutambara/Objective 2/SPI_1_Umvumvumvu_1989_2018.csv"
)

cat("\nSPI data exported successfully.\n")


# ============================================================
# 28. EXPORT SUMMARY / SANITY CHECK
# ============================================================

cat("\n--- EXPORT SUMMARY ---\n")
cat("Rows exported:  ", nrow(SPI_output), "\n")
cat("Columns exported:", ncol(SPI_output), "\n")

cat("\nFirst 6 rows:\n")
print(head(SPI_output))

cat("\nSPI-1 summary:\n")
print(summary(SPI_output$SPI_1))

cat("\nClassification counts:\n")
print(table(SPI_output$SPI_Classification, useNA = "ifany"))

cat("\nZero-precipitation months:\n")
print(sum(SPI_output$Precipitation == 0, na.rm = TRUE))

cat("\nAny remaining NA in SPI_1?\n")
print(sum(is.na(SPI_output$SPI_1)))

# ============================================================
# END OF SCRIPT
# ============================================================