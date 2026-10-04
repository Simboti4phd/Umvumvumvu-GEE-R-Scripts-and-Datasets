```r
# ============================================================
# SPI-12 CALCULATION - UMVUMVUMVU WATERSHED
# January 1989 - December 2018
# Gamma distribution + Center of Mass zero handling
# ============================================================

library(readr)
library(dplyr)
library(lubridate)
library(zoo)

# 1. Import data
data <- read_csv("Monthly_PrecipitationCM1.csv",
                 show_col_types = FALSE)

# 2. Prepare data
data <- data %>%
  rename(Precipitation = Precipitation_mm) %>%
  mutate(
    Year = as.numeric(Year),
    Month = as.numeric(Month),
    Precipitation = as.numeric(Precipitation),
    Date = make_date(Year, Month, 1)
  ) %>%
  arrange(Date)

# 3. Create 12-month precipitation accumulation
data$Precip_12M <- rollapply(
  data$Precipitation,
  width = 12,
  FUN = sum,
  align = "right",
  fill = NA
)

# 4. Create 12-month time series
valid_idx <- which(!is.na(data$Precip_12M))

precip_ts_12m <- ts(
  data$Precip_12M[valid_idx],
  start = c(
    year(data$Date[valid_idx[1]]),
    month(data$Date[valid_idx[1]])
  ),
  frequency = 12
)

# 5. Function for SPI calculation
calculate_spi <- function(precip_ts) {
  
  x <- as.numeric(precip_ts)
  months <- cycle(precip_ts)
  spi <- rep(NA_real_, length(x))
  
  for (m in 1:12) {
    
    idx <- which(months == m)
    xm <- x[idx]
    
    q <- mean(xm == 0)
    positive <- xm[xm > 0]
    
    if (length(positive) < 3) next
    
    mu <- mean(positive)
    sigma2 <- var(positive)
    
    if (sigma2 <= 0 || mu <= 0) next
    
    alpha <- mu^2 / sigma2
    beta <- sigma2 / mu
    
    H <- ifelse(
      xm == 0,
      q / 2,
      q + (1 - q) *
        pgamma(xm, shape = alpha, scale = beta)
    )
    
    H <- pmin(pmax(H, 1e-4), 1 - 1e-4)
    
    spi[idx] <- qnorm(H)
  }
  
  spi
}

# 6. Calculate SPI-12
spi12_values <- calculate_spi(precip_ts_12m)

data$SPI_12 <- NA_real_
data$SPI_12[valid_idx] <- spi12_values

# 7. Round SPI values
data$SPI_12 <- round(data$SPI_12, 3)

# 8. Classify SPI-12
data <- data %>%
  mutate(
    SPI_12_Classification = case_when(
      SPI_12 >= 2.00  ~ "Extremely Wet",
      SPI_12 >= 1.50  ~ "Very Wet",
      SPI_12 >= 1.00  ~ "Moderately Wet",
      SPI_12 > -1.00  ~ "Near Normal",
      SPI_12 > -1.50  ~ "Moderately Dry",
      SPI_12 > -2.00  ~ "Severely Dry",
      SPI_12 <= -2.00 ~ "Extremely Dry",
      TRUE ~ NA_character_
    )
  )

# 9. Create final dataset
SPI12_output <- data %>%
  select(
    Date,
    Year,
    Month,
    Precipitation,
    Precip_12M,
    SPI_12,
    SPI_12_Classification
  )

# 10. Display results
print(head(SPI12_output, 12))
print(tail(SPI12_output, 12))
summary(SPI12_output$SPI_12)

# 11. Export SPI-12 results
write_csv(
  SPI12_output,
  "C:/Users/Chidzaa/Documents/Mutambara/Objective 2/SPI_12_Umvumvumvu_1989_2018.csv"
)

cat("\nSPI-12 calculation completed successfully.\n")
```
