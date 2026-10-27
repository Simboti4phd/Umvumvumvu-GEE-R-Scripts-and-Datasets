# ============================================================
# INNOVATIVE TREND ANALYSIS (ITA)
# ERA5-LAND TEMPERATURE: 1989–2018
# ============================================================

# ------------------------------------------------------------
# 1. Install packages if necessary
# ------------------------------------------------------------

install.packages("readr")
install.packages("ggplot2")

# Load packages
library(readr)
library(ggplot2)


# ------------------------------------------------------------
# 2. Import the CSV file
# ------------------------------------------------------------

# Change this to the location of your CSV file
temp <- read_csv("Temp_GEE_1989.csv")

# View the data
head(temp)

# Check the column names
names(temp)


# ------------------------------------------------------------
# 3. Select Year and AnnualTemperature
# ------------------------------------------------------------

# IMPORTANT:
# Replace column names in your CSV file if they are different.

Year <- temp$year
Temperature<- temp$temperature

# ------------------------------------------------------------
# 4. Keep only the period 1989–2018
# ------------------------------------------------------------

data <- temp[temp$year >= 1989 & temp$year <= 2018,]

Year <- data$year
Temperature<- data$temperature


# Check the number of observations
length(Year)
length(Temperature)


# ------------------------------------------------------------
# 5. Divide the data into two halves
# ------------------------------------------------------------

# There are 30 years from 1989 to 2018.
# To obtain equal halves, use 15 years in each half:
#
# First half:  1989–2003
# Second half: 2004–2018
#
# The middle year, 2003, is excluded from the ITA split.

first_half <- Temperature[Year >= 1989 & Year <= 2003]

second_half <- Temperature[Year >= 2004 & Year <= 2018]


# Check the number of observations
length(first_half)
length(second_half)


# ------------------------------------------------------------
# 7. Sort both halves in ascending order
# ------------------------------------------------------------

first_sorted <- sort(first_half)
second_sorted <- sort(second_half)


# ------------------------------------------------------------
# 8. Calculate the ITA slope
# ------------------------------------------------------------

# Linear regression between the sorted first and second halves

ita_model <- lm(second_sorted ~ first_sorted)

# Display model results
summary(ita_model)


# Extract slope
ita_slope <- coef(ita_model)[2]

cat("ITA slope =", ita_slope, "\n")


# ------------------------------------------------------------
# 9. Calculate correlation coefficient
# ------------------------------------------------------------

ita_correlation <- cor(first_sorted, second_sorted)

cat("Correlation coefficient =", ita_correlation, "\n")


# ------------------------------------------------------------
# 10. Create data frame for ITA plot
# ------------------------------------------------------------

ita_data <- data.frame(
  First_Half = first_sorted,
  Second_Half = second_sorted
)


# ------------------------------------------------------------
# 11. Create ITA plot
# ------------------------------------------------------------

ggplot(ita_data, aes(x = First_Half, y = Second_Half)) +
  
  # 1:1 reference line
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 1
  ) +
  
  # ITA regression/trend line
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1
  ) +
  
  # Data points
  geom_point(
    size = 3
  ) +
  
  labs(
    title = "Umvumvumvu Innovative Trend Analysis of Annual Temperature (1989-2018)",
    x = "Temperature - first half (°C)",
    y = "Temperature - second half (°C)"
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    plot.subtitle = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text = element_text(
      face = "bold",
      size = 10
    )
  )


# ------------------------------------------------------------
# 12. Calculate percentage change in the mean
# ------------------------------------------------------------

mean_first <- mean(first_half, na.rm = TRUE)
mean_second <- mean(second_half, na.rm = TRUE)

percentage_change <- 
  ((mean_second - mean_first) / mean_first) * 100

cat("Mean Temperature, first half =", mean_first, "mm\n")
cat("Mean Temperature, second half =", mean_second, "mm\n")
cat("Percentage change =", percentage_change, "%\n")


# ------------------------------------------------------------
# 13. Basic interpretation
# ------------------------------------------------------------

if (ita_slope > 1) {
  cat("ITA indicates an increasing trend in temperature.\n")
} else if (ita_slope < 1) {
  cat("ITA indicates a decreasing trend in temperature.\n")
} else {
  cat("ITA indicates no overall trend in temperature.\n")
}
