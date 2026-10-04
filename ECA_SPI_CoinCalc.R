# Load the package
library(CoinCalc)

# Check the package version
packageVersion("CoinCalc")  # Should return ‘2.0’

# Read the CSV file
data <- read.csv("Min_Flow_SPI.csv")

# Check the structure
str(data)

# Adjust date format (adjust the format string as needed)
data$Period <- as.Date(data$Period, format = "%Y-%m-%d")

#Binarize Events
# Event B: River intermittence onset
data$intermittence_event <- ifelse(data$Min_Flow <= 0.01, 1, 0)

# Event A: Meteorological drought (do this for each SPI timescale)
data$spi1_event <- ifelse(data$SPI_1 <= -1.0, 1, 0)
data$spi3_event <- ifelse(data$SPI_3 <= -1.0, 1, 0)
data$spi6_event <- ifelse(data$SPI_6 <= -1.0, 1, 0)
data$spi12_event <- ifelse(data$SPI_12 <= -1.0, 1, 0)

eca_spi1_test <- CC.eca.ts(
  seriesA = data$spi1_event,
  seriesB = data$intermittence_event,
  delT = 1,
  tau = 1,
  sigtest = "poisson"
)

eca_spi3_test <- CC.eca.ts(
  seriesA = data$spi3_event,
  seriesB = data$intermittence_event,
  delT = 1,
  tau = 1,
  sigtest = "poisson"
)

eca_spi6_test <- CC.eca.ts(
  seriesA = data$spi6_event,
  seriesB = data$intermittence_event,
  delT = 1,
  tau = 1,
  sigtest = "poisson"
)
eca_spi12_test <- CC.eca.ts(
  seriesA = data$spi12_event,
  seriesB = data$intermittence_event,
  delT = 1,
  tau = 1,
  sigtest = "poisson"
)


