install.packages("lubridate")
library("data.table")
library("lubridate")
library("dplyr")
library("forecast")
library("ggplot2")
library("tseries")


####
# Read in the data

visits_number <- read.csv("Visits_Number.csv", 
                          na.strings = "", 
                          stringsAsFactors = FALSE)

# View(visits_number)
altnagelvin   <- visits_number[, c("Year", "MthAndYrCode", "New_Hospital", "Type", "Total_sum")] %>% 
  filter((New_Hospital == "Altnagelvin Area") & (Type == "Type 1"))

altnagelvin$date <- parse_date_time2(x = altnagelvin$MthAndYrCode, orders = "ym")
names(altnagelvin) <- tolower(names(altnagelvin))


# Convert to a ts class
altnagelvin$total_sum_ts <- ts(altnagelvin$total_sum, 
                               start = year(altnagelvin$date[1]), 
                               frequency = 12)
tot_visits_ts <- altnagelvin$total_sum_ts  # independent ts Easier to deal with


# Exploratory analysis 
options(repr.plot.width=5, repr.plot.height=4.5)
ggplot(data = altnagelvin, aes(x = date, y = total_sum))+
  geom_line(color = "#00AFBB", size = 2)

# Look at the Acf
ggAcf(tot_visits_ts, lag = 30) # Slowly decaying acf could be a sign of non-stationarity

# The series plot and a slowly decaying and cyclic acf suggests that not only is the data 
# likely to have seasonality but may also be non-stationary

plot(tot_visits_ts)
abline(reg=lm(tot_visits_ts~time(tot_visits_ts)))

ggseasonplot(tot_visits_ts, year.labels=TRUE, year.labels.left=TRUE, polar = TRUE) +
  ylab("Seasonal plot of total visits") +
  ggtitle("")

# The boxplot gives us an idea of seasonal effect
boxplot(tot_visits_ts~cycle(tot_visits_ts))

ggsubseriesplot(tot_visits_ts) +
  ylab("Visits") +
  ggtitle("Seasonal plot of total visits")

# Clearly there is both a trend and seasonality
# But it does look proportional to the series so an additive model may suffice
# From the diagrams we can infer that
# 1. that the trend has been dealing with more and more patients on a yearly basis
# 2. The variance and the mean is much higher in the winter months

# Lets try a augmented Dickey Fuller test and a kpss test
adf.test(tot_visits_ts)
# H0: unit root
# Ha: stationary

kpss.test(tot_visits_ts)
# H0: stationary
# Ha: has unit root

# From ADF test reject H0 at 5% significance level, the series are likely stationary. 
# From the KPSS test reject H0 under 5% significance level, the series are likely to 
# have a unit root. Two tests condradict each other. 
#  One of the reasons is explained here 
# https://stats.stackexchange.com/questions/239360/contradictory-results-of-adf-and-kpss-unit-root-tests 

# Let's apply the tests to a seasonally adjusted data and see if the results will change.
total_visits_decomposed <- tot_visits_ts %>% mstl() 
total_visits_decomposed %>% autoplot()

# Adjust the data
# u_t = y_t - S_t where y_t is the original data and S_t is the seasonal component
total_visits_seasonally_adjusted <- total_visits_decomposed[, "Data"] - total_visits_decomposed[, "Seasonal12"]
total_visits_seasonally_adjusted %>% autoplot()

# Test the adjusted data 
adf.test(total_visits_seasonally_adjusted)
# H0: unit root
# Ha: stationary

kpss.test(total_visits_seasonally_adjusted)
# H0: stationary
# Ha: has unit root

# Now both test agree that there is a unit root under 5% significance level

### Let's find out the order of non-seasonal differencing
ndiffs(total_visits_seasonally_adjusted)


tot_visits_diff <- diff(tot_visits_ts)

tot_visits_diff %>% autoplot()

ggAcf(tot_visits_diff, lag.max = 60)

ggPacf(tot_visits_diff, lag.max = 60)

# the acf shows some slow decay at seasonal lags.
# We can take a seasonal difference at lag 12 of the
# first differenced data to stabilize it
tot_visits_diff_seasonal_diff <- tot_visits_diff %>% diff(lag = 12)
tot_visits_diff_seasonal_diff %>% autoplot()

ggAcf(tot_visits_diff_seasonal_diff, lag.max = 60)

# The significant negative correlation at lag 1 suggests a non-seasonal MA component of order 1. 
# This also looks similar to the MA signature described here under Rule 7 
# https://people.duke.edu/~rnau/411arim3.htm#signatures. 
# A significant negative correlation at lag 12 suggests a seasonal MA component of order 1.

ggPacf(tot_visits_diff_seasonal_diff, lag.max = 60)

# The PACF plot above suggests a non-seasonal AR component of order 2 due to significant correlation 
# at lags 1 and 2. Futhermore, significant correlation at lags 12 and 24 suggests addition of a 
# seasonal AR term of order 2.

######################################################################################################
# Fit the models 

# The analysis suggests trying two models
# ARIMA(0,1,1)(0,1,1)12
# ARIMA(2,1,0)(2,1,0)12

# Now lets fit the arima model 

######## First Test ###########
# Applied to the original data
sarima_1 <- tot_visits_ts %>% Arima(order=c(0, 1, 1), seasonal=c(0, 1, 1))

sarima_1_residuals <- sarima_1 %>% residuals()
qqnorm(sarima_1_residuals)
qqline(sarima_1_residuals, col = "red")

sarima_1 %>% checkresiduals # Includes a Ljung-Box test

shapiro.test(sarima_1_residuals)
# There is no significant correlation left in the residuals. 
# From the Ljung-Box test we fail to reject H0 under 5% significance level meaning it is likely 
# that the residuals are independently distributed. 
# Similar conclusion follows from the ACF plot; there is no autocorrelation left in the series. 
# From Shapiro-Wilk test we reject H0 under 5% significance level, the residuals are likely not 
# normally distributed. 
# From the residual plot it is hard to tell if the variance has some patterns or not, although it 
# looks like the variance gets slightly larger going towards 2018

######## Second Test ###########
# Applied to the original data
sarima_2 <- tot_visits_ts %>% Arima(order=c(2, 1, 0), seasonal=c(2, 1, 0))

sarima_2_residuals <- sarima_2 %>% residuals()
qqnorm(sarima_2_residuals)
qqline(sarima_2_residuals, col = "red")

sarima_2 %>% checkresiduals  # Includes a Ljung-Box test

shapiro.test(sarima_2_residuals)
# There is no significant correlation left in the residuals. 
# From the Ljung-Box test we fail to reject H0 under 5% significance level meaning it is likely 
# that the residuals are independently distributed. 
# Similar conclusion follows from the ACF plot; there is no autocorrelation left in the series. 
# From Shapiro-Wilk test we reject H0 under 5% significance level, the residuals are likely not 
# normally distributed. 
# From the residual plot it is hard to tell if the variance has some patterns or not, although it 
# looks like the variance gets slightly larger going towards 2018

# Selecting the best fitting model
summary(sarima_1)
summary(sarima_2)
# The mean absolute percentage error (MAPE)
# measures prediction of accuracy
# Here the MAPE is 2.6% of the number of patient visits
# so this is the forecast accuracy of the error

#  Use auto arima to search and validate manual selection
best_sarima <- auto.arima(y=tot_visits_ts, 
                          d = 1,
                          D = 1,
                          max.p =3,
                          max.q =3,
                          max.P =3,
                          max.Q =3,
                          max.d =1,
                          max.D =1,
                          stationary=FALSE,
                          seasonal=TRUE,
                          stepwise=FALSE)
best_sarima

######################################################################
### Forecast

sarima_1 %>% forecast(h=12)
sarima_1 %>% forecast(h=12) %>% autoplot()


