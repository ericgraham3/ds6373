library(tswge)

# ============================================================================
# MIDTERM CODE TEMPLATE - TIME SERIES ANALYSIS
# ============================================================================
# Remember: Use = for assignment, not <-
# All code assumes tswge library is loaded

# ============================================================================
# UNIT 1: STATIONARITY
# ============================================================================

## Load and plot data
data(dataset_name)
x = dataset_name
plotts.wge(x)

## Check stationarity conditions
# Condition 1: Constant mean (visual inspection of plot)
# Condition 2: Constant variance (visual inspection of plot)
# Condition 3: ACF depends only on lag h

## Generate ACF
acf(x, lag.max = 30, plot = TRUE)

## Plot realization, ACF, and spectral density together
plotts.sample.wge(x)

## For comparing two sections of data (check condition 3)
x_first_half = x[1:floor(length(x)/2)]
x_second_half = x[(floor(length(x)/2)+1):length(x)]
acf(x_first_half, lag.max = 20)
acf(x_second_half, lag.max = 20)

# ============================================================================
# UNIT 2: FREQUENCY DOMAIN AND SPECTRAL DENSITY
# ============================================================================

## Generate spectral density (Parzen window)
parzen.wge(x)
parzen.wge(x, trunc = 50)  # Adjust truncation point if needed

## Calculate frequency from period (or vice versa)
period = 12
frequency = 1/period
cat("Period:", period, "Frequency:", frequency, "\n")

## Nyquist frequency (maximum detectable frequency)
nyquist_freq = 0.5
cat("Nyquist Frequency:", nyquist_freq, "\n")

## Generate signal with specific frequencies
x_signal = gen.sigplusnoise.wge(n = 200, coef = c(2, 3), 
                                freq = c(0.1, 0.25), 
                                psi = c(0, 0), phi = 0, sn = 123)
plotts.sample.wge(x_signal)

## Apply Butterworth filter (high-pass)
x_highpass = butterworth.wge(x, order = 4, type = "high", cutoff = 0.1)
plotts.sample.wge(x_highpass$x.filt)

## Apply Butterworth filter (low-pass)
x_lowpass = butterworth.wge(x, order = 4, type = "low", cutoff = 0.1)
plotts.sample.wge(x_lowpass$x.filt)

## Apply moving average filter
x_smooth = ma.smooth.wge(x, order = 5, plot = FALSE)
plotts.wge(x_smooth$smooth)

# ============================================================================
# UNIT 3: AR(1) MODELS
# ============================================================================

## Generate AR(1) realization
phi1 = 0.9
x_ar1 = gen.arma.wge(n = 200, phi = phi1, sn = 123)
plotts.sample.wge(x_ar1)

## Check factor table for AR(1)
factor.wge(phi = phi1)

## Plot true ACF and spectral density
plotts.true.wge(phi = phi1)

## General Linear Process (GLP) - get psi weights
psi.weights.wge(phi = phi1, lag.max = 10)

# ============================================================================
# UNIT 4: AR(2) AND AR(P) MODELS
# ============================================================================

## Generate AR(2) realization with real roots
phi_ar2 = c(0.5, 0.3)
x_ar2 = gen.arma.wge(n = 200, phi = phi_ar2, sn = 456)
plotts.sample.wge(x_ar2)

## Generate AR(2) with complex roots
phi_complex = c(1.6, -0.9)
x_complex = gen.arma.wge(n = 200, phi = phi_complex, sn = 789)
plotts.sample.wge(x_complex)

## Check stationarity with factor table
factor.wge(phi = phi_ar2)
factor.wge(phi = phi_complex)
# Stationary if all Abs Recip < 1 (roots outside unit circle)

## Calculate system frequency for complex roots
# For AR(2): f0 = (1/(2*pi)) * acos(phi1 / (2*sqrt(-phi2)))
phi1 = -0.5
phi2 = -0.6
system_freq = (1 / (2 * pi)) * acos(phi1 / (2 * sqrt(abs(phi2))))
cat("System Frequency:", system_freq, "\n")

## Generate AR(p) realization
phi_arp = c(0.6, -0.2, 0.1)
x_arp = gen.arma.wge(n = 200, phi = phi_arp, sn = 111)
plotts.sample.wge(x_arp)
factor.wge(phi = phi_arp)

## Plot true ACF and spectral density for AR(p)
plotts.true.wge(phi = phi_arp)

# ============================================================================
# UNIT 5: MA AND ARMA MODELS
# ============================================================================

## Generate MA(1) realization
theta1 = 0.8
x_ma1 = gen.arma.wge(n = 200, theta = theta1, sn = 222)
plotts.sample.wge(x_ma1)

## Check invertibility for MA
factor.wge(phi = theta1)
# Invertible if all Abs Recip < 1

## Generate MA(q) realization
theta_maq = c(0.8, -0.5)
x_maq = gen.arma.wge(n = 200, theta = theta_maq, sn = 333)
plotts.sample.wge(x_maq)
factor.wge(phi = theta_maq)

## Calculate rho_1 for MA(q) by hand
# rho_1 = (-theta1 + theta1*theta2) / (1 + theta1^2 + theta2^2)
theta1 = 0.8
theta2 = -0.5
rho_1 = (-theta1 + theta1*theta2) / (1 + theta1^2 + theta2^2)
cat("rho_1:", rho_1, "\n")

## Get psi weights (GLP form) for MA
psi.weights.wge(theta = theta_maq, lag.max = 10)

## Generate ARMA realization
phi_arma = c(1.6, -0.9)
theta_arma = 0.9
x_arma = gen.arma.wge(n = 200, phi = phi_arma, theta = theta_arma, 
                      mu = 20, sn = 444)
plotts.sample.wge(x_arma)

## Check stationarity and invertibility
factor.wge(phi = phi_arma)  # AR part - check stationarity
factor.wge(phi = theta_arma)  # MA part - check invertibility

## Plot true ACF and spectral density for ARMA
plotts.true.wge(phi = phi_arma, theta = theta_arma)

## Get psi weights for ARMA
psi.weights.wge(phi = phi_arma, theta = theta_arma, lag.max = 10)

## Use AIC to compare models
aic5.wge(x)
# Returns top 5 models by AIC (lowest is best)

# ============================================================================
# UNIT 6: ARIMA AND SEASONAL MODELS
# ============================================================================

## Generate signal plus noise with linear trend
x_trend = gen.sigplusnoise.wge(n = 200, b0 = 2, b1 = 4, vara = 100, sn = 555)
plotts.wge(x_trend)

## Generate signal plus noise with cyclic component
x_cyclic = gen.sigplusnoise.wge(n = 200, coef = c(2, 1), 
                                freq = c(0.05, 0.28), 
                                psi = c(1.2, 2), phi = 0.9, sn = 666)
plotts.sample.wge(x_cyclic)

## Generate ARIMA(0,1,0) - random walk
x_arima010 = gen.arima.wge(n = 200, d = 1, sn = 777)
plotts.sample.wge(x_arima010)

## Generate ARIMA(p,d,q)
phi_arima = c(0.7)
theta_arima = c(0.4)
x_arima = gen.arima.wge(n = 200, phi = phi_arima, theta = theta_arima, 
                        d = 1, sn = 888)
plotts.sample.wge(x_arima)

## Apply first difference using artrans.wge
x_diff = artrans.wge(x, phi.tr = 1)
plotts.sample.wge(x_diff)

## Apply second difference
x_diff2 = artrans.wge(x, phi.tr = c(2, -1))
plotts.sample.wge(x_diff2)

## Generate seasonal model (1-B^s)
# For s=12 (monthly data with annual seasonality)
x_seasonal = gen.arima.wge(n = 200, s = 12, sn = 999)
plotts.sample.wge(x_seasonal)

## Apply seasonal difference
x_seasonal_diff = artrans.wge(x, phi.tr = c(rep(0, 11), 1))
plotts.sample.wge(x_seasonal_diff)

## Factor table for seasonal factor
factor.wge(phi = c(rep(0, 6), 1))  # s=7 example

## Multiply factor models
parms = mult.wge(c(0.975), c(0.2, 0.45), c(-0.53))
parms$model.coef
x_mult = gen.arma.wge(n = 200, phi = parms$model.coef, vara = 1, sn = 101)
plotts.sample.wge(x_mult)

## Airline model: (1-B)(1-B^12)
phi_airline = c(rep(0, 10), 1, 0, -1)
x_airline = gen.arima.wge(n = 200, phi = phi_airline, d = 0, sn = 202)

# ============================================================================
# UNIT 7: FORECASTING
# ============================================================================

## Forecast from ARMA model
phi_forecast = c(0.7)
theta_forecast = c(0.4)
x_data = gen.arma.wge(n = 100, phi = phi_forecast, theta = theta_forecast, 
                      sn = 303)

# Generate forecasts
forecasts = fore.arma.wge(x_data, phi = phi_forecast, theta = theta_forecast, 
                          n.ahead = 10, lastn = TRUE, limits = TRUE)

## Forecast from ARIMA model
x_arima_data = gen.arima.wge(n = 100, phi = c(0.9), d = 1, sn = 404)
forecasts_arima = fore.arima.wge(x_arima_data, phi = c(0.9), d = 1, 
                                 n.ahead = 10, lastn = TRUE, limits = TRUE)

## Forecast from seasonal/airline model
x_seasonal_data = gen.arima.wge(n = 100, s = 12, sn = 505)
forecasts_seasonal = fore.arima.wge(x_seasonal_data, s = 12, 
                                    n.ahead = 12, lastn = TRUE, limits = TRUE)

## Calculate ASE (Average Squared Error) manually
actual = x_data[91:100]
predicted = forecasts$f[1:10]
ase = mean((actual - predicted)^2)
cat("ASE:", ase, "\n")

# ============================================================================
# MODEL IDENTIFICATION AND COMPARISON
# ============================================================================

## Match realizations to models
# Generate multiple models and compare visually
model_a = gen.arma.wge(n = 100, phi = c(0.9), sn = 21)
model_b = gen.arma.wge(n = 100, phi = c(1.6, -0.9), sn = 229)
model_c = gen.arma.wge(n = 100, phi = c(0.28, 0.134, -0.8924), sn = 26)

plotts.sample.wge(model_a)
plotts.sample.wge(model_b)
plotts.sample.wge(model_c)

## Compare with factor tables
factor.wge(phi = c(0.9))
factor.wge(phi = c(1.6, -0.9))
factor.wge(phi = c(0.28, 0.134, -0.8924))

## Model selection with AIC
aic5.wge(x, p = 0:10, q = 0:5)
# Returns top 5 models with lowest AIC
# BIC penalizes complexity more than AIC

# ============================================================================
# COMMON PATTERNS FOR MIDTERM
# ============================================================================

## Identifying model type from ACF:
# - Exponential decay: AR model
# - Cuts off after lag q: MA(q) model
# - Damped sine wave: Complex roots (AR with complex roots)
# - Slow decay: Wandering behavior, may need differencing (ARIMA)
# - Alternating signs: Negative root or oscillatory behavior
# - Repeating pattern at lag s: Seasonal component

## Identifying model type from Spectral Density:
# - Peak at f=0: Wandering/trending behavior
# - Peak at f=0.5: Period-2 oscillations (alternating pattern)
# - Peak at specific f: Cyclic behavior with period = 1/f
# - Multiple peaks: Multiple frequencies present
# - Flat: White noise

## Identifying model type from Realization:
# - Smooth wandering: High positive phi or (1-B) factor
# - Rapid oscillations: Negative root or high frequency
# - Clear trend: Signal plus noise or (1-B) factor
# - Repeating pattern: Seasonal component or specific frequency

## Quick reference for system frequencies:
# f = 0: Wandering/trending (positive real root near 1)
# f = 0.5: Period-2 oscillations (negative real root near -1)
# f = other: Complex roots, period = 1/f

# ============================================================================
# USEFUL DEBUGGING AND CHECKING
# ============================================================================

## Check length of data
length(x)

## Check structure of data
str(x)

## View first few observations
head(x, 10)

## Summary statistics
summary(x)
mean(x)
var(x)

## Check for missing values
sum(is.na(x))

# ============================================================================
# END OF TEMPLATE
# ============================================================================