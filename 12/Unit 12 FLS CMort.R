### Unit 12 FLS

#Question 1

#MLR WITH CORRELATED ERRORS

# ARIMA 1 MLR with Cor Errors (no lag, no seasonl categorical variable)

library(tidyverse)
library(GGally)
library(astsa)
library(tswge)

CM = read.csv(file.choose(),header = TRUE)

head(CM)
ggpairs(CM[2:4]) #matrix of scatter plots


# Start all Cardiac Mortality, Temp and Particulate Data

#Find ASE  Need to forecast last 52 of known series.  
CMsmall = CM[1:456,]
ksfit = lm(cmort~temp + part + Week, data = CMsmall)
phi = aic.wge(ksfit$residuals)
phi$p
fit = arima(CMsmall$cmort,order = c(phi$p,0,0), xreg = cbind(CMsmall$temp, CMsmall$part, CMsmall$Week))

AIC(fit) #2849.62

last52 = data.frame(temp = CM$temp[457:508], part = CM$part[457:508], CMWeek = seq(457,508,1))
#get predictions
predsCMort = predict(fit,newxreg = last52)

ASE = mean((CM$cmort[457:508] - predsCMort$pred)^2)
ASE # 66.55999

plot(seq(1,508,1), CM$cmort, type = "l",xlim = c(0,528), ylab = "Cardiac Mortality", main = "Last 30 Week Cardiac Mortality Forecast")
lines(seq(457,508,1), predsCMort$pred, type = "l", col = "red")

# ... is there an issue?



# Naive with temp and part and Week

# Need for forecast temp and part since we would not know that information... we would know Week

CMsmall = CM[1:456,]

#forecast Particles
plotts.sample.wge(CMsmall$part) #freq near .0192 (annual)
acf(CMsmall$part,lag.max = 100)
CM_52 = artrans.wge(CMsmall$part, c(rep(0,51),1))
plotts.sample.wge(CM_52) #looks like some low freq?
aic5.wge(CM_52) #picks ARMA(2,1) assume stationary
aic5.wge(CM_52,type = "bic") #picks ARMA(0,0) 
ljung.wge(CM_52)$pval #FTR Ho
ljung.wge(CM_52, K = 48)$pval #FTR Ho
#Going with white noise despite peak at 0 in Spec D. 
#CM_52_AR2_MA1 = artrans.wge(CM_52,est$phi)
predsPart = fore.aruma.wge(CMsmall$part,s = 52, n.ahead = 52)
plot(predsPart$f, type = "l")
plot(seq(1,456,1), CMsmall$part, type = "l",xlim = c(0,508), ylab = "Particles", main = "52 Week Particulate Forecast")
lines(seq(457,508,1), predsPart$f, type = "l", col = "red")


#forecast Temp
plotts.sample.wge(CMsmall$temp) #freq near .0192 (annual)
CM_52 = artrans.wge(CMsmall$temp, c(rep(0,51),1))
plotts.sample.wge(CM_52) #looks like some low freq?
aic5.wge(CM_52) #picks ARMA(0,0)
aic5.wge(CM_52,type = "bic") #picks ARMA(0,0) 
ljung.wge(CM_52)$pval
ljung.wge(CM_52, K = 48)$pval #barely rejects
acf(CM_52,lag.max = 48) # acf looks consistent with white noise
predsTemp = fore.aruma.wge(CMsmall$temp,s = 52, n.ahead = 52)
plot(predsTemp$f, type = "l")
plot(seq(1,456,1), CMsmall$temp, type = "l",xlim = c(0,508), ylab = "Temperature", main = "52 Week Temperature Forecast")
lines(seq(457,508,1), predsTemp$f, type = "l", col = "red")




#Find ASE  Need to forecast last 52 of known series.  
CMsmall = CM[1:456,]
ksfit = lm(cmort~temp + part + Week, data = CMsmall)
phi = aic.wge(ksfit$residuals)
phi$p
fit = arima(CMsmall$cmort,order = c(phi$p,0,0), xreg = cbind(CMsmall$temp, CMsmall$part, CMsmall$Week))

AIC(fit) #2849.62

last52 = data.frame(temp = predsTemp$f, part = predsPart$f, CMWeek = seq(457,508,1))
#get predictions
predsCMort = predict(fit,newxreg = last52)

ASE = mean((CM$cmort[457:508] - predsCMort$pred)^2)
ASE # 77.06643

plot(seq(1,508,1), CM$cmort, type = "l",xlim = c(0,528), ylab = "Cardiac Mortality", main = "Last 30 Week Cardiac Mortality Forecast")
lines(seq(457,508,1), predsCMort$pred, type = "l", col = "red")




# Forecast Ahead 52

#forecast Particles
plotts.sample.wge(CM$part) #freq near .0192 (annual)
CM_52 = artrans.wge(CM$part, c(rep(0,51),1))
plotts.sample.wge(CM_52) #looks like some low freq?
aic5.wge(CM_52) #picks ARMA(2,1) assume stationary
aic5.wge(CM_52,type = "bic") #picks ARMA(0,0) 
ljung.wge(CM_52)$pval #FTR Ho
ljung.wge(CM_52, K = 48)$pval #FTR Ho
#Going with white noise despite peak at 0 in Spec D. 
#CM_52_AR2_MA1 = artrans.wge(CM_52,est$phi)
predsPart = fore.aruma.wge(CM$part,s = 52, n.ahead = 52)
plot(predsPart$f, type = "l")
plot(seq(1,508,1), CM$part, type = "l",xlim = c(0,560), ylab = "Temperature", main = "52 Week Particulate Forecast")
lines(seq(509,560,1), predsPart$f, type = "l", col = "red")


#forecast Temp
plotts.sample.wge(CM$temp) #freq near .0192 (annual)
CM_52 = artrans.wge(CM$temp, c(rep(0,51),1))
plotts.sample.wge(CM_52) #looks like some low freq?
aic5.wge(CM_52) #picks ARMA(0,0)
aic5.wge(CM_52,type = "bic") #picks ARMA(0,0) 
ljung.wge(CM_52)$pval
ljung.wge(CM_52, K = 48)$pval #barely rejects
acf(CM_52,lag.max = 48) # acf looks consistent with white noise
predsTemp = fore.aruma.wge(CM$temp,s = 52, n.ahead = 52)
plot(predsTemp$f, type = "l")
plot(seq(1,508,1), CM$temp, type = "l",xlim = c(0,560), ylab = "Temperature", main = "20 Week Temperature Forecast")
lines(seq(509,560,1), predsTemp$f, type = "l", col = "red")


# Model cmort based on predicted part and temp using MLR with Cor Erros
#assuming data is loaded in dataframe CM
ksfit = lm(cmort~temp+part+Week, data = CM)
summary(ksfit)
plotts.sample.wge(ksfit$residuals)
phi = aic.wge(ksfit$residuals) #AR(2)
phi

fit = arima(CM$cmort,order = c(phi$p,0,phi$q), xreg = cbind(CM$temp, CM$part, CM$Week))

AIC(fit) #3167.096


# Check for whiteness of residuals
plotts.sample.wge(fit$residuals)

acf(fit$residuals)
ljung.wge(fit$residuals, p = 2) # pval = .048
ljung.wge(fit$residuals, K = 48, p = 2) # pval = .0004

#load the forecasted Part and Temp in a data frame
next52 = data.frame(temp = predsTemp$f, part = predsPart$f, Week = seq(509,560,1))
#get predictions
predsCMort = predict(fit,newxreg = next52)
#plot next 20 cmort wrt time
plot(seq(1,508,1), CM$cmort, type = "l",xlim = c(0,560), ylab = "Cardiac Mortality", main = "20 Week Cardiac Mortality Forecast")
lines(seq(509,560,1), predsCMort$pred, type = "l", col = "red")












#Start With Lag Temp and Part


#Find ASE  Need to forecast last 56 of known series.  With temp lag 1 and part lag 7
ccf(CM$temp,CM$cmort) #lag 1
ccf(CM$part,CM$cmort) #lag 7

CMsmall = CM[1:456,]
CMsmall$temp_1 = dplyr::lag(CMsmall$temp,1)
CM$temp_1 = dplyr::lag(CM$temp,1)

CMsmall = CM[1:456,]
CMsmall$part_7 = dplyr::lag(CMsmall$part,7)
CM$part_7 = dplyr::lag(CM$part,7)

ggpairs(CMsmall[2:6]) #matrix of scatter plots

ksfit = lm(cmort~temp_1 + part_7+Week, data = CMsmall)
summary(ksfit)
aic5.wge(ksfit$residuals)
phi = aic.wge(ksfit$residuals)
phi$p
phi$q
fit = arima(CMsmall$cmort,order = c(4,0,0), xreg = cbind(CMsmall$temp_1, CMsmall$part_7, CMsmall$Week))

AIC(fit)  #2821.321


lag1PredTemp = dplyr::lag(predsTemp$f,1)
lag7PredPart = dplyr::lag(predsPart$f,7)

Lag1Temp = c(CMsmall$temp[456],lag1PredTemp[2:52])
Lag7Part = c(CMsmall$part[450:456],lag7PredPart[8:52])


last52 = data.frame(temp_1 = Lag1Temp, part_7 = Lag7Part, Week = seq(457,508,1))
#get predictions
predsCMort = predict(fit,newxreg = last52)

plot(seq(1,508,1), CM$cmort, type = "l",xlim = c(0,528), ylab = "Cardiac Mortality", main = "Last 52 Week Cardiac Mortality Forecast")
lines(seq(457,508,1), predsCMort$pred, type = "l", col = "red")


ASE = mean((CM$cmort[457:508] - predsCMort$pred)^2)
ASE #38.93088







# Forecast Ahead 52

#forecast Particles
plotts.sample.wge(CM$part) #freq near .0192 (annual)
CM_52 = artrans.wge(CM$part, c(rep(0,51),1))
plotts.sample.wge(CM_52) #looks like some low freq?
aic5.wge(CM_52) #picks ARMA(2,1) assume stationary
aic5.wge(CM_52,type = "bic") #picks ARMA(0,0) 
ljung.wge(CM_52)$pval #FTR Ho
ljung.wge(CM_52, K = 48)$pval #FTR Ho
#Going with white noise despite peak at 0 in Spec D. 
#CM_52_AR2_MA1 = artrans.wge(CM_52,est$phi)
predsPart = fore.aruma.wge(CM$part,s = 52, n.ahead = 52)
plot(predsPart$f, type = "l")
plot(seq(1,508,1), CM$part, type = "l",xlim = c(0,560), ylab = "Temperature", main = "52 Week Particulate Forecast")
lines(seq(509,560,1), predsPart$f, type = "l", col = "red")


#forecast Temp
plotts.sample.wge(CM$temp) #freq near .0192 (annual)
CM_52 = artrans.wge(CM$temp, c(rep(0,51),1))
plotts.sample.wge(CM_52) #looks like some low freq?
aic5.wge(CM_52) #picks ARMA(0,0)
aic5.wge(CM_52,type = "bic") #picks ARMA(0,0) 
ljung.wge(CM_52)$pval
ljung.wge(CM_52, K = 48)$pval #barely rejects
acf(CM_52,lag.max = 48) # acf looks consistent with white noise
predsTemp = fore.aruma.wge(CM$temp,s = 52, n.ahead = 52)
plot(predsTemp$f, type = "l")
plot(seq(1,508,1), CM$temp, type = "l",xlim = c(0,560), ylab = "Temperature", main = "20 Week Temperature Forecast")
lines(seq(509,560,1), predsTemp$f, type = "l", col = "red")


# Model cmort based on predicted part and temp using MLR with Cor Erros
#assuming data is loaded in dataframe CM
ksfit = lm(cmort~temp_1+part_7+Week, data = CM)
summary(ksfit)
plotts.sample.wge(ksfit$residuals)
phi = aic.wge(ksfit$residuals) #AR(2)
phi

fit = arima(CM$cmort,order = c(phi$p,0,phi$q), xreg = cbind(CM$temp_1, CM$part_7, CM$Week))

AIC(fit) #3095.895

  





 #Adding Seasonal Dummy and lag Temp and part


#Find ASE  Need to forecast last 52 of known series.  With temp lag 1 and part lag 7

CMsmall = CM[1:456,]
CMsmall$temp_1 = dplyr::lag(CMsmall$temp,1)
CM$temp_1 = dplyr::lag(CM$temp,1)

CMsmall = CM[1:456,]
CMsmall$part_7 = dplyr::lag(CMsmall$part,7)
CM$part_7 = dplyr::lag(CM$part,7)

#ggpairs(CMsmall[2:6]) #matrix of scatter plots

CMsmall$FWeek = as.factor(CMsmall$Week%%52)
ksfit = lm(cmort~temp_1+part_7+Week+FWeek, data = CMsmall)

summary(ksfit)
aic5.wge(ksfit$residuals)
phi = aic.wge(ksfit$residuals)
phi$p
phi$q


lag1PredTemp = dplyr::lag(predsTemp$f,1)
lag7PredPart = dplyr::lag(predsPart$f,7)

Lag1Temp = c(CMsmall$temp[456],lag1PredTemp[2:52])
Lag7Part = c(CMsmall$part[450:456],lag7PredPart[8:52])


#load the forecasted Part and Temp in a data frame
next52 = data.frame(temp_1 = Lag1Temp, part_7 = Lag7Part, Week = seq(457,508,1),FWeek = as.factor(seq(457,508,1)%%52))

#predict residuals manually
plotts.sample.wge(ksfit$residuals)
phi = aic.wge(ksfit$residuals)
resids = fore.arma.wge(ksfit$residuals,phi = phi$phi,n.ahead = 52)
#predict trend manually
preds = predict(ksfit, newdata = next52)

#Combine forecast trend and residual
predsFinal = preds + resids$f

predsCMort = predsFinal

plot(seq(1,508,1), CM$cmort, type = "l",xlim = c(0,528), ylab = "Cardiac Mortality", main = "Last 52 Week Cardiac Mortality Forecast")
lines(seq(457,508,1), predsCMort, type = "l", col = "red")


ASE = mean((CM$cmort[457:508] - predsCMort)^2)
ASE #39.4


### I didn't forecast ahead 52 with the seasonal ... higher ASE







#VAR
# With VAR there is no need to forecast the temp and particles ahead since they are responses now
# and will be forecast as part of the model predictions. 





#Find ASE using last 52
library(vars)
CMsmall = CM[1:456,]  # 456 = 508-52
VARselect(cbind(CMsmall$cmort[2:456], CMsmall$part[2:456], CMsmall$temp[2:456]),lag.max = 10, type = "both")

CMortVAR = VAR(cbind(CMsmall$cmort[2:456], CMsmall$part[2:456], CMsmall$temp[2:456]), type = "both",lag.max = 10)
preds=predict(CMortVAR,n.ahead=52)

#Plot
plot(seq(1,508,1), CM$cmort, type = "l",xlim = c(0,508), ylab = "Cardiac Mortality", main = "52 Week Cardiac Mortality Forecast")
lines(seq(457,508,1), preds$fcst$y1[,1], type = "l", col = "red")


ASE = mean((CM$cmort[457:508] - preds$fcst$y1[,1])^2)
ASE #35.07703



### Forecast Ahead 52


VARselect(cbind(CM$cmort, CM$part, CM$temp),lag.max = 10, type = "both")

CMortVAR = VAR(cbind(CM$cmort, CM$part, CM$temp), type = "both",p = 2)
preds=predict(CMortVAR,n.ahead=152)

#Plot
plot(seq(1,508,1), CM$cmort, type = "l",xlim = c(0,660), ylab = "Cardiac Mortality", main = "52 Week Cardiac Mortality Forecast")
lines(seq(509,660,1), preds$fcst$y1[,1], type = "l", col = "red")












#Find ASE using last 52 with Seasonal Dummy Variables
library(vars)
CMsmall = CM[1:456,]  # 456 = 508-52
VARselect(cbind(CMsmall$cmort[2:456], CMsmall$part[2:456], CMsmall$temp[2:456]),lag.max = 10, season = 52, type = "both")

CMortVAR = VAR(cbind(CMsmall$cmort[2:456], CMsmall$part[2:456], CMsmall$temp[2:456]), season = 52, type = "both",p = 2)
preds=predict(CMortVAR,n.ahead=52)

#Plot
plot(seq(1,508,1), CM$cmort, type = "l",xlim = c(0,508), ylab = "Cardiac Mortality", main = "52 Week Cardiac Mortality Forecast")
lines(seq(457,508,1), preds$fcst$y1[,1], type = "l", col = "red")


ASE = mean((CM$cmort[457:508] - preds$fcst$y1[,1])^2)
ASE #40.67298



### Forecast Ahead 52


VARselect(cbind(CM$cmort, CM$part, CM$temp),lag.max = 10, season = 52, type = "both")

CMortVAR = VAR(cbind(CM$cmort, CM$part, CM$temp),season = 52, type = "both",p = 2)
preds=predict(CMortVAR,n.ahead=52)

#Plot
plot(seq(1,508,1), CM$cmort, type = "l",xlim = c(0,560), ylab = "Cardiac Mortality", main = "52 Week Cardiac Mortality Forecast")
lines(seq(509,560,1), preds$fcst$y1[,1], type = "l", col = "red")






