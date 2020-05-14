#####################################################################################
###################################  P R O J E C T  #################################
#####################################################################################

#Load packages
library(fBasics)
library(tseries)
library(astsa)
library(fUnitRoots)
library(TSA)
library(lmtest)
library(fpp2)
library(forecast)
library(GGally)

#load parking data
park=read.csv('2018.06.01-2019.04.13_ParkingDataFix.csv',header=T)
head(park)
str(park)

###Look at Revenue Variables 1 at a time
#set each revenue variable
xT=(park$T2PP_Revenue)
xR=(park$Reservations_Revenue)
xH=(park$Hub_Revenue)
xV=(park$Valet_Rev)

xTT=(park$T2PP_Trans)
xRR=(park$Reservations_Bookings)
xHH=(park$Hub_Trans)
xVV=(park$Valet_Trans)

#plot time series with all revenue
#Convert Date from Factor to Date
park$Date=as.Date(park$Date,format="%m/%d/%Y")
str(park)

plot(xT~park$Date,type='l',col='skyblue2',main='Parking Revenue for June 1, 2018 - April 13, 2019',ylim=c(0,250000),xlab='Date',ylab='$')
points(xR~park$Date,type='l',col='seagreen2')
points(xH~park$Date,type='l',col='cyan4')
points(xV~park$Date,type='l',col='midnightblue')
legend('topright',legend=c("T2 Parking Plaza", "Reservations", "HUB (5 Lots)","Valet"),col=c('skyblue2','seagreen2','cyan4','midnightblue'),lty=c(1,1,1))

plot(xTT~park$Date,type='l',col='skyblue2',main='Parking Transactions for June 1, 2018 - April 13, 2019',ylim=c(0,8000),xlab='Date',ylab='# of Transactions')
points(xRR~park$Date,type='l',col='seagreen2')
points(xHH~park$Date,type='l',col='cyan4')
points(xVV~park$Date,type='l',col='midnightblue')
legend('topright',legend=c("T2 Parking Plaza", "Reservations", "HUB (5 Lots)","Valet"),col=c('skyblue2','seagreen2','cyan4','midnightblue'),lty=c(1,1,1))


#Compute Mean
mean(xT)
mean(xR)
mean(xH)
mean(xV)

# Compute median and range of all 3
median(xT)
median(xR)
median(xH)
median(xV)
range(xT)
range(xR)
range(xH)
range(xV)

# Compute sample variation and sample std dev of all 3
var(xT)
var(xR)
var(xH)
var(xV)
sd(xT)
sd(xR)
sd(xH)
sd(xV)
#Plot histograms of all 3
par(mfrow=c(1,1))
hist(xT, main = paste("T2PP Revenue"))
hist(xR, main = paste("Reservations Revenue"))
hist(xH, main = paste("T1,T2W,LTL1,ECON,Lot7 Revenue"))
hist(xV, main = paste("Valet Revenue"))

GGally::ggpairs(as.data.frame(park))

############### T2PP Revenue ####################
par(mfrow=c(1,3))
plot(xT~park$Date,type='l',main='T2PP Revenue',xlab='Date',ylab='$')
acf(xT,lag.max=60,main='T2PP Revenue ACF')
pacf(xT,lag.max=60,main='T2PP Revenue PACF')
adf.test(xT)

#The variace seems pretty stable. Will not take the log.
#We will take the difference of the T2PP data, replot, and conduct DF test
dxT=diff(xT)
plot(dxT,type='l', main='Diff T2PP Revenue')
acf(dxT,lag.max=60,main='Diff T2PP Revenue ACF')
pacf(dxT,lag.max=60,main='Diff T2PP Revenue PACF')
adf.test(dxT)

#Mean looks much more stabilized.
#P-value = 0.01; data is stationary
#We can tell that the data shows strong seasonality at period 7.
#Autocorrelation is strong at lags 7, 14, 21, 28, etc.
#PACF is significant around lag 7.
#Because ACF is decaying very slowly, we will take seasonal diff first
outT1=arima(dxT,order=c(0,0,0),seasonal=list(order=c(0,1,0),period=7))
outT1
par(mfrow=c(1,2))
acf(outT1$residuals,lag.max=60)
pacf(outT1$residuals,lag.max=60)

#aic = 6447.63
#For seasonal, ACF shows something at Lag 7, but not 14, 21, etc.
#PACF, shows some significance at lag 7.
#Will consider SMA1, rerun, and evaluate.
outT2=arima(xT,order=c(0,1,0),seasonal=list(order=c(0,1,1),period=7))
outT2
coeftest(outT2)
acf(outT2$residuals,lag.max=60)
pacf(outT2$residuals,lag.max=60)

#All are significant.
#aic = 6278.44
#ACF and PACF are not significant at seasonal lags.
#However, still seeing significant after period 7.
#Will increase D by 1, rerun, and evaluate.
outT3=arima(xT,order=c(0,1,0),seasonal=list(order=c(0,2,1),period=7))
outT3
coeftest(outT3)
par(mfrow=c(1,2))
acf(outT3$residuals,lag.max=60)
pacf(outT3$residuals,lag.max=60)

#SMA1 still significant.
#aic = 6344.92
#Now showing seasonality again at period 7, but not 14, 21.
#Will increase SMA by 1, rerun, and evaluate.
outT4=arima(xT,order=c(0,1,0),seasonal=list(order=c(0,2,2),period=7))
outT4
coeftest(outT4)
par(mfrow=c(1,2))
acf(outT4$residuals,lag.max=60)
pacf(outT4$residuals,lag.max=60)

#All are significant.
#aic = 6201.55
#Will now specify p and q orders
eacf(outT4$residuals)

#Based on EACF, we evaluate ARMA(2,4)
#Though MA8 and AR8 seem to reflect in ACF/PACF, we may have to disregard these possibilities.
#p and q values must not be over 7 because our period is 7.

### ARIMA(2,1,4) x SARIMA (0,2,2) ###
model1T1=arima(xT,order=c(2,1,4),seasonal=list(order=c(0,2,2),period=7))
model1T1
acf(model1T1$residuals,lag.max=60)
pacf(model1T1$residuals,lag.max=60)
Box.test(model1T1$residuals,lag=7,type='Ljung')
coeftest(model1T1)

#Most of the residuals are close to 0 = white noise process
#All are significant
#AIC = 6125.24
#Check Roots
polyroot(c(1,-model1T1$coef[1:2]))
abs(polyroot(c(1,-model1T1$coef[1:2])))
polyroot(c(1,model1T1$coef[3:6]))

#|AR Roots| are above 1.
#No matching roots.
#Will account for intercept/drift, rerun, and evaluate.
nT=length(xT)
model1T2=arima(xT,order=c(2,1,4),xreg=1:nT,seasonal=list(order=c(0,2,2),period=7))
model1T2
coeftest(model1T2)

#Intercept not significant.
#All are significant.
#Residuals are white.
#aic = 6198.14

#model1T1 is the only possible model in this case.

############### Reservations Revenue ####################
par(mfrow=c(1,3))
plot(park$Date,xR,type='l',main='Reservations Revenue',xlab='Date',ylab='$')
acf(xR,lag.max=60,main='Reservations ACF')
pacf(xR,lag.max=60,main='Reservations PACF')
adf.test(xR)

#Dickey Fuller test's p-value = 0.045; data is stationary
#We can tell that the data shows strong seasonality at period 7.
#Autocorrelation is strong at lags 7, 14, 21, 28, etc.
#PACF is significant around lag 7.
#We consider Seasonal AR1 model.
outR1=arima(xR,order=c(0,0,0),seasonal=list(order=c(1,0,0),period=7))
outR1
coeftest(outR1)
acf(outR1$residuals,lag.max=60)
pacf(outR1$residuals,lag.max=60)

#SAR1 and Intercept are significant.
#AIC = 5992.31
#We will check the ACF and PACF of the residuals.
#Still some seasonality at period 7; will increase SMA by 1.
outR2=arima(xR,order=c(0,0,0),seasonal=list(order=c(1,0,1),period=7))
outR2
coeftest(outR2)
acf(outR2$residuals,lag.max=60)
pacf(outR2$residuals,lag.max=60)

#All are significant.
#We will check the ACF and PACF of the residuals.
#AIC improves to 5947.86
#Looks like no more seasonality at 7.
#Lets evaluate AR and MA orders
eacf(outR2$residuals)

#Based on the EACF, ARMA(1,5) is a possibility
#Based on ACF and EACF, MA4 is a possibility.
#Based on the PACF, AR1 is a possibility.

### ARIMA(1,0,5) x SARIMA (1,0,1) ###
model1R1=arima(xR,order=c(1,0,5),seasonal=list(order=c(1,0,1),period=7))
model1R1

#NANs were produced.
#Reduce AR and MA by 1.
model1R2=arima(xR,order=c(0,0,4),seasonal=list(order=c(1,0,1),period=7))
model1R2
acf(model1R2$residuals,lag.max=60)
pacf(model1R2$residuals,lag.max=60)
Box.test(model1R2$residuals,lag=7,type='Ljung')
coeftest(model1R2)

#All are significant.
#AIC = 5783.96
#Conduct a roots test on SAR1 and SMA1.
polyroot(c(1,-model1R2$coef[5]))
polyroot(c(1,model1R2$coef[6]))

#SAR1 and SMA1 are not matching;
#however, SAR1 coeffient is over 1
#Will rerun model with SAR1 decrease and random walk increase by 1
nR=length(xR)
model1R3=arima(xR,order=c(0,0,4),xreg=1:nR,seasonal=list(order=c(0,1,1),period=7))
model1R3
acf(model1R3$residuals,lag.max=60)
pacf(model1R3$residuals,lag.max=60)
Box.test(model1R3$residuals,lag=7,type='Ljung')
coeftest(model1R3)

#All are significant.
#Residuals are white
#AIC improves to 5648.97

### ARIMA(1,0,0) x SARIMA (1,0,1) ###
model2R1=arima(xR,order=c(1,0,0),seasonal=list(order=c(1,0,1),period=7))
model2R1
acf(model2R1$residuals,lag.max=60)
pacf(model2R1$residuals,lag.max=60)
Box.test(model2R1$residuals,lag=7,type='Ljung')
coeftest(model2R1)

#Residuals are white.
#All are significant.
#AIC = 5779.64
#Lets test roots
polyroot(c(1,-model2R1$coef[1]))
abs(polyroot(c(1,-model2R1$coef[1])))
polyroot(c(1,model2R1$coef[2:3]))
abs(polyroot(c(1,-model2R1$coef[2:3])))

#|AR Root| is greater than 1.
#SAR1 coefficient is close 1, so there maybe a seasonal unit root.
#Will rerun model; decrease SAR1 and increase random walk by 1
model2R2=arima(xR,order=c(1,0,0),xreg=1:nR,seasonal=list(order=c(0,1,1),period=7))
model2R2
acf(model2R2$residuals,lag.max=60)
pacf(model2R2$residuals,lag.max=60)
Box.test(model2R2$residuals,lag=7,type='Ljung')
coeftest(model2R2)

#All are significant
#Residuals are white.
#AIC improves to 5644.98
#Evaluate the best performing and select the best 1.
model1R1$aic
model1R2$aic
model1R3$aic
model2R1$aic
model2R2$aic

#model2R2 is the best model for reservations.

############### HUB Revenue ####################

par(mfrow=c(1,3))
plot(park$Date,xH,type='l',main='Hub Revenue',xlab='Date',ylab='$')
acf(xH,lag.max=60,main='Hub Revenue ACF')
pacf(xH,lag.max=60,main='Hub Revenue PACF')
adf.test(xH)

#Dickey Fuller test's p-value = 0.01; data is stationary.
#We can tell that the data shows strong seasonality at period 7.
#Autocorrelation is strong at lags 7, 14, 21, 28, etc.
#PACF is significant around lag 7.
#We consider Seasonal D to the model.
outH1=arima(xH,order=c(0,0,0),seasonal=list(order=c(0,1,0),period=7))
outH1
acf(outH1$residuals,lag.max=60)
pacf(outH1$residuals,lag.max=60)

#AIC = 7220.6
#Per ACF, we still see some seasonality at 7, but not 14, 21, etc.
#Per PACF, we will see some at 7.
#We will increase SMA by 1, rerun, and evaluate.
outH2=arima(xH,order=c(0,0,0),seasonal=list(order=c(0,1,1),period=7))
outH2
coeftest(outH2)
acf(outH2$residuals,lag.max=60)
pacf(outH2$residuals,lag.max=60)

#SMA1 is significant.
#AIC improves to 7064.82
#PACF does not show anymore seasonality at 7.
#From ACF, we do some something at 7, but not at 14, 21, etc.
#We will increase SAR1, rerun, and evaluate.
outH3=arima(xH,order=c(0,0,0),seasonal=list(order=c(1,1,1),period=7))
outH3
coeftest(outH3)
acf(outH3$residuals,lag.max=60)
pacf(outH3$residuals,lag.max=60)

#All are significant
#aic = 7062.55
#No more seasonality at lags 7.
#We now evaluate AR and MA orders
eacf(outH3$residuals)

#Based on EACF and ACF, we will evaluate MA1.
#Based on PACF, AR2.

### ARIMA(2,0,0) x SARIMA (1,1,1) ###
model1H1=arima(xH,order=c(2,0,0),seasonal=list(order=c(1,1,1),period=7))
model1H1
acf(model1H1$residuals,lag.max=60)
pacf(model1H1$residuals,lag.max=60)
Box.test(model1H1$residuals,lag=7,type='Ljung')
coeftest(model1H1)

#All are significant.
#Residuals are white
#AIC = 6992.01
#Lets look over the roots.
polyroot(c(1,-model1H1$coef[1:2]))
abs(polyroot(c(1,-model1H1$coef[1:2])))
polyroot(c(1,model1H1$coef[3:4]))

#|AR Roots| are greater than 1
#Will add intercept/drift, rerun, and evaluate
nH=length(xH)
model1H2=arima(xH,order=c(2,0,0),xreg=1:nH,seasonal=list(order=c(1,1,1),period=7))
model1H2
acf(model1H2$residuals,lag.max=60)
pacf(model1H2$residuals,lag.max=60)
Box.test(model1H2$residuals,lag=7,type='Ljung')
coeftest(model1H2)

#SAR is not significant.
#Residuals are white.
#AIC = 6986.36
#Will remove SAR, rerun, and evaluate.
model1H3=arima(xH,order=c(2,0,0),xreg=1:nH,seasonal=list(order=c(0,1,1),period=7))
model1H3
acf(model1H3$residuals,lag.max=60)
pacf(model1H3$residuals,lag.max=60)
Box.test(model1H3$residuals,lag=7,type='Ljung')
coeftest(model1H3)

#All are significant.
#Residuals are white.
#aic = 6984.36

### ARIMA(0,0,1) x SARIMA (1,1,1) ###
model2H1=arima(xH,order=c(0,0,1),seasonal=list(order=c(1,1,1),period=7))
model2H1
acf(model2H1$residuals,lag.max=60)
pacf(model2H1$residuals,lag.max=60)
Box.test(model2H1$residuals,lag=7,type='Ljung')
coeftest(model2H1)

#SAR is not significant
#Residuals are white.
#aic = 6994.33
#Will add intercept/drift, rerun, and evaluate.
model2H2=arima(xH,order=c(0,0,1),xreg=1:nH,seasonal=list(order=c(1,1,1),period=7))
model2H2
acf(model2H2$residuals,lag.max=60)
pacf(model2H2$residuals,lag.max=60)
Box.test(model2H2$residuals,lag=7,type='Ljung')
coeftest(model2H2)

#SAR is still not significant.
#Residuals are white.
#aic = 6987.73
#Will remove SAR, rerun, and evaluate.
model2H3=arima(xH,order=c(0,0,1),xreg=1:nH,seasonal=list(order=c(0,1,1),period=7))
model2H3
acf(model2H3$residuals,lag.max=60)
pacf(model2H3$residuals,lag.max=60)
Box.test(model2H3$residuals,lag=7,type='Ljung')
coeftest(model2H3)

#Lets evaluate the better models and select the best.
model1H1$aic
model1H2$aic
model1H3$aic
model2H1$aic
model2H2$aic
model2H3$aic

#model1H3 is the better model for Hub (T1, T2W, LTL1, ECON, Lot7) Revenue.

############### VALET Revenue ####################
par(mfrow=c(1,3))
plot(park$Date,xV,type='l',main='VALET Revenue',xlab='Date',ylab='$')
acf(xV,lag.max=60,main='VALET Revenue ACF')
pacf(xV,lag.max=60,main='VALET Revenue PACF')
adf.test(xV)

#Dickey Fuller test's p-value = 0.99; data is not stationary.
#The variance isn't too bad, so we will not take the log.
#However, we will take the difference.
dxV=diff(xV)
plot(dxV,type='l',main='VALET Revenue DIFF',xlab='Date',ylab='$')
acf(dxV,lag.max=60,main='VALET Revenue ACF Diff')
pacf(dxV,lag.max=60,main='VALET Revenue PACF Diff')
adf.test(dxV)

#P-value is 0.01 = the data is now stationary.
#We can tell that the data shows strong seasonality at lag 7.
#ACF is very strong at lags 7, 14, 21, 28, etc.
#PACF is significant around lag 7.
#We consider Seasonal D first.
outV1=arima(dxV,order=c(0,0,0),seasonal=list(order=c(0,1,0),period=7))
outV1
acf(outV1$residuals,lag.max=60)
pacf(outV1$residuals,lag.max=60)

#AIC = 6092.63
#Per ACF, there is some seasonality at lag 7, but not 14, 21, etc.
#With PACF, there is some at Lag 7.
#Will add MA1, rerun, and evaluate.
outV2=arima(xV,order=c(0,1,0),seasonal=list(order=c(0,1,1),period=7))
outV2
coeftest(outV2)
acf(outV2$residuals,lag.max=60)
pacf(outV2$residuals,lag.max=60)

#SMA is significant.
#AIC = 5956.66
#ACF and PACF show no more seasonality lags at 7.
#We now evaluate AR and MA orders
eacf(outV2$residuals)

#Based on EACF and ACF, we will evaluate MA2
#Based on PACF, we will evaluate AR4

### ARIMA(4,1,0) x SARIMA (0,1,1) ###
model1V1=arima(xV,order=c(4,1,0),seasonal=list(order=c(0,1,1),period=7))
model1V1
acf(model1V1$residuals,lag.max=60)
pacf(model1V1$residuals,lag.max=60)
Box.test(model1V1$residuals,lag=7,type='Ljung')
coeftest(model1V1)

#All are significant.
#Residuals are white
#AIC = 5890.48
#Lets look over the roots.
polyroot(c(1,-model1V1$coef[1:4]))
abs(polyroot(c(1,-model1V1$coef[1:4])))
polyroot(c(1,model1V1$coef[5]))

#|AR roots| > 1
#Will add intercept/drift, rerun, and evaluate.
nV=length(xV)
model1V2=arima(xV,order=c(4,1,0),xreg=1:nV,seasonal=list(order=c(0,1,1),period=7))

#Error message with intercept/drift. Will omit going forward.

### ARIMA(0,1,2) x SARIMA (0,1,1) ###
model2V1=arima(xV,order=c(0,1,2),seasonal=list(order=c(0,1,1),period=7))
model2V1
acf(model2V1$residuals,lag.max=60)
pacf(model2V1$residuals,lag.max=60)
Box.test(model2V1$residuals,lag=7,type='Ljung')
coeftest(model2V1)

#All are significant.
#Residuals are white.
#aic = 5881.86
#Because this is a  MA model, no need to evaluate the roots.

#Lets evaluate all better models and select the best.
model1V1$aic
model2V1$aic

#model2V1 is the better model.

##load new data
newdata=read.csv("2019.04.14-2019.04.30_Parking_NewData.csv",header=T)
head(newdata)
str(newdata)

#Convert Date from Factor to Date
newdata$Date=as.Date(newdata$Date,format="%m/%d/%Y")
str(newdata)

###Look at Revenue Variables 1 at a time
#set each revenue variable
newxT=(newdata$T2PP_Revenue)
newxR=(newdata$Reservations_Revenue)
newxH=(newdata$Hub_Revenue)
newxV=(newdata$Valet_Rev)

###ROLLING FORECAST###

#model2V1 is the best model to predict Valet Revenue.

print(c(model1R1$aic, model1R2$aic, model1R3$aic, model2R1$aic, model2R2$aic))
print(c(model1H1$aic, model1H2$aic, model2H1$aic, model2H2$aic, model2H3$aic))
print(c(model1V1$aic, model2V1$aic))

source("rolling.forecast.R")

## T2PP
#No rolling forecast.
#Only 1 model.

## Reservations
print(c(model1R2$aic, model1R3$aic, model2R1$aic, model2R2$aic))

print(rolling.forecast(xR,6,250,c(0,0,4),seasonal=list(order=c(1,0,1),period=7)))
print(rolling.forecast.drift(xR,6,250,c(0,0,4),seasonal=list(order=c(0,1,1),period=7)))
print(rolling.forecast(xR,6,250,c(1,0,0),seasonal=list(order=c(1,0,1),period=7)))
print(rolling.forecast.drift(xR,6,250,c(1,0,0),seasonal=list(order=c(0,1,1),period=7)))
errorR1=rolling.forecast(xR,6,250,c(0,0,4),seasonal=list(order=c(1,0,1),period=7))
errorR2=rolling.forecast.drift(xR,6,250,c(0,0,4),seasonal=list(order=c(0,1,1),period=7))
errorR3=rolling.forecast(xR,6,250,c(1,0,0),seasonal=list(order=c(1,0,1),period=7))
errorR4=rolling.forecast.drift(xR,6,250,c(1,0,0),seasonal=list(order=c(0,1,1),period=7))

errorR=c(errorR1,errorR2,errorR3,errorR4)

#Reservation Error Plot
par(mfrow=c(1,1))
plot(errorR1,type='l',ylim=c(min(errorR),max(errorR)),main='Reservations - Rolling Forecasting Errors',xlab='Forecast horizon',ylab='Error')
lines(errorR2,col=2)
lines(errorR3,col=3)
lines(errorR4,col=4)
legend.text=c("Model 1 - model1R2","Model 2 - model1R3","Model 3 - model2R1","Model 4 - model2R2")
legend("bottomright",legend.text,lty=rep(1,4),col=1:4)

#Based on rolling forecast, model 4 (model2R2)
#Based on AIC, model2R2 is also better.

##HUB
print(c(model1H1$aic, model1H2$aic, model1H3$aic, model2H1$aic, model2H2$aic, model2H3$aic))

print(rolling.forecast(xH,6,250,c(2,0,0),seasonal=list(order=c(1,1,1),period=7)))
print(rolling.forecast.drift(xH,6,250,c(2,0,0),seasonal=list(order=c(1,1,1),period=7)))
print(rolling.forecast.drift(xH,6,250,c(2,0,0),seasonal=list(order=c(0,1,1),period=7)))
print(rolling.forecast(xH,6,250,c(0,0,1),seasonal=list(order=c(1,1,1),period=7)))
print(rolling.forecast.drift(xH,6,250,c(0,0,1),seasonal=list(order=c(1,1,1),period=7)))
print(rolling.forecast.drift(xH,6,250,c(0,0,1),seasonal=list(order=c(0,1,1),period=7)))
errorH1=rolling.forecast(xH,6,250,c(2,0,0),seasonal=list(order=c(1,1,1),period=7))
errorH2=rolling.forecast.drift(xH,6,250,c(2,0,0),seasonal=list(order=c(1,1,1),period=7))
errorH3=rolling.forecast.drift(xH,6,250,c(2,0,0),seasonal=list(order=c(0,1,1),period=7))
errorH4=rolling.forecast(xH,6,250,c(0,0,1),seasonal=list(order=c(1,1,1),period=7))
errorH5=rolling.forecast.drift(xH,6,250,c(0,0,1),seasonal=list(order=c(1,1,1),period=7))
errorH6=rolling.forecast.drift(xH,6,250,c(0,0,1),seasonal=list(order=c(0,1,1),period=7))

errorH=c(errorH1,errorH2,errorH3,errorH4,errorH5,errorH6)

#Hub Error Plot
par(mfrow=c(1,1))
plot(errorH1,type='l',ylim=c(min(errorH),max(errorH)),main='Hub - Rolling Forecasting Errors',xlab='Forecast horizon',ylab='Error')
lines(errorH2,col=2)
lines(errorH3,col=3)
lines(errorH4,col=4)
lines(errorH5,col=5)
lines(errorH6,col=6)
legend.text=c("Model 1 - model1H1","Model 2 - model1H2","Model 3 - model1H3","Model 4 - model2H1","Model 5 - model2H2","Model 6 - model2H3")
legend("bottomright",legend.text,lty=rep(1,7),col=1:7)

##Per Rolling Forecast, Model 3 (model2H1) is the better model.
##Per AIC, Model 5 (model1H3) is the better model.

##VALET

print(c(model1V1$aic, model2V1$aic))

print(rolling.forecast(xV,6,250,c(4,0,1),seasonal=list(order=c(0,1,1),period=7)))
print(rolling.forecast(xV,6,250,c(0,1,2),seasonal=list(order=c(0,1,1),period=7)))
errorV1=rolling.forecast(xV,6,250,c(4,0,1),seasonal=list(order=c(0,1,1),period=7))
errorV2=rolling.forecast(xV,6,250,c(0,1,2),seasonal=list(order=c(0,1,1),period=7))

errorV=c(errorV1,errorV2)

par(mfrow=c(1,1))
plot(errorV1,type='l',ylim=c(min(errorV),max(errorV)),main='Valet Models - Rolling Forecasting Errors',xlab='Forecast horizon',ylab='Error')
lines(errorV2,col=2)
legend.text.V=c("Model 1 - model1V1","Model 2 - model2V1")
legend("bottomright",legend.text.V,lty=rep(1,2),col=1:2)

##Per Rolling Forecast, Model 2 (model2V1) is better.
##Per AIC, Model 2 (model2V1) is also better.

########## GRAPHS

### T2PP
###model1T1 (only 1 model)

ppTT=predict(model1T1,17)
ppTT$pred
ppTT$se

#Plot = data + predictions

par(mfrow=c(1,1))
nnT=length(xT)	#length of data
ntT=17	#forecast horizon
nbT=316	#number of data points to plot
ttT=(nnT-nbT):nnT	#indexes of data points to plot
xxxT=xT[ttT]		#data you want to plot
rrT=range(c(xxxT,ppTT$pred+2*ppTT$se,ppTT$pred-2*ppTT$se))	#minimum and maximum y values in the plot
plot(ttT,xxxT,pch=3,xlim=c(nnT-nbT,nnT+ntT),ylim=rrT,main='17-Day Prediction for T2PP')	
lines(ttT,xxxT)	#observed values
points(nnT+1:ntT,ppTT$pred,pch=2,col='red',type='o')	#predicted values
lines(nnT+1:ntT,ppTT$pred+2*ppTT$se,lty=2,col='red')	#upper bound of predicted interval
lines(nnT+1:ntT,ppTT$pred-2*ppTT$se,lty=2,col='red')	#lower bound of predicted interval

#Plot = close up of predictions
nnT=length(xT)	#length of data
ntT=17	#forecast horizon
nbT=10	#number of data points to plot
ttT=(nnT-nbT):nnT	#indexes of data points to plot
xxxT=xT[ttT]		#data you want to plot
rrT=range(c(xxxT,ppTT$pred+2*ppTT$se,ppTT$pred-2*ppTT$se))	#minimum and maximum y values in the plot
plot(ttT,xxxT,pch=3,xlim=c(nnT-nbT,nnT+ntT),ylim=rrT,main='17-Day Prediction for T2PP')	
lines(ttT,xxxT)	#observed values
points(nnT+1:ntT,ppTT$pred,pch=2,col='red',type='o')	#predicted values
lines(nnT+1:ntT,ppTT$pred+2*ppTT$se,lty=2,col='red')	#upper bound of predicted interval
lines(nnT+1:ntT,ppTT$pred-2*ppTT$se,lty=2,col='red')	#lower bound of predicted interval

head(newxT)
head(ppTT$pred)

### Reservations
### model2R2

hR=17
nR=length(xR)
ppRR=predict(model2R2,17,newxreg=(nR+1):(nR+hR))
ppRR$pred
ppRR$se

head(newxR)

#Plot = data + predictions
par(mfrow=c(1,1))
nnRR=length(xR)	#length of data
ntRR=17	#forecast horizon
nbRR=316	#number of data points to plot
ttRR=(nnRR-nbRR):nnRR	#indexes of data points to plot
xxxRR=xR[ttRR]		#data you want to plot
rrRR=range(c(xxxRR,ppRR$pred+2*ppRR$se,ppRR$pred-2*ppRR$se))	#minimum and maximum y values in the plot
plot(ttRR,xxxRR,pch=3,xlim=c(nnRR-nbRR,nnRR+ntRR),ylim=rrRR,main='17-Day Prediction for Reservations')	
lines(ttRR,xxxRR)	#observed values
points(nnRR+1:ntRR,ppRR$pred,pch=2,col='red',type='o')	#predicted values
lines(nnRR+1:ntRR,ppRR$pred+2*ppRR$se,lty=2,col='red')	#upper bound of predicted interval
lines(nnRR+1:ntRR,ppRR$pred-2*ppRR$se,lty=2,col='red')	#lower bound of predicted interval

#Plot = close up of predictions
nnRR=length(xR)	#length of data
ntRR=17	#forecast horizon
nbRR=10	#number of data points to plot
ttRR=(nnRR-nbRR):nnRR	#indexes of data points to plot
xxxRR=xR[ttRR]		#data you want to plot
rrRR=range(c(xxxRR,ppRR$pred+2*ppRR$se,ppRR$pred-2*ppRR$se))	#minimum and maximum y values in the plot
plot(ttRR,xxxRR,pch=3,xlim=c(nnRR-nbRR,nnRR+ntRR),ylim=rrRR,main='17-Day Prediction for Reservations')	
lines(ttRR,xxxRR)	#observed values
points(nnRR+1:ntRR,ppRR$pred,pch=2,col='red',type='o')	#predicted values
lines(nnRR+1:ntRR,ppRR$pred+2*ppRR$se,lty=2,col='red')	#upper bound of predicted interval
lines(nnRR+1:ntRR,ppRR$pred-2*ppRR$se,lty=2,col='red')	#lower bound of predicted interval

### HUB
model2H1
model1H3

hH=17
nH=length(xH)
ppH=predict(model2H1,17)
ppH$pred
ppH$se

ppH=predict(model1H3,17,newxreg=(nH+1):(nH+hH))
ppH$pred
ppH$se

#SE for model1H3 is the lowest.  We will use this model.

head(newxH)


#Plot = data + predictions
par(mfrow=c(1,1))
nnH=length(xH)	#length of data
ntH=17	#forecast horizon
nbH=316	#number of data points to plot
ttH=(nnH-nbH):nnH	#indexes of data points to plot
xxxH=xH[ttH]		#data you want to plot
rrH=range(c(xxxH,ppH$pred+2*ppH$se,ppH$pred-2*ppH$se))	#minimum and maximum y values in the plot
plot(ttH,xxxH,pch=3,xlim=c(nnH-nbH,nnH+ntH),ylim=rrH,main='17-Day Prediction for Hub')	
lines(ttH,xxxH)	#observed values
points(nnH+1:ntH,ppH$pred,pch=2,col='red',type='o')	#predicted values
lines(nnH+1:ntH,ppH$pred+2*ppH$se,lty=2,col='red')	#upper bound of predicted interval
lines(nnH+1:ntH,ppH$pred-2*ppH$se,lty=2,col='red')	#lower bound of predicted interval

#Plot = close up of predictions
nnH=length(xH)	#length of data
ntH=17	#forecast horizon
nbH=10	#number of data points to plot
ttH=(nnH-nbH):nnH	#indexes of data points to plot
xxxH=xH[ttH]		#data you want to plot
rrH=range(c(xxxH,ppH$pred+2*ppH$se,ppH$pred-2*ppH$se))	#minimum and maximum y values in the plot
plot(ttH,xxxH,pch=3,xlim=c(nnH-nbH,nnH+ntH),ylim=rrH,main='17-Day Prediction for Hub')	
lines(ttH,xxxH)	#observed values
points(nnH+1:ntH,ppH$pred,pch=2,col='red',type='o')	#predicted values
lines(nnH+1:ntH,ppH$pred+2*ppH$se,lty=2,col='red')	#upper bound of predicted interval
lines(nnH+1:ntH,ppH$pred-2*ppH$se,lty=2,col='red')	#lower bound of predicted interval

### VALET
### model2V1

ppV=predict(model2V1,17)
ppV$pred
ppV$se

head(newxV)

#Plot = data + predictions
par(mfrow=c(1,1))
nnV=length(xV)	#length of data
ntV=17	#forecast horizon
nbV=316	#number of data points to plot
ttV=(nnV-nbV):nnV	#indexes of data points to plot
xxxV=xV[ttV]		#data you want to plot
rrV=range(c(xxxV,ppV$pred+2*ppV$se,ppV$pred-2*ppV$se))	#minimum and maximum y values in the plot
plot(ttV,xxxV,pch=3,xlim=c(nnV-nbV,nnV+ntV),ylim=rrV,main='17-Day Prediction for Valet')	
lines(ttV,xxxV)	#observed values
points(nnV+1:ntV,ppV$pred,pch=2,col='red',type='o')	#predicted values
lines(nnV+1:ntV,ppV$pred+2*ppV$se,lty=2,col='red')	#upper bound of predicted interval
lines(nnV+1:ntV,ppV$pred-2*ppV$se,lty=2,col='red')	#lower bound of predicted interval

#Plot = close up of predictions
nnV=length(xV)	#length of data
ntV=17	#forecast horizon
nbV=10	#number of data points to plot
ttV=(nnV-nbV):nnV	#indexes of data points to plot
xxxV=xV[ttV]		#data you want to plot
rrV=range(c(xxxV,ppV$pred+2*ppV$se,ppV$pred-2*ppV$se))	#minimum and maximum y values in the plot
plot(ttV,xxxV,pch=3,xlim=c(nnV-nbV,nnV+ntV),ylim=rrV,main='17-Day Prediction for Valet')	
lines(ttV,xxxV)	#observed values
points(nnV+1:ntV,ppV$pred,pch=2,col='red',type='o')	#predicted values
lines(nnV+1:ntV,ppV$pred+2*ppV$se,lty=2,col='red')	#upper bound of predicted interval
lines(nnV+1:ntV,ppV$pred-2*ppV$se,lty=2,col='red')	#lower bound of predicted interval

########################################################################
################### PARKING TRANSACTION ANALYSIS #######################
#######################################################################

#project r

data=read.csv("project.csv",header=T)

plot(x,type='l',col='skyblue2',main='Parking transactions',ylim=c(0,7500),ylab='$')
points(r,type='l',col='seagreen2')
points(H,type='l',col='cyan4')
points(v,type='l',col='midnightblue')
legend('topright',legend=c("T2 Parking Plaza", "Reservations", "Hub (5 Lots)","Valet"),col=c('skyblue2','seagreen2','cyan4','midnightblue'),lty=c(1,1,1))



mean(x)
mean(r)
mean(H)
mean(v)


median(x)
median(r)
median(H)
median(v)



range(x)
range(r)
range(H)
range(v)


var(x)
var(r)
var(H)
var(v)


sd(x)
sd(r)
sd(H)
sd(v)

####T2


plot(data[,3])
plot(data[,3],type='l')
x=data[,3]
plot(x)
plot(x,type='l')
plot(data[,3],type='l')
plot(x,type='l')
acf(x)
pacf(x)
eacf(x)
plot(x,type='l')
logx=log(x)
plot(logx,type='l')
diffx=diff(logx)
plot(diffx,type='l')
acf(diffx)
pacf(diffx)
eacf(diffx)
Box.test(diffx)
out1=arima(diffx,order=c(0,0,0),seasonal=list(order=c(1,0,0),period=7))
out1
coeftest(out1)
acf(out1$residuals,lag.max=25)
pacf(out1$residuals,lag.max=25)
out2=arima(logx,order=c(0,1,0),seasonal=list(order=c(1,0,1),period=7))
out2
coeftest(out2)
acf(out2$residuals,lag.max=25)
pacf(out2$residuals,lag.max=25)
eacf(out2$residuals)
### arma 2,5
out3=arima(logx,order=c(2,1,5),seasonal=list(order=c(1,0,1),period=7))
out3
acf(out3$residuals)
pacf(out3$residuals)
Box.test(out3$residuals,lag=12,type='Ljung')
coeftest(out3)

#AR4 for T2 transactions

out4=arima(logx,order=c(4,1,0),seasonal=list(order=c(1,0,1),period=7))
out4
coeftest(out4)
out41=arima(logx,order=c(3,1,0),seasonal=list(order=c(1,0,1),period=7))
out41
coeftest(out41)
acf(out41$residuals)
pacf(out41$residuals)
Box.test(out41$residuals,type='Ljung')

#MA4 for T2 transactions

out5=arima(logx,order=c(0,1,4),seasonal=list(order=c(1,0,1),period=7))
out5

coeftest(out5)
out51=arima(logx,order=c(0,1,4),seasonal=list(order=c(0,1,1),period=7))
out51
coeftest(out51)
acf(out51$residuals)
pacf(out51$residuals)
Box.test(out51$residuals)

### reservations

r=data[,5]
plot(r,type='l')
logr=log(r)
plot(logr,type='l')
diffr=diff(logr)

plot(diffr,type='l')
par(mfrow=c(1,2))
acf(diffr)
pacf(diffr)




outr1=arima(diffr,order=c(0,0,0),seasonal=list(order=c(1,0,0),period=7))
outr1
coeftest(outr1)
## intercept not significant 
acf(outr1$residuals)
pacf(outr1$residuals)

outr2=arima(logr,order=c(0,1,0),seasonal=list(order=c(1,0,1),period=7))
outr2
coeftest(outr2)
acf(outr2$residuals)
pacf(outr2$residuals)
eacf(outr2$residuals)

outr3=arima(logr,order=c(0,1,0),seasonal=list(order=c(1,1,1),period=7))
outr3
coeftest(outr3)
eacf(outr3$residuals)
acf(outr3$residuals)
pacf(outr3$residuals)

outr4=arima(diffr,order=c(0,0,0),seasonal=list(order=c(0,1,0),period=7))
outr4
coeftest(outr4)
acf(outr4$residuals)
pacf(outr4$residuals)

outr5=arima(logr,order=c(0,1,0),seasonal=list(order=c(0,1,1),period=7))
outr5
coeftest(outr5)
acf(outr5$residuals)
pacf(outr5$residuals)
eacf(outr5$residuals)

outr6=arima(logr,order=c(0,1,0),seasonal=list(order=c(0,1,1),period=7))
outr6
coeftest(outr6)
eacf(outr6$residuals)


arma 2,5


arma25=arima(logr,order=c(2,1,5),seasonal=list(order=c(1,1,1),period=7))
arma25
coeftest(arma25)
Box.test(arma25$residuals)
acf(arma25$residuals)

arma26=arima(logr,order=c(2,1,5),seasonal=list(order=c(0,1,1),period=7))
arma26
coeftest(arma26)

#### ARMA 2,4

arma24=arima(logr,order=c(2,1,4),seasonal=list(order=c(1,1,1),period=7))
arma24
acf(arma24$residuals)
pacf(arma24$residuals)
Box.test(arma24$residuals,lag=7,type='Ljung')

coeftest(arma24)

we remove the not significant coeficcients 

arma241=arima(logr,order=c(2,1,4),seasonal=list(order=c(1,1,1),period=7),fixed=c(0,NA,0,NA,0,0,NA,NA))
arma241
coeftest(arma241)

polyroot(c(1,-arma24$coef[1:2]))
abs(polyroot(c(1,-arma24$coef[1:2])))
polyroot(c(1,arma24$coef[3:8]))

acf(arma241$residuals)
pacf(arma241$residuals)
Box.test(arma241$residuals,lag=7,type='Ljung')

residuals not clean and it does not pass the box test, aic is worse so we stay with model 2,4

arma242=arima(logr,order=c(2,1,4),seasonal=list(order=c(0,1,1),period=7))
arma242
coeftest(arma242)
Box.test(arma242$residuals,lag=7,type='Ljung')
acf(arma242$residuals)

arma243=arima(logr,order=c(2,1,4),seasonal=list(order=c(0,1,1),period=7),fixed=c(0,NA,0,NA,0,0,NA))
Box.test(arma243$residuals,type='Ljung')

arma244=arima(logr,order=c(2,1,2),seasonal=list(order=c(0,1,1),period=7))
arma244
coeftest(arma244)
Box.test(arma244$residuals,type='Ljung')

arma245=arima(logr,order=c(1,1,1),seasonal=list(order=c(0,1,1),period=7))
arma245
coeftest(arma245)

Box.test(arma245$residuals,type='Ljung')

acf(arma245$residuals)


AR6



Rar6=arima(logr,order=c(6,1,0),seasonal=list(order=c(1,1,1),period=7))
Rar6

coeftest(Rar6)

acf(Rar6$residuals)
pacf(Rar6$residuals)

Box.test(Rar6$residuals,lag=7,type='Ljung')

polyroot(c(1,-Rar6$coef[1:6]))
abs(polyroot(c(1,-Rar6$coef[1:6])))

polyroot(c(1,Rar6$coef[7:8]))



MA8

Rma8=arima(logr,order=c(0,1,8),seasonal=list(order=c(1,1,1),period=7))

Rma8

coeftest(Rma8)
Box.test(Rma8$residuals,lag=7,type='Ljung')
acf(Rma8$residuals,lag.max=60)
pacf(Rma8$residuals,lag.max=60)


polyroot(c(1,Rma8$coef[1:10]))

abs(polyroot(c(1,Rma8$coef[1:10])))

Rma81=arima(logr,order=c(0,1,8),seasonal=list(order=c(1,1,1),period=7),fixed=c(NA,NA,NA,0,0,0,NA,NA,NA,NA))
Rma81

coeftest(Rma81)
Box.test(Rma81$residuals,lag=7,type='Ljung')
acf(Rma81$residuals)
pacf(Rma81$residuals)

polyroot(c(1,Rma81$coef[1:5]))


abs(polyroot(c(1,Rma81$coef[1:8])))

polyroot(c(1,Rma81$coef[9:10]))

abs(polyroot(c(1,Rma81$coef[9:10])))

Rma82=arima(logr,order=c(0,1,8),seasonal=list(order=c(1,1,1),period=7),fixed=c(NA,NA,NA,0,0,0,NA,NA,NA,0))


#### HUB

plot(data[,7],type='l')
H=data[,7]
acf(H)
pacf(H)
adf.test(H)

data not stationary so we take the log
logH=log(H)
plot(logH,type='l')

still not stationary so we take differance 

dH=diff(log(H))
plot(dH,type='l')

acf(dH)
pacf(dH)

H1=arima(dH,order=c(0,0,0),seasonal=list(order=c(1,0,0),period=7))
H1

acf(H1$residuals)
pacf(H1$residuals)

H2=arima(dH,order=c(0,0,0),seasonal=list(order=c(1,0,1),period=7))
H2


coeftest(H1)

acf(H2$residuals)
pacf(H2$residuals)
eacf(H2$residuals)

we will try AR6, MA6 and ARMA 3,4

AR6

Har6=arima(logH,order=c(6,1,0),seasonal=list(order=c(1,0,1),period=7))
Har6
coeftest(Har6)

polyroot(c(1,-Har6$coef[1:6]))
abs(polyroot(c(1,-Har6$coef[1:6])))
polyroot(c(1,-Har6$coef[7:8]))
abs(polyroot(c(1,-Har6$coef[7:8])))

Har60=arima(logH,order=c(6,1,0),seasonal=list(order=c(0,1,1),period=7))
Har60


acf(Har6$residuals)
pacf(Har6$residuals)
Box.test(Har6$residuals,lag=7,type='Ljung')
model looks good, SAR1 close to 1 we will take out and re evaluate

Har61=arima(logH,order=c(6,1,0),seasonal=list(order=c(1,0,1),period=7),fixed=c(NA,NA,NA,NA,NA,NA,0,NA))

aic is worse so we stay with ar6

Har62=arima(logH,order=c(6,1,0),seasonal=list(order=c(0,1,1),period=7))

##### MA6

Hma6=arima(logH,order=c(0,1,6),seasonal=list(order=c(1,0,1),period=7))
Hma6

coeftest(Hma6)
acf(Hma6$residuals)
pacf(Hma6$residuals)

polyroot(c(1,Hma6$coef[1:8]))
abs(polyroot(c(1,Hma6$coef[1:8])))

we remove non significant coefficients
Hma61=arima(logH,order=c(0,1,6),seasonal=list(order=c(1,0,1),period=7),fixed=c(NA,NA,NA,NA,0,0,NA,NA))

acf(Hma61$residuals)
pacf(Hma61$residuals)
Box.test(Hma61$residuals,type='Ljung')

aic decreases a little bit Hma6 better model

arma 3,4

Harma34=arima(logH,order=c(3,1,4),seasonal=list(order=c(1,0,1),period=7))
Harma34

coeftest(Harma34)
acf(Harma34$residuals)
pacf(Harma34$residuals)
Box.test(Harma34$residuals)
polyroot(c(1,-Harma34$coef[1:3]))
abs(polyroot(c(1,-Harma34$coef[1:3])))
polyroot(c(1,Harma34$coef[4:9]))
abs(polyroot(c(1,Harma34$coef[4:9])))

no matching roots

Harma341=arima(logH,order=c(3,1,4),seasonal=list(order=c(1,0,1),period=7),fixed=c(NA,0,NA,NA,0,NA,NA,NA,NA))

Harma341
coeftest(Harma341)

we stay with Harma34


Har6$aic
Hma6$aic
Harma34$aic

arma 3,4 better model 


####### valet


v=data[,9]
plot(v,type='l')

acf(v)
pacf(v)
adf.test(v)

data not stationary 

logv=log(v)
plot(logv,type='l')

adf.test(logv)

dv=diff(log(v))
plot(dv,type='l')
adf.test(dv)

acf(dv)
pacf(dv)


V1=arima(dv,order=c(0,0,0),seasonal=list(order=c(1,0,0),period=7))
V1
coeftest(V1)
acf(V1$residuals,lag.max=60)
pacf(V1$residuals)

still seasonality
we add SMA1

V2=arima(logv,order=c(0,1,0),seasonal=list(order=c(1,0,1),period=7))
V2
coeftest(V2)

acf(V2$residuals)
pacf(V2$residuals)
eacf(V2$residuals)


we will try MA2 and AR6

####MA2

vma2=arima(logv,order=c(0,1,2),seasonal=list(order=c(1,0,1),period=7))
vma2

acf(vma2$residuals)
pacf(vma2$residuals)
Box.test(vma2$residuals,type='Ljung')
coeftest(vma2)

sar1 is really close to 1 so we take it out 

vma3=arima(logv,order=c(0,1,2),seasonal=list(order=c(0,1,1),period=7))
vma3

acf(vma3$residuals)
pacf(vma3$residuals)
Box.test(vma3$residuals,type='Ljung')
coeftest(vma3)

better AIC


##### AR6


var6=arima(logv,order=c(6,1,0),seasonal=list(order=c(1,0,1),period=7))
var6
coeftest(var6)

we take Sar1 out because it is really close to 1


var61=arima(logv,order=c(6,1,0),seasonal=list(order=c(0,1,1),period=7))
var61

better AIC

coeftest(var61)
polyroot(c(1,-var61$coef[1:6]))
abs(polyroot(c(1,-var61$coef[1:6])))
polyroot(c(1,var61$coef[7]))

no matching roots, all ar roots over 1


acf(var61$residuals)
pacf(var61$residuals)
Box.test(var61$residuals)


AIC for models

vma2$aic
vma3$aic
var6$aic
var61$aic

vma3 best model


#####   rolling forecast 



## t2

source("rolling.forecast.R")
print(rolling.forecast(logx,7,250,c(2,1,5),seasonal=list(order=c(1,0,1))))
print(rolling.forecast(logx,7,250,c(3,1,0),seasonal=list(order=c(1,0,1))))
print(rolling.forecast(logx,7,250,c(0,1,4),seasonal=list(order=c(0,1,1))))

errort1=rolling.forecast(logx,7,250,c(2,1,5),seasonal=list(order=c(1,0,1)))
errort2=rolling.forecast(logx,7,250,c(3,1,0),seasonal=list(order=c(1,0,1)))
errort3=rolling.forecast(logx,7,250,c(0,1,4),seasonal=list(order=c(0,1,1)))

errort=c(errort1,errort2,errort3)

plot(errort1,type='l',ylim=c(min(errort),1))
lines(errort2,col=2)
lines(errort3,col=3)

legend.text=c("ARMA2,5","AR4","MA4")
legend("topright",legend.text,lty=rep(1,3),col=1:3)


model out41


##### reservations


error1=rolling.forecast(logr,7,250,c(2,1,4),seasonal=list(order=c(1,1,1)))
error2=rolling.forecast(logr,7,250,c(6,1,0),seasonal=list(order=c(1,1,1)))
error3=rolling.forecast(logr,7,250,c(0,1,8),seasonal=list(order=c(1,1,1),fixed=c(NA,NA,NA,0,0,0,NA,NA,NA,NA)))
error4=rolling.forecast(logr,7,250,c(2,1,2),seasonal=list(order=c(0,1,1)))
error5=rolling.forecast(logr,7,250,c(2,1,5),seasonal=list(order=c(1,1,1)))
error6=rolling.forecast(logr,7,250,c(2,1,5),seasonal=list(order=c(0,1,1)))
error=c(error1,error2,error3,error4,error5,error6)

plot(error1,type='l',ylim=c(min(error),3))
lines(error2,col=2)
lines(error3,col=3)
lines(error4,col=4)
lines(error5,col=5)


plot(error4,type='l')
plot(error5,type='l')
plot(error6,type='l')

print(rolling.forecast(logr,7,250,c(2,1,4),seasonal=list(order=c(1,1,1))))

Rar6 model



####### HUB


errorh1=rolling.forecast(logH,7,250,c(6,1,0),seasonal=list(order=c(1,0,1)))
errorh2=rolling.forecast(logH,7,250,c(0,1,6),seasonal=list(order=c(1,0,1)))
errorh3=rolling.forecast(logH,7,250,c(3,1,4),seasonal=list(order=c(0,1,1)))

errorh=c(errorh1,errorh2,errorh3)

plot(errorh1,type='l',ylim=c(min(errorh),2))
lines(errorh2,col=2)
lines(errorh3,col=3)

Har6 model


######### valet




errorv1=rolling.forecast(logv,7,250,c(0,1,2),seasonal=list(order=c(0,1,1)))
errorv2=rolling.forecast(logv,7,250,c(6,1,0),seasonal=list(order=c(0,1,1)))


errorv=c(errorv1,errorv2)

plot(errorv1,type='l',ylim=c(min(errorv),15))
lines(errorv2,col=2)

var61 model


##### predictions 

T2

out41
pp=predict(out41,17)
nn=length(logx)
nt=17
nb=17
tt=(nn-nb):nn
xxx=logx[tt]
rr=range(c(exp(xxx),exp(pp$pred+2*pp$se),exp(pp$pred-2*pp$se)))
plot(tt,exp(xxx),pch=3,xlim=c(nn-nb,nn+nt),ylim=rr)
lines(tt,exp(xxx))
points(nn+1:nt,exp(pp$pred),pch=2,col='red',type='o')
lines(nn+1:nt,exp(pp$pred+2*pp$se),lty=2,col='red')
lines(nn+1:nt,exp(pp$pred-2*pp$se),lty=2,col='red')

pp$pred
pp$se

##### Reservations

Rar6
ppr=predict(Rar6,17)
nnr=length(logr)
ntr=17
nbr=17
ttr=(nnr-nbr):nnr
xxxr=logr[ttr]
rrr=range(c(exp(xxxr),exp(ppr$pred+2*ppr$se),exp(ppr$pred-2*pp$se)))
plot(ttr,exp(xxxr),pch=3,xlim=c(nnr-nbr,nnr+ntr),ylim=rrr)
lines(ttr,exp(xxxr))
points(nnr+1:ntr,exp(ppr$pred),pch=2,col='red',type='o')
lines(nnr+1:ntr,exp(ppr$pred+2*ppr$se),lty=2,col='red')
lines(nnr+1:ntr,exp(ppr$pred-2*ppr$se),lty=2,col='red')



ppr$pred
ppr$se

ppr1=predict(x,17)
ppr1$pred



###### HUB 

Har6
pph=predict(Har6,17)
nnh=length(logH)
nth=17
nbh=17
tth=(nnh-nbh):nnh
xxxh=logH[tth]
rrh=range(c(exp(xxxh),exp(pph$pred+2*pph$se),exp(pph$pred-2*pph$se)))
plot(tth,exp(xxxh),pch=3,xlim=c(nnh-nbh,nnh+nth),ylim=rrh)
lines(tth,exp(xxxh))
points(nnh+1:nth,exp(pph$pred),pch=2,col='red',type='o')
lines(nnh+1:ntr,exp(pph$pred+2*pph$se),lty=2,col='red')
lines(nnh+1:nth,exp(pph$pred-2*ppr$se),lty=2,col='red')

pph$pred
pph$se


######## valet

var61
ppv=predict(var61,17)
nnv=length(logv)
ntv=17
nbv=17
ttv=(nnv-nbv):nnv
xxxv=logv[ttv]
rrv=range(c(exp(xxxv),exp(ppv$pred+2*ppv$se),exp(ppv$pred-2*ppv$se)))
plot(ttv,exp(xxxv),pch=3,xlim=c(nnv-nbv,nnh+ntv),ylim=rrv)
lines(ttv,exp(xxxv))
points(nnv+1:ntv,exp(ppv$pred),pch=2,col='red',type='o')
lines(nnv+1:ntv,exp(ppv$pred+2*ppv$se),lty=2,col='red')
lines(nnv+1:ntv,exp(ppv$pred-2*ppv$se),lty=2,col='red')