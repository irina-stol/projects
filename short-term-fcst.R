
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(MASS)
library(caret)
library(mgcv)


set.seed(1433)
# Read in St Austell Data
data_18<- read.csv('SAB_14-18new_Xymatrix_OA.csv')
data_14<- read.csv('SAB_08-15_Xymatrix_OA.csv')
data_18$sst <- rowMeans(data_18[,c('ECMWFsst','Eosst')])
data_14$sst <- rowMeans(data_14[,c('ECMWFsst','EOsst')])
data_18 <- data_18[,!names(data_18) %in% c('ECMWFsst','Eosst','ECMWFssrd')]
data_14 <- data_14[,!names(data_14) %in% c('ECMWFsst','EOsst','ECMWFssrd',
                                           "sumSun","hrsun","parsum",'logOA')]
names(data_18) <- c('Date','windspeed','winddirection','river','lagriver','rain',
                    'lagrain','ssrd','OA','sst')
names(data_14) <- c('Date','rain','lagrain','river','lagriver','windspeed','winddirection',
                    'ssrd','OA','sst')
data <- rbind(data_14,data_18)
abis <- strptime(data$Date,format="%d/%m/%Y")
data$Date <- as.Date(abis,format="%d-%m-%Y")
data$OA <- NULL
# Read in Count St Austell Data
phyto <- read.xlsx("Biotoxin and Phytoplankton data 2010-2019.xlsx",
                   sheet = 1, detectDates = TRUE,
                   startRow = 2)[-1,] # Read in the phytoplankton
data.phyto <- phyto[,-tail(1:dim(phyto)[2],3)]
phyto <- as.data.frame(phyto)
phyto[,21] <- NULL
phyto <- filter(phyto,Bed.ID=='B70AE') # Filter to St Austell Bay.
phyto <- rename(phyto,count=`Dinophysis.spp.`)%>%filter((!count%in%c("NOT TESTED SUBMITTED
OUTSIDE ROUTINE TESTING FREQUENCY","TEST NOT REQUIRED"))&!is.na(count))
# Call the phytoplankton count 'count'.
phyto_min <- min(as.numeric(phyto$count),na.rm=T)-1
# Replace the 'not detected' values with a random value between 1 and the lowest observed value.
phyto$count[phyto$count=='ND']=sample(1:phyto_min,sum(phyto$count=='ND'),replace=TRUE)

# Read in OA St Austell Data
toxins <- read.xlsx("Biotoxin and Phytoplankton data 2010-2019.xlsx",
                    sheet = 4, detectDates = TRUE,
                    startRow = 1)[-1,] # Read in the toxin data.
toxins <- filter(toxins,Bed.ID=='B70AE') # Filter to St Austell Bay.
toxins <- rename(toxins,OA=X13)%>%filter((!OA%in%c('Test not required',
                                                   'Not tested','Not Tested'))&!is.na(OA))
# Rename the toxin column 'OA'.
toxin_min <- min(as.numeric(toxins$OA),na.rm=T)-1
toxins$OA[toxins$OA=='<RL']=sample(1:toxin_min,sum(toxins$OA=='<RL'),replace=TRUE)
toxins <- rename(toxins,Date=Date.Sample.Collected)
phyto <- rename(phyto, Date=Date.Sample.Collected)
toxins$Date <- as.Date(toxins$Date)
toxins$OA <- as.numeric(toxins$OA)
dfc <- phyto[c('Date','count')]

# Fix messed up dates
x <- strptime(dfc$Date, format="%Y-%m-%d")
x[171:189,] <- phyto[171:189,]$Date.Sample.Arrived.at.Lowestoft
dfc$Date <- as.Date(format(x, "%Y-%m-%d"))
dfc$Time <- as.numeric(dfc$Date)
dfc$count <- as.numeric(dfc$count)

# code for Figure 1. plot
par(mar = c(5, 4, 4, 4) + 0.3)
plot(dfc$Date, dfc$count,type='o', ylab="OA concentration level (100 cells/litre)",
     xlab="Date",main = "Ocadic Acid Toxicity and Concentration Levels",
     cex.lab=1.5, cex.axis=1.5)
par(new = TRUE)
plot(toxins$Date, toxins$OA, type = "o",col='red', axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(toxins$OA)),cex.axis=1.5)
mtext("OA toxicity level (micrograms per liter)", side=4, line=3,cex=1.5)

# code for St Austell Histogram of OA values Figure 2.
hist(toxins$OA, xlab="OA toxicity level (micrograms per liter)",
     main = "OA toxcity distribution",cex.lab=1.5, cex.axis=1.5)

# fit time GAM to 'counts' to interpolate missing values and add them to the data set
m <- gam(count ~ s(Time,k=70), data=dfc,family = Gamma(link='log'))
par(mfrow = c(2, 2))
gam.check(m)
AIC(m)
summary(m)


dates <- seq(as.Date("2011-08-09"), as.Date("2019-12-03"),by='day')
dates <- as.data.frame(as.numeric(dates))
names(dates) <- 'Time'
preds <- gamPredictions <- predict(m, newdata=dates, type="response", se.fit = T)

dates$count <- preds$fit
dates$Time <- as.Date(dates$Time)
names(dates) <- c('Date','count')
plot(dates$Date,dates$count,type='o')

toxins <- toxins[,c('Date','OA')]
x <- strptime(toxins$Date, format="%Y-%m-%d")
toxins$Date <- as.Date(format(x, "%Y-%m-%d"))
toxins$OA <- as.numeric(toxins$OA)

# Preserve original count data
common <- dfc$Date
indx <- c()
for (i in 1:length(common)){
  indx[i] <- which(dates$Date==common[i])
}
dates[indx,]$count <- dfc$count

# megre interpolated counts data and OA toxins data
df_tox <- merge(toxins,dates,by='Date')

# calculate cross-correlation between counts and OA
zz <- ccf(df_tox$count,df_tox$OA)
zz

# Check the difference in weeks between Dates
diffs <- c()
for (i in 1:length(df_tox$Date)){
  diffs[i] <-round(difftime(df_tox$Date[i+1], df_tox$Date[i] , units = c("weeks")),3)
}
diffs <- diffs[1:192]
mean(diffs)
median(diffs)
max(diffs)
min(diffs)

# Create lagged count variables
for (i in 1:2){
  df_tox[[paste0('lag_',i)]] <- sapply(1:nrow(df_tox), function(x) df_tox$count[x-i])
}

#clean up messed up rows
df_tox[1,]$lag_1 <- 0
3
df_tox$lag_1 <- as.numeric(df_tox[,'lag_1'])
df_tox[1:2,]$lag_2 <- 0
df_tox$lag_2 <- as.numeric(df_tox[,'lag_2'])

# check correlation between OA and lagged count
cor(df_tox$OA, df_tox[,c('count','lag_1','lag_2')])

# megre envriomental data and OA toxicity and count data
newdf <- merge(df_tox,data, by='Date')
newdf <- aggregate(.~Date, newdf, median)
newdf$day <- as.numeric(format(newdf$Date, format="%j"))
newdf$time <- as.numeric(newdf$Date)

# Plot of enviromental variables and OA toxicity Figure 3.
par(mfrow = c(2, 4))
for (i in 6:13){
  plot(newdf[,i],newdf$OA,main = names(newdf[i]),xlab=names(newdf[i]),
       ylab = 'OA toxicity',cex=1,cex.lab=2, cex.axis=2)
}


# Split data for testing cnd training
n <- nrow(newdf)-5
xTrain <- newdf[1:n,]
xTest <- newdf[as.numeric(n+1):nrow(newdf),]
dim(xTrain)
dim(xTest)

# fit GAM1 to St Austell Data
m <- gam(OA ~ s(lag_1,k=6) +
           s(day,bs='cc',k=6)+
           s(winddirection,k=6)+
           s(lagriver,k=4),
         data=xTrain,family = Gamma (link = "log"),method = 'REML')
AIC(m)
par(mfrow = c(2,2))
gam.check(m)
summary(m)
plot(m)

# Predict OA toxicity for the next 5 observations
gamPredictions <- predict(m, newdata=xTest, type="response", se.fit = T)
gamPredictions$fit
mpred <- gamPredictions
RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

#Calculate RMSE
RMSE(gamPredictions$fit,xTest$OA)

# Plot predictions against observed with confidence and prediction intervals
# Figure 5. in report
plot(xTest$Date,xTest$OA,ylim=c(0,7000),xlab='Date',ylab='OA toxicity',
     main=' St Austell Bay OA toxicity predictions',pch = 19,
     cex.lab=1.6, cex.axis=1.6, cex=1.6)
points(xTest$Date,gamPredictions$fit,col='red',pch=19,cex=1.6)
lines(xTest$Date,gamPredictions$fit+2*gamPredictions$se.fit,lty=2,col="blue")
lines(xTest$Date,gamPredictions$fit-2*gamPredictions$se.fit,lty=2,col="blue")
abline(h=160)
rmvn <- function(n,mu,sig) { ## function to simulate multivariate-Normal random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n))
}
betas <- rmvn(1000,coef(m),m$Vp) ## 1000 replicate param. vectors
newd <- xTest
Xpred <- predict(m,newd,type="lpmatrix")
LPpreds <- Xpred %*% t(betas)
MUpreds <- apply(LPpreds,1,exp)
shape <- 1/m$scale # and
Ypreds <- MUpreds
for(i in 1:length(xTest$Date)){
  Ypreds[,i] <- rgamma(1000,shape=shape,rate=shape/MUpreds[,i])
}
lines(xTest$Date,apply(Ypreds,2,quantile,probs=0.025),lty=2,col="red")
lines(xTest$Date,apply(Ypreds,2,quantile,probs=0.975),lty=2,col="red")
# Perform a rolling re-train
dfs <- list()
dfs_m <- list()
fcst <- list()
for (i in 30:length(newdf$Date)){
  dfs[[i]] <- newdf[1:i,]
  dfs_m[[i]]<- gam(OA ~ s(lag_2,k=6) +
                     s(day,bs='cc',k=6)+
                     s(winddirection,k=6)+
                     s(lagriver,k=4),
                   data=dfs[[i]],family = Gamma(link='log'))
  future <- as.data.frame(newdf[c(i+1,i+2,i+3,i+4,i+5),])
  names(future) <- names(newdf)
  fcst[[i]] <- predict(dfs_m[[i]], newdata=future, type="response", se.fit = T)
}

# Calculate RMSE from rolling re-train models
RMSE_ <- c()
for (i in 30:length(newdf$Date)){
  RMSE_[i] <- RMSE(fcst[[i]]$fit,dfs[[i+5]]$OA[c(i+1,i+2,i+3,i+4,i+5)])
}

# Plot RMSE socres against No. of observations
# Figure 7 in report
plot(30:length(newdf$Date),RMSE_[30:length(newdf$Date)],type='o',ylim=c(0,8000),
     main='RMSE measures for GAM1', ylab='RMSE',xlab='No. of observations',cex=1.8,
     cex.lab=1.6,cex.axis=2)
RMSE_ <- RMSE_[30:length(RMSE_)]

# False positives flase negatives for farm closures (OA above 160 = 1)
# for last 5 predicted observations
newdf$closure <- ifelse(newdf$OA >= 160,1,0)
confy <- newdf[as.numeric(n+1):nrow(newdf),]
confy$pred <- ifelse(gamPredictions$fit >=160,1,0)
confy$closure <- factor(confy$closure,levels=c(0,1))
confy$pred <- factor(confy$pred,levels=c(0,1))
confusionMatrix(confy$closure,confy$pred)

# Average accuracy for classification of threshold breaches
# based on rolling re-train from 30 data points onwards
confy <- c()
for (i in 30:length(newdf$Date)){
  confy[[i]] <- as.data.frame(newdf$closure[c(i+1,i+2,i+3,i+4,i+5)])
  confy[[i]]$pred <- ifelse(fcst[[i]]$fit >=160,1,0)
  names(confy[[i]]) <-c('y','pred')
  confy[[i]]$y <- factor(confy[[i]]$y,levels=c(0,1))
  confy[[i]]$pred <- factor(confy[[i]]$pred,levels=c(0,1))
}
acc <-c()
for(i in 30:length(newdf$Date)){
  acc[[i]] <- confusionMatrix(confy[[i]]$y,confy[[i]]$pred)$overall[1]
}
acc <- as.data.frame(acc)
acc <- acc[30:as.numeric(nrow(newdf)-1),]
mean(acc)
plot(1:length(acc),acc, type='o')


############## St Austell Bay time series model GAM2 ######################
#####################################################################

# Create time variable
newdf <- newdf[,c('Date','OA')]
newdf$Time <- as.numeric(newdf$Date)

# Split data for training and testing
n <- nrow(newdf)-5
xTrain <- newdf[1:n,]
xTest <- newdf[as.numeric(n+1):nrow(newdf),]
dim(xTrain)
dim(xTest)

# Fit GAM2 for St Austell OA time series
m <- gam(OA ~ s(Time,k=40), data=xTrain,family = Gamma(link='log'))

par(mfrow = c(2, 2))
gam.check(m)
AIC(m)
plot(m)
summary(m)

# Predict next 5 observations
gamPredictions <- predict(m, newdata=xTest, type="response", se.fit = T)
gamPredictions$fit
tpred <- gamPredictions

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

# Calculate RMSE
RMSE(gamPredictions$fit,xTest$OA)

# Figure 6 in report
plot(xTest$Date,xTest$OA,ylim=c(0,7000),xlab='Date',ylab='OA toxicity',
     main=' St Austell Bay OA toxicity predictions',pch = 19,
     cex.lab=1.6, cex.axis=1.6, cex=1.6)
points(xTest$Date,gamPredictions$fit,col='red',pch=19,cex=1.6)
lines(xTest$Date,gamPredictions$fit+2*gamPredictions$se.fit,lty=2,col="blue")
lines(xTest$Date,gamPredictions$fit-2*gamPredictions$se.fit,lty=2,col="blue")
abline(h=160)
lines(xTest$Date,gamPredictions$fit+2*gamPredictions$se.fit,lty=2,col="blue")
lines(xTest$Date,gamPredictions$fit-2*gamPredictions$se.fit,lty=2,col="blue")
abline(h=160)
rmvn <- function(n,mu,sig) { ## function to simulate multivariate-Normal random deviates
  L <- mroot(sig);m <- ncol(L);
  7
  t(mu + L%*%matrix(rnorm(m*n),m,n))
}
betas <- rmvn(1000,coef(m),m$Vp) ## 1000 replicate param. vectors
newd <- xTest
Xpred <- predict(m,newd,type="lpmatrix")
LPpreds <- Xpred %*% t(betas)
MUpreds <- apply(LPpreds,1,exp)
shape <- 1/m$scale # and
Ypreds <- MUpreds
for(i in 1:length(xTest$Date)){
  Ypreds[,i] <- rgamma(1000,shape=shape,rate=shape/MUpreds[,i])
}
lines(xTest$Date,apply(Ypreds,2,quantile,probs=0.025),lty=2,col="red")
lines(xTest$Date,apply(Ypreds,2,quantile,probs=0.975),lty=2,col="red")


# GAM2 Rolling re-train from 50 observations onwards
dfs <- list()
dfs_m <- list()
fcst <- list()
for (i in 50:length(newdf$Date)){
  dfs[[i]] <- newdf[1:i,]
  dfs_m[[i]]<- gam(OA ~ s(Time,k=40), data=dfs[[i]],family = Gamma(link='log'))
  future <- as.data.frame(newdf$Time[c(i+1,i+2,i+3,i+4,i+5)])
  names(future) <- 'Time'
  fcst[[i]] <- predict(dfs_m[[i]], newdata=future, type="response", se.fit = T)
}

# Plots of some of the models during re-trains
seqm <- seq(50,length(newdf$Date), 5)
par(mfrow = c(2, 3))
for (i in seqm){
  plot(dfs[[i+5]]$Time[c(i+1,i+2,i+3,i+4,i+5)],fcst[[i]]$fit,type='o',col='red',ylim=c(0,5000))
  lines(dfs[[i+5]]$Time[c(i+1,i+2,i+3,i+4,i+5)],dfs[[i+5]]$OA[c(i+1,i+2,i+3,i+4,i+5)], type='o')
}

# Calculate RMSe based on rolling re-train
RMSE_T <- c()
for (i in 50:length(newdf$Date)){
  RMSE_T[i] <- RMSE(fcst[[i]]$fit,dfs[[i+5]]$OA[c(i+1,i+2,i+3,i+4,i+5)])
}

# Figure 7 in report
plot(50:length(newdf$Date),RMSE_T[50:length(newdf$Date)],type='o',ylim=c(0,8000),
     main='RMSE measures for GAM2', ylab='RMSE',xlab='No. of observations',cex=1.8,
     cex.lab=1.6,cex.axis=2)


# Accuracy socres for true positives and true negatives
newdf$closure <- ifelse(newdf$OA >= 160,1,0)
confy <- newdf[as.numeric(n+1):nrow(newdf),]
confy$pred <- ifelse(gamPredictions$fit >=160,1,0)
confy$closure <- factor(confy$closure,levels=c(0,1))
confy$pred <- factor(confy$pred,levels=c(0,1))
confusionMatrix(confy$closure,confy$pred)


# Avergaed accuracy socres for classification of threshold breaches
# based on the rolling re-trains
confy <- c()
for (i in 50:length(newdf$Date)){
  confy[[i]] <- as.data.frame(newdf$closure[c(i+1,i+2,i+3,i+4,i+5)])
  confy[[i]]$pred <- ifelse(fcst[[i]]$fit >=160,1,0)
  names(confy[[i]]) <-c('y','pred')
  confy[[i]]$y <- factor(confy[[i]]$y,levels=c(0,1))
  confy[[i]]$pred <- factor(confy[[i]]$pred,levels=c(0,1))
}
acc <-c()
for(i in 50:length(newdf$Date)){
  acc[[i]] <- confusionMatrix(confy[[i]]$y,confy[[i]]$pred)$overall[1]
}
acc <- as.data.frame(acc)
acc <- acc[50:as.numeric(nrow(newdf)-1),]
mean(acc)


### Averaged predictions
mpred
tpred
avg_pred <- c()
for(i in 1:5){
  avg_pred[i] <- (2/3)*tpred$fit[i] + (1/3)*mpred$fit[i]
}
RMSE(avg_pred,xTest$OA)
Code for Lyme Bay - Same methodology as for St Austell Bay code with minor data manipulation changes
set.seed(1433)

####################### CODE FOR LYME BAY ######################
################################################################

# Read in Count Lyme Bay
phyto <- read.xlsx("Biotoxin and Phytoplankton data 2010-2019.xlsx",
                   sheet = 1, detectDates = TRUE,
                   9
                   startRow = 2)[-1,] # Read in the phytoplankton
data.phyto <- phyto[,-tail(1:dim(phyto)[2],3)]
phyto <- as.data.frame(phyto)
phyto[,21] <- NULL
phyto <- filter(phyto,Production.Area=='Lyme Bay') # Filter to St Austell Bay.
phyto <- rename(phyto,count=`Dinophysis.spp.`)%>%filter((!count%in%c("NOT TESTED SUBMITTED
OUTSIDE ROUTINE TESTING FREQUENCY","TEST NOT REQUIRED"))&!is.na(count))
# Call the phytoplankton count 'count'.
phyto_min <- min(as.numeric(phyto$count),na.rm=T)-1
#Replace the 'not detected' values with a random value between 1 and the lowest observed value.
phyto$count[phyto$count=='ND']=sample(1:phyto_min,sum(phyto$count=='ND'),replace=TRUE)
# Read in OA Lyme Bay Data
toxins <- read.xlsx("Biotoxin and Phytoplankton data 2010-2019.xlsx",
                    sheet = 4, detectDates = TRUE,
                    startRow = 1)[-1,] # Read in the toxin data.
toxins <- filter(toxins,Production.Area=='Lyme Bay') # Filter to St Austell Bay.
toxins <- rename(toxins,OA=X13)%>%filter((!OA%in%c('Test not required',
                                                   'Not tested','Not Tested'))&!is.na(OA))
# Rename the toxin column 'OA'.
toxin_min <- min(as.numeric(toxins$OA),na.rm=T)-1
toxins$OA[toxins$OA=='<RL']=sample(1:toxin_min,sum(toxins$OA=='<RL'),replace=TRUE)
toxins <- rename(toxins,Date=Date.Sample.Collected)
phyto <- rename(phyto, Date=Date.Sample.Collected)

toxins$OA <- as.numeric(toxins$OA)
toxins$Date <- as.Date(toxins$Date)
dfc <- phyto[,c('Date','count')]

x <- strptime(dfc$Date, format="%Y-%m-%d")
x[81:105,] <- phyto[81:105,]$Date.Sample.Arrived.at.Lowestoft
dfc$Date <- as.Date(format(x, "%Y-%m-%d"))
dfc$Time <- as.numeric(dfc$Date)
dfc$count <- as.numeric(dfc$count)

# remove two observations from 2010
dfc <- dfc[3:105,]

# code for Figure 10
hist(toxins$OA, xlab="OA toxicity level (micrograms per liter)",
     main = "OA toxcity distribution",cex.lab=1.7, cex.axis=1.7)
# code for Figure 12
par(mar = c(5, 4, 4, 4) + 0.3)
plot(dfc$Date, dfc$count,type='o', ylab="OA concentration level (100 cells/litre)",
     xlab="Date",main = "Ocadic Acid Toxicity and Concentration Levels",
     cex.lab=1.5, cex.axis=1.7)
par(new = TRUE)
plot(toxins$Date, toxins$OA, type = "o",col='red', axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(toxins$OA)),cex.axis=1.5)
mtext("OA toxicity level (micrograms per liter)", side=4, line=3,cex=1.5)


# fit time GAM to counts to interpolate missing values
m <- gam(count ~ s(Time,k=70), data=dfc,family = Gamma(link='log'))

par(mfrow = c(2, 2))
gam.check(m)
AIC(m)
plot(m)
summary(m)

dates <- seq(as.Date("2015-07-22"), as.Date("2019-12-05"),by='day')
dates <- as.data.frame(as.numeric(dates))
names(dates) <- 'Time'
preds <- gamPredictions <- predict(m, newdata=dates, type="response", se.fit = T)

dates$count <- preds$fit
dates$Time <- as.Date(dates$Time)
names(dates) <- c('Date','count')
plot(dates$Date,dates$count,type='o')

toxins <- toxins[,c('Date','OA')]
x <- strptime(toxins$Date, format="%Y-%m-%d")
toxins$Date <- as.Date(format(x, "%Y-%m-%d"))
toxins$OA <- as.numeric(toxins$OA)

common <- dfc$Date
indx <- c()
for (i in 1:length(common)){
  indx[i] <- which(dates$Date==common[i])
}
dates[indx,]$count <- dfc$count
df_tox <- merge(toxins,dates,by='Date')

zz <- ccf(df_tox$count,df_tox$OA)
zz

# Check the difference in weeks between Dates
diffs <- c()
for (i in 1:length(df_tox$Date)){
  diffs[i] <-round(difftime(df_tox$Date[i+1], df_tox$Date[i] , units = c("weeks")),2)
}

diffs <- diffs[1:88]
mean(diffs)
median(diffs)
max(diffs)
min(diffs)

# Create lagged count variable
for (i in 1:2){
  df_tox[[paste0('lag_',i)]] <- sapply(1:nrow(df_tox), function(x) df_tox$count[x-i])
}

#remove messed up rows
df_tox[1,]$lag_1 <- 0
df_tox$lag_1 <- as.numeric(df_tox[,'lag_1'])
df_tox[1:2,]$lag_2 <- 0
df_tox$lag_2 <- as.numeric(df_tox[,'lag_2'])

# check correlation between OA and lagged counts
cor(df_tox$OA, df_tox[,c('count','lag_1','lag_2')])

# cretae day of year variable
df_tox$day <- as.numeric(format(df_tox$Date,format="%j"))
newdf <- df_tox

# Plot time series of OA toxicity
plot(newdf$Date,newdf$OA,type='o',col='red',main='Lyme Bay')
lines(newdf$Date,newdf$count,type='o')

# Split data for testing cnd training
n <- nrow(newdf)-5
xTrain <- newdf[1:n,]
xTest <- newdf[as.numeric(n+1):nrow(newdf),]
dim(xTrain)
dim(xTest)


# Fit GAM1 model to Lyme Bay data
m <- gam(OA ~ s(lag_1,k=6) +
           s(day,bs='cc',k=6),
         data=xTrain,family = Gamma (link = "log"),method = 'REML')
AIC(m)
par(mfrow = c(2,2))
gam.check(m)
summary(m)
plot(m)

gamPredictions <- predict(m, newdata=xTest, type="response", se.fit = T)
gamPredictions$fit
mpred <- gamPredictions

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

RMSE(gamPredictions$fit,xTest$OA)

# code for Figure 14
plot(xTest$Date,xTest$OA,ylim=c(0,800),xlab='Date',ylab='OA toxicity',
     main=' St Austell Bay OA toxicity predictions',pch = 19,
     cex.lab=2, cex.axis=2, cex=1.6)
points(xTest$Date,gamPredictions$fit,col='red',pch=19,cex=1.6)
lines(xTest$Date,gamPredictions$fit+2*gamPredictions$se.fit,lty=2,col="blue")
lines(xTest$Date,gamPredictions$fit-2*gamPredictions$se.fit,lty=2,col="blue")
abline(h=160)
rmvn <- function(n,mu,sig) { ## function to simulate multivariate-Normal random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n))
}
betas <- rmvn(1000,coef(m),m$Vp) ## 1000 replicate param. vectors
newd <- xTest
Xpred <- predict(m,newd,type="lpmatrix")
LPpreds <- Xpred %*% t(betas)
MUpreds <- apply(LPpreds,1,exp)
shape <- 1/m$scale # and
Ypreds <- MUpreds
for(i in 1:length(xTest$Date)){
  Ypreds[,i] <- rgamma(1000,shape=shape,rate=shape/MUpreds[,i])
}
lines(xTest$Date,apply(Ypreds,2,quantile,probs=0.025),lty=2,col="red")
lines(xTest$Date,apply(Ypreds,2,quantile,probs=0.975),lty=2,col="red")

# Rolling re-train
dfs <- list()
dfs_m <- list()
fcst <- list()
for (i in 30:length(newdf$Date)){
  dfs[[i]] <- newdf[1:i,]
  dfs_m[[i]]<- gam(OA ~ s(lag_2,k=6) + s(day,bs='cc',k=6),
                   data=dfs[[i]],family = Gamma(link='log'))
  future <- as.data.frame(newdf[c(i+1,i+2,i+3,i+4,i+5),])
  names(future) <- names(newdf)
  fcst[[i]] <- predict(dfs_m[[i]], newdata=future, type="response", se.fit = T)
}


seqm <- seq(30, length(newdf$Date), 5)
par(mfrow = c(2, 3))
for (i in seqm){
  plot(dfs[[i+5]]$Date[c(i+1,i+2,i+3,i+4,i+5)],fcst[[i]]$fit,type='o',col='red',ylim=c(0,1000))
  lines(dfs[[i+5]]$Date[c(i+1,i+2,i+3,i+4,i+5)],dfs[[i+5]]$OA[c(i+1,i+2,i+3,i+4,i+5)], type='o')
}


RMSE_ <- c()
for (i in 30:length(newdf$Date)){
  RMSE_[i] <- RMSE(fcst[[i]]$fit,dfs[[i+5]]$OA[c(i+1,i+2,i+3,i+4,i+5)])
}


plot(30:length(newdf$Date),RMSE_[30:length(newdf$Date)],type='o')


# False positives flase negatives for farm closures (OA above 160 = 1)
# for last 5 predicted observations
newdf$closure <- ifelse(newdf$OA >=160, 1, 0)
confy <- newdf[as.numeric(n+1):nrow(newdf),]
confy$pred <- ifelse(gamPredictions$fit >=160,1,0)
confy$closure <- factor(confy$closure,levels=c(0,1))
confy$pred <- factor(confy$pred,levels=c(0,1))
confusionMatrix(confy$closure,confy$pred)


# Average accuracy from 50 data points onwards
dfs <- list()
dfs_m <- list()
fcst <- list()
for (i in 30:length(newdf$Date)){
  dfs[[i]] <- newdf[1:i,]
  dfs_m[[i]]<- gam(OA ~ s(lag_2,k=6) + s(day,bs='cc',k=6),
                   data=dfs[[i]],family = Gamma(link='log'))
  future <- as.data.frame(newdf[c(i+1,i+2,i+3,i+4,i+5),])
  names(future) <- names(newdf)
  fcst[[i]] <- predict(dfs_m[[i]], newdata=future, type="response", se.fit = T)
}


confy <- c()
for (i in 30:length(newdf$Date)){
  confy[[i]] <- as.data.frame(newdf$closure[c(i+1,i+2,i+3,i+4,i+5)])
  confy[[i]]$pred <- ifelse(fcst[[i]]$fit >=160,1,0)
  names(confy[[i]]) <-c('y','pred')
  14
  confy[[i]]$y <- factor(confy[[i]]$y,levels=c(0,1))
  confy[[i]]$pred <- factor(confy[[i]]$pred,levels=c(0,1))
}


acc <-c()
for(i in 30:length(newdf$Date)){
  acc[[i]] <- confusionMatrix(confy[[i]]$y,confy[[i]]$pred)$overall[1]
}

acc <- as.data.frame(acc)
acc <- acc[30:as.numeric(nrow(newdf)-1),]
mean(acc)

############## Lyme Bay time series model GAM2 ######################
#####################################################################
toxins_ts <- read.xlsx("Biotoxin and Phytoplankton data 2010-2019.xlsx",
                       sheet = 4, detectDates = TRUE,
                       startRow = 1)[-1,] # Read in the toxin data.
toxins_ts <- filter(toxins_ts,Production.Area=='Lyme Bay') # Filter to St Austell Bay.
toxins_ts <- rename(toxins_ts,OA=X13)%>%filter((!OA%in%c('Test not required',
                                                         'Not tested','Not Tested'))&!is.na(OA))
# Rename the toxin column 'OA'.
toxin_min <- min(as.numeric(toxins_ts$OA),na.rm=T)-1
toxins_ts$OA[toxins_ts$OA=='<RL']=sample(1:toxin_min,sum(toxins_ts$OA=='<RL'),replace=TRUE)
toxins_ts <- rename(toxins_ts,Date=Date.Sample.Collected)
toxins_ts <- toxins_ts[,c(7,13)]
toxins_ts$OA <- as.numeric(toxins_ts$OA)

x <- strptime(toxins_ts$Date, format="%Y-%m-%d")
toxins_ts$Date <- as.Date(format(x, "%Y-%m-%d"))

df <- toxins_ts
df$Time <- as.numeric(df$Date)
df$day <- as.numeric(format(df$Date, format="%j"))


newdf <- df
n <- nrow(newdf)-5
xTrain <- newdf[1:n,]
xTest <- newdf[as.numeric(n+1):nrow(newdf),]
dim(xTrain)
dim(xTest)


# Fit GAM to time
m <- gam(OA ~ s(Time,k=45), data=xTrain,family = Gamma(link='log'))

par(mfrow = c(2, 2))
gam.check(m)
AIC(m)
plot(m)
summary(m)


gamPredictions <- predict(m, newdata=xTest, type="response", se.fit = T)
gamPredictions$fit
tpred <- gamPredictions

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

RMSE(gamPredictions$fit,xTest$OA)


# code for Figure 14
plot(xTest$Date,xTest$OA,ylim=c(0,800),xlab='Date',ylab='OA toxicity',
     main=' St Austell Bay OA toxicity predictions',pch = 19,
     cex.lab=2, cex.axis=2, cex=1.6)
points(xTest$Date,gamPredictions$fit,col='red',pch=19,cex=1.6)
lines(xTest$Date,gamPredictions$fit+2*gamPredictions$se.fit,lty=2,col="blue")
lines(xTest$Date,gamPredictions$fit-2*gamPredictions$se.fit,lty=2,col="blue")
abline(h=160)
lines(xTest$Date,gamPredictions$fit+2*gamPredictions$se.fit,lty=2,col="blue")
lines(xTest$Date,gamPredictions$fit-2*gamPredictions$se.fit,lty=2,col="blue")
abline(h=160)
rmvn <- function(n,mu,sig) { ## function to simulate multivariate-Normal random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n))
}
betas <- rmvn(1000,coef(m),m$Vp) ## 1000 replicate param. vectors
newd <- xTest
Xpred <- predict(m,newd,type="lpmatrix")
LPpreds <- Xpred %*% t(betas)
MUpreds <- apply(LPpreds,1,exp)
shape <- 1/m$scale # and
Ypreds <- MUpreds
for(i in 1:length(xTest$Date)){
  Ypreds[,i] <- rgamma(1000,shape=shape,rate=shape/MUpreds[,i])
}
lines(xTest$Date,apply(Ypreds,2,quantile,probs=0.025),lty=2,col="red")
lines(xTest$Date,apply(Ypreds,2,quantile,probs=0.975),lty=2,col="red")


# Rolling re-train
dfs <- list()
dfs_m <- list()
fcst <- list()
for (i in 50:length(newdf$Date)){
  dfs[[i]] <- df[1:i,]
  dfs_m[[i]]<- gam(OA ~ s(Time,k=40), data=dfs[[i]],family = Gamma(link='log'))
  future <- as.data.frame(df$Time[c(i+1,i+2,i+3,i+4,i+5)])
  names(future) <- 'Time'
  fcst[[i]] <- predict(dfs_m[[i]], newdata=future, type="response", se.fit = T)
}

seqm <- seq(50,length(newdf$Date), 5)
par(mfrow = c(2, 3))
for (i in seqm){
  plot(dfs[[i+5]]$Time[c(i+1,i+2,i+3,i+4,i+5)],fcst[[i]]$fit,type='o',col='red',ylim=c(0,5000))
  lines(dfs[[i+5]]$Time[c(i+1,i+2,i+3,i+4,i+5)],dfs[[i+5]]$OA[c(i+1,i+2,i+3,i+4,i+5)], type='o')
}

RMSE_T <- c()
for (i in 50:length(newdf$Date)){
  RMSE_T[i] <- RMSE(fcst[[i]]$fit,dfs[[i+5]]$OA[c(i+1,i+2,i+3,i+4,i+5)])
}

plot(50:90,RMSE_T[50:90],type='o')


# Accuracy socres for true positives and true negatives
newdf$closure <- ifelse(newdf$OA >= 160,1,0)
confy <- newdf[as.numeric(n+1):nrow(newdf),]
confy$pred <- ifelse(gamPredictions$fit >=160,1,0)
confy$closure <- factor(confy$closure,levels=c(0,1))
confy$pred <- factor(confy$pred,levels=c(0,1))
confusionMatrix(confy$closure,confy$pred)


# Avergaed accuracy socres for classification of threshold breaches
dfs <- list()
dfs_m <- list()
fcst <- list()
for (i in 50:length(newdf$Date)){
  dfs[[i]] <- newdf[1:i,]
  dfs_m[[i]]<- dfs_m[[i]]<- gam(OA ~ s(Time,k=40), data=dfs[[i]],family = Gamma(link='log'))
  17
  future <- as.data.frame(newdf[c(i+1,i+2,i+3,i+4,i+5),])
  names(future) <- names(newdf)
  fcst[[i]] <- predict(dfs_m[[i]], newdata=future, type="response", se.fit = T)
}


confy <- c()
for (i in 50:length(newdf$Date)){
  confy[[i]] <- as.data.frame(newdf$closure[c(i+1,i+2,i+3,i+4,i+5)])
  confy[[i]]$pred <- ifelse(fcst[[i]]$fit >=160,1,0)
  names(confy[[i]]) <-c('y','pred')
  confy[[i]]$y <- factor(confy[[i]]$y,levels=c(0,1))
  confy[[i]]$pred <- factor(confy[[i]]$pred,levels=c(0,1))
}


acc <-c()
for(i in 50:length(newdf$Date)){
  acc[[i]] <- confusionMatrix(confy[[i]]$y,confy[[i]]$pred)$overall[1]
}
acc <- as.data.frame(acc)
acc <- acc[50:as.numeric(nrow(newdf)-1),]
mean(acc)
Code for Lantivet Bay - Same methodology as for St Austell Bay code but with minor data manipulation
changes
set.seed(1433)

############## CODE FOR LANTIVET BAY ########################
###########################################################
# Read in Count Lantivet Bay
phyto <- read.xlsx("Biotoxin and Phytoplankton data 2010-2019.xlsx",
                   sheet = 1, detectDates = TRUE,
                   startRow = 2)[-1,] # Read in the phytoplankton
data.phyto <- phyto[,-tail(1:dim(phyto)[2],3)]
phyto <- as.data.frame(phyto)
phyto[,21] <- NULL
phyto <- filter(phyto,Production.Area=='Lantivet Bay') # Filter to St Austell Bay.
phyto <- rename(phyto,count=`Dinophysis.spp.`)%>%filter((!count%in%c("NOT TESTED SUBMITTED
OUTSIDE ROUTINE TESTING FREQUENCY","TEST NOT REQUIRED"))&!is.na(count))
# Call the phytoplankton count 'count'.
phyto_min <- min(as.numeric(phyto$count),na.rm=T)-1
phyto$count[phyto$count=='ND']=sample(1:phyto_min,sum(phyto$count=='ND'),replace=TRUE) # Replace the 'not # a random value between 1 and the lowest observed value.
# Read in OA Lantivet Bay Data
toxins <- read.xlsx("Biotoxin and Phytoplankton data 2010-2019.xlsx",
                    sheet = 4, detectDates = TRUE,
                    startRow = 1)[-1,] # Read in the toxin data.
toxins <- filter(toxins,Production.Area=='Lantivet Bay') # Filter to St Austell Bay.
toxins <- rename(toxins,OA=X13)%>%filter((!OA%in%c('Test not required',
                                                   'Not tested','Not Tested'))&!is.na(OA))
# Rename the toxin column 'OA'.
toxin_min <- min(as.numeric(toxins$OA),na.rm=T)-1
toxins$OA[toxins$OA=='<RL']=sample(1:toxin_min,sum(toxins$OA=='<RL'),replace=TRUE)

# Rename Date columns
toxins <- rename(toxins,Date=Date.Sample.Collected)
phyto <- rename(phyto, Date=Date.Sample.Collected)
toxins$OA <- as.numeric(toxins$OA)
toxins$Date <- as.Date(toxins$Date)
toxins <- toxins[c('Date','OA')]
toxins <- na.omit(toxins)

#Select relevant columns
dfc <- phyto[c('Date','count')]

# Fix Date formats for currupted dates
x <- strptime(dfc$Date, format="%Y-%m-%d")
x[83:99,] <- phyto[83:99,]$Date.Sample.Arrived.at.Lowestoft
dfc$Date <- as.Date(format(x, "%Y-%m-%d"))

# Create time variable and make numeric
dfc$Time <- as.numeric(dfc$Date)
dfc$count <- as.numeric(dfc$count)

# Histogram for Lantivet Bay Figure 11
hist(toxins$OA, xlab="OA toxicity level (micrograms per liter)",
     main = "OA toxcity distribution",cex.lab=1.7, cex.axis=1.7)

# code for Figure 13
par(mar = c(5, 4, 4, 4) + 0.3)
plot(dfc$Date, dfc$count,type='o', ylab="OA concentration level (100 cells/litre)",
     xlab="Date",main = "Ocadic Acid Toxicity and Concentration Levels",
     cex.lab=1.5, cex.axis=2)
par(new = TRUE)
plot(toxins$Date, toxins$OA, type = "o",col='red', axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(toxins$OA)),cex.axis=1.5)
mtext("OA toxicity level (micrograms per liter)", side=4, line=3,cex=1.5)


# fit time GAM to counts data to interpolate missing values
m <- gam(count ~ s(Time,k=70), data=dfc,family = Gamma(link='log'))

par(mfrow = c(2, 2))
gam.check(m)
AIC(m)
plot(m)
summary(m)

# Create date sequence starting from the first date of the dfc data set
dates <- seq(as.Date("2015-04-29"), as.Date("2019-12-03"),by='day')
dates <- as.data.frame(as.numeric(dates))
names(dates) <- 'Time'

# Interpolate count data
preds <- gamPredictions <- predict(m, newdata=dates, type="response", se.fit = T)
plot(dates$Time,preds$fit,type='o')
dates$count <- preds$fit
dates$Time <- as.Date(dates$Time)
names(dates) <- c('Date','count')
plot(dates$Date,dates$count,type='o')

# Select relevant columns from OA toxicity data set
toxins <- toxins[,c('Date','OA')]

# Make dates as.Date format
x <- strptime(toxins$Date, format="%Y-%m-%d")
toxins$Date <- as.Date(format(x, "%Y-%m-%d"))

# Make OA as numeric format
toxins$OA <- as.numeric(toxins$OA)

# Preserve count data from the dfc dataset
common <- dfc$Date
indx <- c()
for (i in 1:length(common)){
  indx[i] <- which(dates$Date==common[i])
}

# Replace modelled count data with available real count data
dates[indx,]$count <- dfc$count

# Merge OA toxins and dates(with counts) dataframes by Date
df_tox <- merge(toxins,dates,by='Date')
df_tox <-na.omit(df_tox)

# Calculate cross-correlation of counts and OA
zz <- ccf(df_tox$count,df_tox$OA)
zz

#Check the difference in weeks between Dates
diffs <- c()
for (i in 1:length(df_tox$Date)){
  diffs[i] <-round(difftime(df_tox$Date[i+1], df_tox$Date[i] , units = c("weeks")),2)
}
diffs <- diffs[1:119]
mean(diffs)
median(diffs)
max(diffs)
min(diffs)

# Create lagged count variables for 2 lags
for (i in 1:2){
  df_tox[[paste0('lag_',i)]] <- sapply(1:nrow(df_tox), function(x) df_tox$count[x-i])
}
#remove messed up rows
df_tox[1,]$lag_1 <- 0
df_tox$lag_1 <- as.numeric(df_tox[,'lag_1'])
df_tox[1:2,]$lag_2 <- 0
df_tox$lag_2 <- as.numeric(df_tox[,'lag_2'])

# Check correlation between OA and lagged counts
cor(df_tox$OA, df_tox[,c('count','lag_1','lag_2')])

# Create day of year variable
df_tox$day <- as.numeric(format(df_tox$Date,format="%j"))
newdf <- df_tox

# Plot time series of OA toxicity and count
plot(newdf$Date,newdf$OA,type='o',col='red')
lines(newdf$Date,newdf$count,type='o')

# Split data for testing cnd training leaving 5 last observations for testing
n <- nrow(newdf)-5
xTrain <- newdf[1:n,]
xTest <- newdf[as.numeric(n+1):nrow(newdf),]
dim(xTrain)
dim(xTest)

# GAM1 for Lantivet Bay as described in the report
m <- gam(OA ~ s(lag_1,k=6) +
           s(day,bs='cc',k=6),
         data=xTrain,family = Gamma (link = "log"),method = 'REML')
AIC(m)
par(mfrow = c(2,2))
gam.check(m)
summary(m)
plot(m)

# Predict values
gamPredictions <- predict(m, newdata=xTest, type="response", se.fit = T)
gamPredictions$fit
mpred <- gamPredictions$fit

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

RMSE(gamPredictions$fit,xTest$OA)

# Plot predictions against observations Figure 15
plot(xTest$Date,xTest$OA,ylim=c(0,2000),xlab='Date',ylab='OA toxicity',
     main=' St Austell Bay OA toxicity predictions',pch = 19,
     cex.lab=2, cex.axis=2, cex=1.6)
points(xTest$Date,gamPredictions$fit,col='red',pch=19,cex=1.6)
lines(xTest$Date,gamPredictions$fit+2*gamPredictions$se.fit,lty=2,col="blue")
lines(xTest$Date,gamPredictions$fit-2*gamPredictions$se.fit,lty=2,col="blue")
abline(h=160)
rmvn <- function(n,mu,sig) { ## function to simulate multivariate-Normal random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n))
}
betas <- rmvn(1000,coef(m),m$Vp) ## 1000 replicate param. vectors
newd <- xTest
Xpred <- predict(m,newd,type="lpmatrix")
LPpreds <- Xpred %*% t(betas)
MUpreds <- apply(LPpreds,1,exp)
shape <- 1/m$scale # and
Ypreds <- MUpreds
for(i in 1:length(xTest$Date)){
  Ypreds[,i] <- rgamma(1000,shape=shape,rate=shape/MUpreds[,i])
}
# Plot confidence intervals
lines(xTest$Date,apply(Ypreds,2,quantile,probs=0.025),lty=2,col="red")
lines(xTest$Date,apply(Ypreds,2,quantile,probs=0.975),lty=2,col="red")



# Rolling re-train
dfs <- list()
dfs_m <- list()
fcst <- list()
for (i in 30:length(newdf$Date)){
  dfs[[i]] <- newdf[1:i,]
  dfs_m[[i]]<- gam(OA ~ s(lag_2,k=6) + s(day,bs='cc',k=6),
                   data=dfs[[i]],family = Gamma(link='log'))
  future <- as.data.frame(newdf[c(i+1,i+2,i+3,i+4,i+5),])
  names(future) <- names(newdf)
  fcst[[i]] <- predict(dfs_m[[i]], newdata=future, type="response", se.fit = T)
}


seqm <- seq(30, length(newdf$Date), 3)
par(mfrow = c(2, 3))
for (i in seqm){
  plot(dfs[[i+5]]$Date[c(i+1,i+2,i+3,i+4,i+5)],fcst[[i]]$fit,type='o',col='red',ylim=c(0,5000))
  lines(dfs[[i+5]]$Date[c(i+1,i+2,i+3,i+4,i+5)],dfs[[i+5]]$OA[c(i+1,i+2,i+3,i+4,i+5)], type='o')
}

RMSE_ <- c()
for (i in 30:length(newdf$Date)){
  RMSE_[i] <- RMSE(fcst[[i]]$fit,dfs[[i+5]]$OA[c(i+1,i+2,i+3,i+4,i+5)])
}

plot(30:length(newdf$Date),RMSE_[30:length(newdf$Date)],type='o')


# False positives flase negatives for farm closures (OA above 160 = 1)
# for last 5 predicted observations
newdf$closure <- ifelse(newdf$OA >= 160,1,0)
confy <- newdf[as.numeric(n+1):nrow(newdf),]
confy$pred <- ifelse(gamPredictions$fit >=160,1,0)
confy$closure <- factor(confy$closure,levels=c(0,1))
confy$pred <- factor(confy$pred,levels=c(0,1))
confusionMatrix(confy$closure,confy$pred)


# Average accuracy from 50 data points onwards
confy <- c()
for (i in 30:length(newdf$Date)){
  confy[[i]] <- as.data.frame(newdf$closure[c(i+1,i+2,i+3,i+4,i+5)])
  confy[[i]]$pred <- ifelse(fcst[[i]]$fit >=160,1,0)
  names(confy[[i]]) <-c('y','pred')
  confy[[i]]$y <- factor(confy[[i]]$y,levels=c(0,1))
  confy[[i]]$pred <- factor(confy[[i]]$pred,levels=c(0,1))
}


acc <-c()
for(i in 30:length(newdf$Date)){
  acc[[i]] <- confusionMatrix(confy[[i]]$y,confy[[i]]$pred)$overall[1]
}
acc <- as.data.frame(acc)
acc <- acc[30:as.numeric(nrow(newdf)-1),]
mean(acc)
plot(1:90,acc, type='o')


############## Lantivet Bay time series model GAM2 ######################
####################################################################
toxins_ts <- read.xlsx("Biotoxin and Phytoplankton data 2010-2019.xlsx",
                       sheet = 4, detectDates = TRUE,
                       startRow = 1)[-1,] # Read in the toxin data.
toxins_ts <- filter(toxins_ts,Production.Area=='Lantivet Bay') # Filter to St Austell Bay.
toxins_ts <- rename(toxins_ts,OA=X13)%>%filter((!OA%in%c('Test not required',
                                                         'Not tested','Not Tested'))&!is.na(OA))
# Rename the toxin column 'OA'.
toxin_min <- min(as.numeric(toxins_ts$OA),na.rm=T)-1
toxins_ts$OA[toxins_ts$OA=='<RL']=sample(1:toxin_min,sum(toxins_ts$OA=='<RL'),replace=TRUE)
toxins_ts <- rename(toxins_ts,Date=Date.Sample.Collected)
toxins_ts <- toxins_ts[,c(7,13)]
toxins_ts$OA <- as.numeric(toxins_ts$OA)

x <- strptime(toxins_ts$Date, format="%Y-%m-%d")
toxins_ts$Date <- as.Date(format(x, "%Y-%m-%d"))
df <- toxins_ts
df$Time <- as.numeric(df$Date)
df$day <- as.numeric(format(df$Date, format="%j"))
df <-na.omit(df)


newdf <- df
n <- nrow(newdf)-5
xTrain <- newdf[1:n,]
xTest <- newdf[as.numeric(n+1):nrow(newdf),]
dim(xTrain)
dim(xTest)


# Fit GAM to time
m <- gam(OA ~ s(Time,k=40), data=xTrain,family = Gamma(link='log'))


par(mfrow = c(2, 2))
gam.check(m)
AIC(m)
plot(m)
summary(m)

gamPredictions <- predict(m, newdata=xTest, type="response", se.fit = T)
gamPredictions$fit
tpred <- gamPredictions

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}
RMSE(gamPredictions$fit,xTest$OA)


# Figure 15
plot(xTest$Date,xTest$OA,ylim=c(0,2000),xlab='Date',ylab='OA toxicity',
     main=' St Austell Bay OA toxicity predictions',pch = 19,
     cex.lab=2.2, cex.axis=2.2, cex=1.6)
points(xTest$Date,gamPredictions$fit,col='red',pch=19,cex=1.6)
lines(xTest$Date,gamPredictions$fit+2*gamPredictions$se.fit,lty=2,col="blue")
lines(xTest$Date,gamPredictions$fit-2*gamPredictions$se.fit,lty=2,col="blue")
abline(h=160)
rmvn <- function(n,mu,sig) { ## function to simulate multivariate-Normal random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n))
}
betas <- rmvn(1000,coef(m),m$Vp) ## 1000 replicate param. vectors
newd <- xTest
Xpred <- predict(m,newd,type="lpmatrix")
LPpreds <- Xpred %*% t(betas)
MUpreds <- apply(LPpreds,1,exp)
shape <- 1/m$scale # and
Ypreds <- MUpreds
for(i in 1:length(xTest$Date)){
  Ypreds[,i] <- rgamma(1000,shape=shape,rate=shape/MUpreds[,i])
}
lines(xTest$Date,apply(Ypreds,2,quantile,probs=0.025),lty=2,col="red")
lines(xTest$Date,apply(Ypreds,2,quantile,probs=0.975),lty=2,col="red")


# GAM Retrain
dfs <- list()
dfs_m <- list()
fcst <- list()
for (i in 50:length(newdf$Date)){
  dfs[[i]] <- df[1:i,]
  dfs_m[[i]]<- gam(OA ~ s(Time,k=40), data=dfs[[i]],family = Gamma(link='log'))
  future <- as.data.frame(df$Time[c(i+1,i+2,i+3,i+4,i+5)])
  names(future) <- 'Time'
  fcst[[i]] <- predict(dfs_m[[i]], newdata=future, type="response", se.fit = T)
}

seqm <- seq(50,length(newdf$Date), 3)
par(mfrow = c(2, 3))
for (i in seqm){
  plot(dfs[[i+5]]$Time[c(i+1,i+2,i+3,i+4,i+5)],fcst[[i]]$fit,type='o',col='red',ylim=c(0,5000))
  lines(dfs[[i+5]]$Time[c(i+1,i+2,i+3,i+4,i+5)],dfs[[i+5]]$OA[c(i+1,i+2,i+3,i+4,i+5)], type='o')
}

RMSE_T <- c()
for (i in 50:length(newdf$Date)){
  RMSE_T[i] <- RMSE(fcst[[i]]$fit,dfs[[i+5]]$OA[c(i+1,i+2,i+3,i+4,i+5)])
}


plot(50:length(newdf$Date),RMSE_T[50:length(newdf$Date)],type='o')


# Accuracy socres for true positives and true negatives
newdf$closure <- ifelse(newdf$OA >= 160,1,0)
confy <- newdf[as.numeric(n+1):nrow(newdf),]
confy$pred <- ifelse(gamPredictions$fit >=160,1,0)
confy$closure <- factor(confy$closure,levels=c(0,1))
confy$pred <- factor(confy$pred,levels=c(0,1))
confusionMatrix(confy$closure,confy$pred)


# Avergaed accuracy socres for classification of threshold breaches
dfs <- list()
dfs_m <- list()
fcst <- list()
for (i in 50:length(newdf$Date)){
  dfs[[i]] <- newdf[1:i,]
  dfs_m[[i]]<- dfs_m[[i]]<- gam(OA ~ s(Time,k=40), data=dfs[[i]],family = Gamma(link='log'))
  future <- as.data.frame(newdf[c(i+1,i+2,i+3,i+4,i+5),])
  names(future) <- names(newdf)
  fcst[[i]] <- predict(dfs_m[[i]], newdata=future, type="response", se.fit = T)
}


confy <- c()
for (i in 50:length(newdf$Date)){
  confy[[i]] <- as.data.frame(newdf$closure[c(i+1,i+2,i+3,i+4,i+5)])
  confy[[i]]$pred <- ifelse(fcst[[i]]$fit >=160,1,0)
  names(confy[[i]]) <-c('y','pred')
  confy[[i]]$y <- factor(confy[[i]]$y,levels=c(0,1))
  confy[[i]]$pred <- factor(confy[[i]]$pred,levels=c(0,1))
}


acc <-c()
for(i in 50:length(newdf$Date)){
  acc[[i]] <- confusionMatrix(confy[[i]]$y,confy[[i]]$pred)$overall[1]
}
acc <- as.data.frame(acc)
acc <- acc[50:as.numeric(nrow(newdf)-1),]
mean(acc)
