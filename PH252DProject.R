#######################################################
#########################ANALYSIS######################
#######################################################

# Clear workspace
rm(list=ls(all=TRUE))

#Check and set working directory to access files for data, CTMLE, and Wrapper code
#getwd()
setwd("/home/ph252d/R/PH252DFinal")

# Import wrapper functions for Super Learner library
source("SL_Wrappers.R")

# Import data from CSV file to dataframe data
data <- read.csv("faces1.csv")

# check number of rows in dataframe
# nrow(data)

# Only select columns from dataframe data needed for analysis
data<- subset(data, select=c(id, TOTFSL1, aP2524L1, gina1, gina2, wheeze, spring, summer, atopy))

# Omit all missing variables
data<- na.omit(data) 

# check number of observations with non-missing data
# nrow(data)

# make A (fungal spores) categorical, based on quantiles
# create quantile cutoff values
TOTFS1_25 <- quantile(data$TOTFSL1, c(0.25, 0.5, 0.75))[1]
TOTFS1_50 <- quantile(data$TOTFSL1, c(0.25, 0.5, 0.75))[2]
TOTFS1_75 <- quantile(data$TOTFSL1, c(0.25, 0.5, 0.75))[3]

# Create treatment variable A with levels 1-4, from low dose to high dose
data$A<-cut(data$TOTFSL1, c(0,TOTFS1_25,TOTFS1_50,TOTFS1_75, max(data$TOTFSL1) ), labels=c(1:4) )

# Make new column numeric
data$A <- factor(data$A, ordered=TRUE)

# Check A levels. We have 1371 observations with A=1 and 1361 with A=4...
# table(data$A)

# check dataframe data
# tail(data)

# table(data$wheeze)

#rename some variables to be more explicit about what goes into Super Learner
data$W_gina1  <- data$gina1
data$W_gina2  <- data$gina2

# data$W_PM25   <- as.numeric(data$aP2524L1 > 12)
# Let's make PM2.5 continuous
data$W_PM25   <- data$aP2524L1

data$W_spring <- data$spring
data$W_summer <- data$summer
data$W_atopy  <- data$atopy
data$Y        <- data$wheeze

#Making a new dataframe that will be used for SL - I left PM2.5 as continuous
FullData <- subset(data, select=c(id, W_PM25, W_gina1, W_gina2, W_spring, W_summer, W_atopy, A, Y))

# Start by loading SuperLearner for estimating Simple Substition Estimator 
library("SuperLearner")

# Library of candidate algorithms.

SL.library <- c("SL.glm", 
                "SL.step",
                "SL.glm.interaction",
                "SL.knn",
                #                "SL.gam",
                #                "SL.randomForest",
                #                "SL.glm.interaction_6",
                #                "SL.rpart",
                #                "SL.ipredbagg",
                #                "SL.polymars",
                "SL.nnet"
                #                "SL.bart"
                #                "SL.glmnet"
                #                "SL.bayesglm"
)

#########################################################################################################################################
## Simple Substitution G-Computation

#Create X dataframe, which is everything but Y
X<-subset(FullData, select=-c(Y, id))

# create data frames with A=4 and A=1
X1 <- X0 <-X

#Here, the high dose is A=4 while the low dose is A=1
X1$A<- 4 # will get pred outcome for ALL children with A=4
X0$A<- 1 # will get pred outcome for ALL children with A=1
# create newdata by stacking
newdata<- rbind(X,X1,X0)

system.time(Qinit<- SuperLearner(Y=FullData$Y, X=X, newX=newdata, SL.library=SL.library, id = FullData$id, family="binomial"))

Qinit

n<-nrow(FullData)
# n # not sure if we specify the total number of observations or total number of unique ids

# Simple Substitution (G-Computation) Estimator
QbarAW <- Qinit$SL.predict[1:n]
Qbar1W <- Qinit$SL.predict[(n+1): (2*n)]
Qbar0W <- Qinit$SL.predict[(2*n+1): (3*n)]

PsiHat.SS<-mean(Qbar1W - Qbar0W)

PsiHat.SS

#########################################################################################################################################
# IPTW Estimator

# estimate the treatment mechanism g_0(A|W) with multinomial logistic regression
library("nnet")
gAW.reg<-multinom(A~W_PM25*W_gina1*W_gina2*W_spring*W_summer*W_atopy, data=FullData)

#show that with ordinal, multgee, related packages, not able to integrate repeated measures... Maybe this
#accounts for why we got estimates slightly off from TMLE package...

# 2. predict each person's probability of getting a set level of exposure P(A=a | W's)
# this will return the predicted probability for each exposure
gAW.pred<- predict(gAW.reg, type="probs")
#head(gAW.pred)

# set up a vector for the predicted probabilities
gAW <- rep(NA,n)

# assign the appropriate predicted probabilities
gAW[FullData$A==1] <- gAW.pred[FullData$A==1, "1"]
gAW[FullData$A==2] <- gAW.pred[FullData$A==2, "2"]
gAW[FullData$A==3] <- gAW.pred[FullData$A==3, "3"]
gAW[FullData$A==4] <- gAW.pred[FullData$A==4, "4"]

summary(gAW)

# here we are the predicted probabilities of observations with A=1 OR A=4
#summary(gAW[FullData$A==1 | FullData$A==4])

# calc the weight as inverse of prob
wt = 1/gAW

# look at the distribution with summary and hist
summary(wt)

# Implement IPTW
# First the weighted average of outcomes, where kids with A=4 are
# weighted 1/gAW and kids with A !=4 have weight=0
meanY.4<- mean(wt*as.numeric(FullData$A==4)*FullData$Y)
#meanY.4

# then the weighted average of outcomes, where kids with A=1 are weighted
# 1/gAW and kids with A!=1 have weight=0
meanY.1<- mean(wt*as.numeric(FullData$A==1)*FullData$Y)
#meanY.1

PsiHat.IPTW <- meanY.4 - meanY.1
#PsiHat.IPTW

# Stabilized IPTW
PsiHat.IPTW.stab <- mean(wt*as.numeric(FullData$A==4)*FullData$Y)/mean(wt*as.numeric(FullData$A==4)) - mean(wt*as.numeric(FullData$A==1)*FullData$Y)/mean(wt*as.numeric(FullData$A==1))
PsiHat.IPTW.stab

#########################################################################################################################################
##TMLE Package

# Here, I assume that FullData is ALL the data from the CSV that I care about: W, A, Y, id

# Get the observations with A=1 or A=4
ObsData<- FullData[FullData$A==1 | FullData$A==4, ]

# Check to make sure we only have values for A=1 and A=4
#summary(ObsData)
#table(ObsData$A)

# number of children
n<- nrow(ObsData)
#n

# Recode the exposure so that A.binary=1 corresponds to A=4
# and A.binary=0 corresponds to A=1 
A.binary <- rep(0, n)
A.binary[ObsData$A==4] <- 1
# Check
#table(ObsData$A)
#table(A.binary)

# Load the tmle package
library(tmle)
tmle.QGiven.gGiven <- tmle(Y=ObsData$Y, A=A.binary, W=subset(ObsData, select=-c(A,Y,id)), family="binomial", Q.SL.library = SL.library, g.SL.library = SL.library, id = ObsData$id, verbose=TRUE)
#summary(tmle.QGiven.gGiven)
PsiHat.TMLE.package <- tmle.QGiven.gGiven$estimates$ATE$psi

#########################################################################################################################################

# ATE for all Estimators
Psi.Hat.estimates<- c(PsiHat.SS, PsiHat.IPTW, PsiHat.IPTW.stab, PsiHat.TMLE.package)
names(Psi.Hat.estimates)<-c("SimpSubs", "IPTW", "IPTW Stabilized", "TMLE Package")

# Inference Estimates
Psi.Hat.Inference <- c(PsiHat.TMLE.package, tmle.QGiven.gGiven$estimates$ATE$var.psi, tmle.QGiven.gGiven$estimates$ATE$CI[1], tmle.QGiven.gGiven$estimates$ATE$CI[2], tmle.QGiven.gGiven$estimates$ATE$pvalue)
names(Psi.Hat.Inference) <- c("ATE estimate" ,"Variance", "2.5%", "97.5%", "P-Value")

# 12/14/2013 @ 12:05pm

#ATE estimate      Variance          2.5%         97.5%       P-Value 
#0.0100669267  0.0009779673 -0.0512271128  0.0713609661  0.7475207287 

#> Psi.Hat.estimates
#SimpSubs            IPTW IPTW Stabilized    TMLE Package 
#0.00000000      0.02803393      0.01297198      0.01006693 

#########################################################################################################################################
# BOOTSTRAPPING!!!
# set seed
set.seed(252)

n <- nrow(FullData)

# number of bootstrap samples
B=500 # later change to B=500
estimates<- data.frame(matrix(NA, nrow=B, ncol=4)) # we have 4 dif't estimators right now

for(b in 1:B){
  
  # sample the indices 1 to n with replacement
  bootIndices<- sample(1:n, replace=T)
  bootData<- FullData[bootIndices,]
  
  
  ## Simple Substitution G-Computation for bootstrapping ##########################################
  
  #Create X dataframe, which is everything but Y
  X<-subset(bootData, select=-c(Y, id))
  
  # create data frames with A=4 and A=1
  X1 <- X0 <-X
  
  #Here, the high dose is A=4 while the low dose is A=1
  X1$A<- 4 # will get pred outcome for ALL children with A=4
  X0$A<- 1 # will get pred outcome for ALL children with A=1
  # create newdata by stacking
  newdata<- rbind(X,X1,X0)
  
  Qinit<- SuperLearner(Y=bootData$Y, X=X, newX=newdata, SL.library=SL.library, id = bootData$id, family="binomial")
  
  n<-nrow(FullData)
  
  # Simple Substitution (G-Computation) Estimator
  QbarAW <- Qinit$SL.predict[1:n]
  Qbar1W <- Qinit$SL.predict[(n+1): (2*n)]
  Qbar0W <- Qinit$SL.predict[(2*n+1): (3*n)]
  
  PsiHat.SS.b<-mean(Qbar1W - Qbar0W)
  
  ## IPTW for bootstrapping #######################################################################
  
  # estimate the treatment mechanism g_0(A|W) with multinomial logistic regression
  #library("nnet")
  gAW.reg<-multinom(A~W_PM25*W_gina1*W_gina2*W_spring*W_summer*W_atopy, data=bootData)
  
  # 2. predict each person's probability of getting a set level of exposure P(A=a | W's)
  # this will return the predicted probability for each exposure
  gAW.pred<- predict(gAW.reg, type="probs")
  #head(gAW.pred)
  
  # set up a vector for the predicted probabilities
  gAW <- rep(NA,n)
  
  # assign the appropriate predicted probabilities
  gAW[bootData$A==1] <- gAW.pred[bootData$A==1, "1"]
  gAW[bootData$A==2] <- gAW.pred[bootData$A==2, "2"]
  gAW[bootData$A==3] <- gAW.pred[bootData$A==3, "3"]
  gAW[bootData$A==4] <- gAW.pred[bootData$A==4, "4"]
  
  #summary(gAW)
  
  # here we are the predicted probabilities of observations with A=1 OR A=4
  #summary(gAW[bootData$A==1 | bootData$A==4])
  
  # calc the weight as inverse of prob
  wt = 1/gAW
  
  # look at the distribution with summary and hist
  #summary(wt)
  
  # Implement IPTW
  # First the weighted average of outcomes, where kids with A=4 are
  # weighted 1/gAW and kids with A !=4 have weight=0
  meanY.4<- mean(wt*as.numeric(bootData$A==4)*bootData$Y)
  #meanY.4
  
  # then the weighted average of outcomes, where kids with A=1 are weighted
  # 1/gAW and kids with A!=1 have weight=0
  meanY.1<- mean(wt*as.numeric(bootData$A==1)*bootData$Y)
  #meanY.1
  
  PsiHat.IPTW.b <- meanY.4 - meanY.1
  
  # Stabilized IPTW
  PsiHat.IPTW.stab.b <- mean(wt*as.numeric(bootData$A==4)*bootData$Y)/mean(wt*as.numeric(bootData$A==4)) - mean(wt*as.numeric(bootData$A==1)*bootData$Y)/mean(wt*as.numeric(bootData$A==1))
  
  ## TMLE for bootstrapping #######################################################################  
  
  # Get the observations with A=1 or A=4
  ObsData<- bootData[bootData$A==1 | bootData$A==4, ]
  
  # Check to make sure we only have values for A=1 and A=4
  # summary(ObsData)
  # table(ObsData$A)
  
  # number of children
  n_ObsData<- nrow(ObsData)
  #n
  
  # Recode the exposure so that A.binary=1 corresponds to A=4
  # and A.binary=0 corresponds to A=1 
  A.binary <- rep(0, n_ObsData)
  A.binary[ObsData$A==4] <- 1
  # Check
  #table(ObsData$A)
  #table(A.binary)
  
  # Load the tmle package
  #library(tmle)
  tmle.QGiven.gGiven <- tmle(Y=ObsData$Y, A=A.binary, W=subset(ObsData, select=-c(A,Y,id)), family="binomial", Q.SL.library = SL.library, g.SL.library = SL.library, id = ObsData$id, verbose=TRUE)
  #tmle.QGiven.gGiven <- tmle(Y=ObsData$Y, A=A.binary, W=subset(ObsData, select=-c(A,Y,id)), family="binomial", Q.SL.library = SL.library, g.SL.library = SL.library)
  PsiHat.TMLE.package.b <- tmle.QGiven.gGiven$estimates$ATE$psi
  
  estimates[b,]<- c(PsiHat.SS.b, PsiHat.IPTW.b, PsiHat.IPTW.stab.b, PsiHat.TMLE.package.b)
  print(b)
}

colnames(estimates)<-c("SimpSubs", "IPTW", "IPTW Stabilized", "TMLE Package")

write.csv(estimates, file="estimates.csv")

summary(estimates)

colMeans(estimates)

saving the histograms as a pdf
pdf(file="Final_Project_hist_boot.pdf")
par(mfrow=c(4,1))
hist(estimates[,1], main="Histogram of point estimates from the Simple Substitution estimator over 500 bootstrapped samples", xlab="Point Estimates")
hist(estimates[,2], main="Histogram of point estimates from IPTW estimator over 500 bootstrapped samples", xlab="Point Estimates")
hist(estimates[,3], main="Histogram of point estimates from Stabilized IPTW estimator over 500 bootstrapped samples", xlab="Point Estimates")
hist(estimates[,4], main="Histogram of point estimates from TMLE (Package) over 500 bootstrapped samples", xlab="Point Estimates")
dev.off()

diag(var(estimates))

#Variance for each estimator
varHat.SS<- var(estimates[,"SimpSubs"])
varHat.SS

varHat.IPTW<- var(estimates[,"IPTW"])
varHat.IPTW

varHat.IPTW.Stab<- var(estimates[,"IPTW Stabilized"])
varHat.IPTW.Stab

varHat.TMLE.Package<- var(estimates[,"TMLE Package"])
varHat.TMLE.Package

#Confidence Intervals for each estimator
c(PsiHat.SS-1.96*sqrt(varHat.SS), PsiHat.SS +1.96*sqrt(varHat.SS))
quantile(estimates[,"SimpSubs"], prob=c(0.025,0.975))

c(PsiHat.IPTW -1.96*sqrt(varHat.IPTW), PsiHat.IPTW +1.96*sqrt(varHat.IPTW))
quantile(estimates[,"IPTW"], prob=c(0.025,0.975))

c(PsiHat.IPTW.stab -1.96*sqrt(varHat.IPTW.Stab), PsiHat.IPTW.stab +1.96*sqrt(varHat.IPTW.Stab))
quantile(estimates[,"IPTW Stabilized"], prob=c(0.025,0.975))

c(PsiHat.TMLE.package -1.96*sqrt(varHat.TMLE.Package), PsiHat.TMLE.package +1.96*sqrt(varHat.TMLE.Package))
quantile(estimates[,"TMLE Package"], prob=c(0.025,0.975))
