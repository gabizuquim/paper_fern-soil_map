rm(list = ls())
####This script run a Weighted Averaging analysis to estimate species optima of occurrence along 
#an environmental gradient####
###It corresponds to step 1 in the environmental mapping based on species records method described 
#in Zuquim et al. 2019####

#Load data
#Species and environmental data should be entered as a single species abundance table:
#the first column contains a unique ID for each inventory plot, 
#the second column values of the environmental variable of interest, 
#all subsequent values represent the abundance of individual species observed in each plot

plotdata<-read.table("plotdata.txt", head=T)
##Get package
if (!require("rioja", character.only = TRUE)) install.packages("rioja"); library(rioja)

#Create two separate tables: one for species and one for environmental data
species<-subset(plotdata, select=-c(plotcode,logCationContent))
LogSoil<-plotdata$logCationContent

#Predict soil cation concentration (or the environmental variable of interest)
WA_ALL<- WA(species, LogSoil, tolDW=TRUE, mono=TRUE)

#Plot predicted values to see how optima is distributed along the gradient
plot(WA_ALL)

#Evaluate the predictions using cross-validation (see Table 1 in Zuquim et al. 2019)
crossWA<-crossval(WA_ALL)

##Check if the correlations (R2) between predicted and 
#observed values are satisfactory (e.g. >0.7) are high. If so, it means that the species 
#pool chosen are in fact a good set of indicator species.
performance(WA_ALL)

##Provided the species are good indicators, the estimated species optima will be kept 
#(as data frame) to be used in the next steps towards soil mapping
otimosALL<-as.data.frame(coef(WA_ALL))

