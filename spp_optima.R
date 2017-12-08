
####This script run a Weighted Averaging analysis to estimate species optima of occurrence along an environmetal gradient####
###It corresponds to step 1 in the enviromental mapping based on species records method describe in Zuquim et al. XXXX###


#load data
plotdata<-read.table("plotdata.txt", head=T)
head(plotdata)

##get package
install.packages("rioja")
library(rioja)

#separate species and environmetal data in two tables
species<-plotdata[,-c(1,2)]
head(species)
LogSoil<-plotdata$logCationContent
head(LogSoil)

#run the predictions
WA_ALL<- WA(species, LogSoil, tolDW=TRUE, mono=TRUE)
plot(WA_ALL)

#evaluate the predictions using cross-validation
crossWA<-crossval(WA_ALL)
summary(WA_ALL)

###check if the correlations are high and if they do, proceed to save optimas to be associated with ocuurance records for mapping

#assigns an object to the estimated species optima
otimosALL<-as.data.frame(coef(WA_ALL))
head(otimosALL)

#plot
plot(sort(otimosALL$Optima),)
