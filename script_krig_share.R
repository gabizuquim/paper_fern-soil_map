
#load packages
library(sp)
library(raster)
library(gstat)
library(rgdal)

#Set a grid for the Amazon limits

#get the Eva&Huber shapefile
setwd("/myworkingspace")
AMZ<-shape <- readOGR(dsn = ".", layer = "ama_final2t_dissolve")
BB<-bbox(AMZ)
BBdf<-t(BB)
BBdf<-as.data.frame(BBdf)

####Create and plot the grid####
grd <- expand.grid(x = seq(from = min(BBdf$x), to = max(BBdf$x), by = 0.1),
                   y = seq(from = min(BBdf$y), to = max(BBdf$y), by = 0.1))
coordinates(grd) <- ~x + y
proj4string(grd) <- CRS("+proj=longlat +datum=WGS84")
gridded(grd) <- TRUE



####Averaging####
####soil+ferns average####
#setwd("/Volumes/ZUQUIMUTU/workspaceHD_Gabi_Zuquim20DEC2016/Papers em andamento copy/Ancillary data paper/ANALYSIS_NEW/codes shared")
data<-read.csv2("krig_dataSHARE.csv")
inTable3 <- na.omit(data)
head(inTable3)
#check for duplicates
dupl3 <- duplicated(inTable3)
sum(dupl3)

#seems there are none
inTable3 <- inTable3[!dupl3,]
inLon3 <- inTable3$x
inLat3 <- inTable3$y
inSB3 <- inTable3$LogSoilmean

lonMin3 <- min(inLon3)
lonMax3 <- max(inLon3)
latMin3 <- min(inLat3)
latMax3 <- max(inLat3)


res <- 1/60 #1arcmin resolution

r3 <- raster(xmn=lonMin3, ymn=lonMax3, xmx=latMin3, ymx=latMax3, res=res)

cellNrs3 <- cellFromXY(r3, cbind(inLon3, inLat3))
tab3 <- table(cellNrs3)
head(tab3)
r3[as.numeric(names(tab3))] <- tab3

d3 <- data.frame(coordinates(r3), count=r3[])
dOut3 <- na.omit(d3)
head(d3)
Log_soil.meanFS <- numeric(nrow(dOut3))+NA

for(i in 1:nrow(dOut3)){
  lon3 <- dOut3[i,1]
  lat3 <- dOut3[i,2]
  
  locs3 <- inLon3>=(lon3-res/2) & inLon3<=(lon3+res/2) & inLat3>=(lat3-res/2) & inLat3<=(lat3+res/2)
  Log_soil.meanFS[i] <- mean(inSB3[locs3])
}


dOut3$LogSoilmean <- Log_soil.meanFS


####End of Averaging####

#spatial data
soil_fernSpatial<-dOut3
coordinates(soil_fernSpatial) = ~x + y
proj4string(soil_fernSpatial) <- CRS("+proj=longlat +datum=WGS84")
plot(soil_fernSpatial)

##there is a cell with no values, some rounding error probably. We need to remove it
loc3 <- is.na(soil_fernSpatial$LogSoilmean)
#plot(dOut_spatial[loc,], add=TRUE, col="red", cex=2)
#session.zoom()
soil_fernSpatial <- soil_fernSpatial[!loc3,]


#variogram
vSoil_fern = variogram(LogSoilmean~1, soil_fernSpatial)
plot(vSoil_fern)
v.fitSoil_fern<-vgm(psill=0.19, model="Exp", nugget=0.3, range=1550)
plot(vSoil_fern,v.fitSoil_fern)

####correct####
vgm.fittedFS <- fit.variogram(vSoil_fern,v.fitSoil_fern)
plot(vSoil_fern,vgm.fittedFS)
####
KrigeSoil_fern = krige(formula=LogSoilmean~1, locations=soil_fernSpatial, newdata=grd, model = vgm.fittedFS)
#to raster and plot
KrigeSoil_fernR<-raster(KrigeSoil_fern)
#check this KrigeSoil_fernR_var<-raster(KrigeSoil_fern, layer=2)

plot(KrigeSoil_fernR)
plot(AMZ,add=T)
#save
writeRaster(KrigeSoil_fernR, "Krigins/KrigeSoil_fernR.tif", overwrite=T)

