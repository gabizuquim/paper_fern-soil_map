
#load packages
library(sp)
library(raster)
library(gstat)
library(rgdal)

#Set a grid for the Amazon limits

#set your area of interest shapefile
AMZ<-readOGR(dsn = ".", layer = "myAreaOfInterestShapefile")
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
####plant derived soil values average####
inTable<-read.csv2("fern_dataSHARE.csv")
names(inTable)
inLon <- inTable$lon
inLat <- inTable$lat
inSB <- inTable$Optima

lonMin <- min(inLon)
lonMax <- max(inLon)
latMin <- min(inLat)
latMax <- max(inLat)

res <- 1/60 #1arcmin resolution

library(raster)
r <- raster(xmn=lonMin, ymn=lonMax, xmx=latMin, ymx=latMax, res=res)

cellNrs <- cellFromXY(r, cbind(inLon, inLat))
tab <- table(cellNrs)
head(tab)
r[as.numeric(names(tab))] <- tab

d <- data.frame(coordinates(r), count=r[])
dOut <- na.omit(d)
head(d)

LogOptima.mean <- numeric(nrow(dOut))+NA

for(i in 1:nrow(dOut)){
  lon <- dOut[i,1]
  lat <- dOut[i,2]
  
  locs <- inLon>=(lon-res/2) & inLon<=(lon+res/2) & inLat>=(lat-res/2) & inLat<=(lat+res/2)
  LogOptima.mean[i] <- mean(inSB[locs])
}

dOut$LogOptimamean <- LogOptima.mean
#####end of fern-soil points averaging####

########Soils only averaging#######
soil_0_30Spatial<-read.csv2("soil_dataSHARE.csv")
########Soils only averaging#######
coordinates(soil_0_30Spatial) = ~longitude+latitude
proj4string(soil_0_30Spatial) <- CRS("+proj=longlat +datum=WGS84")
soil_0_30crop<- crop(soil_0_30Spatial, limitAMAZ)
plot(soil_0_30crop, col="red")
dim(soil_0_30crop)
soil_0_30crop_df<-as.data.frame(soil_0_30crop)
#check for duplicates
dupl2 <- duplicated(soil_0_30crop_df)
sum(dupl2)

#seems there are none
inTable2 <- soil_0_30crop_df[!dupl2,]
write.csv(inTable2, "/Volumes/ZUQUIMUTU/workspaceHD_Gabi_Zuquim20DEC2016/Papers em andamento copy/Ancillary data paper/ANALYSIS_NEW/SoilMapFern/krigins/soilpoints.csv")
dim(inTable2)

names(inTable2)
inLon2 <- inTable2$lon
inLat2 <- inTable2$lat
inSB2 <- inTable2$logSoil

lonMin2 <- min(inLon2)
lonMax2 <- max(inLon2)
latMin2 <- min(inLat2)
latMax2 <- max(inLat2)

#sP <- SpatialPoints(cbind(inLon, inLat))
#plot(sP, pch=".")

res <- 1/60 #1arcmin resolution

library(raster)
r2 <- raster(xmn=lonMin2, ymn=lonMax2, xmx=latMin2, ymx=latMax2, res=res)

cellNrs2 <- cellFromXY(r2, cbind(inLon2, inLat2))
tab2 <- table(cellNrs2)
head(tab2)
r2[as.numeric(names(tab2))] <- tab2
# plot(r)
# rVals <- values(r)
# sum(!is.na(rVals))

d2 <- data.frame(coordinates(r2), count=r2[])
dOut2 <- na.omit(d2)
head(d2)
#sP <- SpatialPoints(cbind(dOut[,1], dOut[,2]))
#plot(sP, pch=".")

Log_soil.mean <- numeric(nrow(dOut2))+NA

for(i in 1:nrow(dOut2)){
  lon2 <- dOut2[i,1]
  lat2 <- dOut2[i,2]
  
  locs2 <- inLon2>=(lon2-res/2) & inLon2<=(lon2+res/2) & inLat2>=(lat2-res/2) & inLat2<=(lat2+res/2)
  Log_soil.mean[i] <- mean(inSB2[locs2])
}


dOut2$LogSoilmean <- Log_soil.mean
####end of soil averaging####
####End of Averaging####

dOut$source<-"fern"
colnames(dOut)[4]<-"LogSoilmean"
dOut2$source<-"soil"

#spatial data
soil_fernSpatial<-rbind(dOut,dOut2)
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

