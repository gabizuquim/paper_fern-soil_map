rm(list = ls())
####This script will perform steps 3 and 4 of the Methods section presented in Zuquim et al. (2019).
### Step 3: Obtain species-derived estimates and combine with direct environmental measurements
### Step 4: Generate maps by interpolating between data points run an ordinary kriging to 
#interpolate species optima and the values and produces an environmental map covering 
#the whole area of interest

#load packages
x<-c("sp", "raster","gstat", "rgdal")
lapply(x, library, character.only = TRUE)

#set working directory
setwd("C:/Data/R")

##############################################################################################################
##  Specify the extent and spatial resolution of the output product.                                        ##
##  Here, we use a shapefile to Amazonia to define the extent, and a resolution of 0.1 degree.              ##
##  Projection, extent and resolution can be set by the user by specifying values below                     ##
##############################################################################################################

##  Extent from shapefile Amazonia
amzExtent <- extent(readOGR(dsn = ".", layer = "myAreaOfInterestShapefile", verbose=FALSE))
out.xmin <- amzExtent@xmin
out.xmax <- amzExtent@xmax
out.ymin <- amzExtent@ymin
out.ymax <- amzExtent@ymax
out.xres <- 0.1
out.yres <- 0.1
out.proj <- CRS("+proj=longlat +datum=WGS84")

####  Create output grid, as a SpatialPixels object
grd <- expand.grid(x = seq(from = out.xmin, to = out.xmax, by = out.xres),
                   y = seq(from = out.ymin, to = out.ymax, by = out.yres))
coordinates(grd) <- ~x + y
proj4string(grd) <- out.proj
gridded(grd) <- TRUE

##################
##### Step 3: Obtain species-derived estimates and combine with direct environmental measurements
##  Preprocessing (reduction) of input data

########Species-based soil values only averaging#######
##  Fern-derived soil estimates are aggregated over an area of 1 arcsec squared ("out.xres" and "out.yres"), 
#and the average value is assigned to the centroid of that area

refdata.proj <- CRS("+proj=longlat +datum=WGS84")

#get species-based soil data
inTable<-read.csv2("fern_dataSHARE.csv")
inLon <- inTable$lon
inLat <- inTable$lat
inSB <- inTable$Optima

#create an output raster to populate
r <- raster(xmn=min(inLon), ymn=max(inLon), xmx=min(inLat), ymx=max(inLat), res=1/60, crs=refdata.proj)

#append the lat long data to the output raster
cellNrs <- cellFromXY(r, cbind(inLon, inLat))
uniqueCells <- as.numeric(levels(as.factor(cellNrs)))

#calculate and append the mean optmia per pixel
optimaMean <- sapply(uniqueCells, function(i){mean(inSB[cellNrs==i], na.rm=TRUE)})

#create a data frame with the result
fernAvg <- data.frame(x=xFromCell(r, uniqueCells),
                      y=yFromCell(r, uniqueCells),
                      logSoilMean=optimaMean,
                      source="fern")

#####end of fern-soil points averaging####

########Soil measuerements only averaging#######
## Soil measurements are aggregated over an area of 1 arcsec squared, 
#and the average value is assigned to the centroid of that area

#get direct soil measurements data 
inTable2 <- read.table("soil_dataSHARE.txt", head=T)
inLon2 <- inTable2$lon
inLat2 <- inTable2$lat
inSB2 <- inTable2$logSoil

#the following steps are the same as in lines 46 to 61 of this code, but using soils as input points
r2 <- raster(xmn=min(inLon2), ymn=max(inLon2), xmx=min(inLat2), ymx=max(inLat2), res=1/60, crs=refdata.proj)
cellNrs2 <- cellFromXY(r2, cbind(inLon2, inLat2))
uniqueCells2 <- as.numeric(levels(as.factor(cellNrs2)))

optimaMean2 <- sapply(uniqueCells2, function(i){mean(inSB2[cellNrs2==i], na.rm=TRUE)})

soilAvg <- data.frame(x=xFromCell(r2, uniqueCells2),
                      y=yFromCell(r2, uniqueCells2),
                      logSoilMean=optimaMean2,
                      source="soil")

####end of soil averaging####
####End of Averaging####

#Combine the averaged fern-derived and soil measurements pixel values
soil_fernSpatial <- rbind(fernAvg, soilAvg)
##transform to a spatial object
coordinates(soil_fernSpatial) = ~x + y
proj4string(soil_fernSpatial) <- refdata.proj
#crop the spatial files to the original grid extent
soil_fernSpatial<-crop(soil_fernSpatial,grd)

### Step 4: Generate maps by interpolating between data points
#################################################################################
####             Here starts the actual ordinary kriging                     ####
####      1. model fitting, 2. kriging, 3. inspect the resulting map         ####
#################################################################################

#1. Model fitting
#variogram
vSoil_fern = variogram(logSoilMean~1, soil_fernSpatial)
plot(vSoil_fern)
v.fitSoil_fern<-vgm(psill=0.19, model="Exp", nugget=0.3, range=1550)

####adjusts the variogram####
vgm.fittedFS <- fit.variogram(vSoil_fern,v.fitSoil_fern)

#2. Kriging
####Run the kriging using the variogram defined above
# This step may take several hours depending on the amount of input data and pixel size
KrigeSoil_fern = krige(formula=logSoilMean~1, locations=soil_fernSpatial, newdata=grd, model = vgm.fittedFS)

####the map is ready!####

#3. Inspect results
#for better visualization, transform the object to a raster
KrigeSoil_fernR<-raster(KrigeSoil_fern)
#check this KrigeSoil_fernR_var<-raster(KrigeSoil_fern, layer=2)

#you can also plot it.
plot(KrigeSoil_fernR)
#... and save. 
writeRaster(KrigeSoil_fernR, "KrigeSoil_mapR.tif")

#... and extract values of your location points using the extract function from raster package
