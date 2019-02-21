#Open libraries for packages used
library(raster)
library(rgdal)
require(graphics)

#Import models
brt <- "L:/Boreal/DensityModelsV6/brt/cv/"
resultclip <- raster("L:/Boreal/StudyClip/resultclip.asc")

speclist <- read.csv("I:/BAM/BAMData/SpeciesClassesModv5.csv")
speclist <- as.factor(as.character(speclist[1:103,1]))

for (j in 1:length(speclist)) {
	try(dens <- raster(paste(brt,speclist[j],"_currmean.asc",sep="")))
	pres <- (1 - (exp(-1 * dens * pi)))
	writeRaster(pres, file=paste(brt,speclist[j],"_pres_currmean.asc",sep=""), format="ascii", overwrite=TRUE)
	png(paste(brt,speclist[j],"_pres_currmean.png",sep=""))
	plot(pres)
	dev.off()
	}
	
for (j in 1:length(speclist)) {
	try(dens <- raster(paste(brt,speclist[j],"_2080mean.asc",sep="")))
	pres <- (1 - (exp(-1 * dens * pi)))
	writeRaster(pres, file=paste(brt,speclist[j],"_pres_2080mean.asc",sep=""), format="ascii", overwrite=TRUE)
	png(paste(brt,speclist[j],"_pres_2080mean.png",sep=""))
	plot(pres)
	dev.off()
	}
	
for (j in 1:length(speclist)) {
	try(dens <- raster(paste(brt,speclist[j],"_2050mean.asc",sep="")))
	pres <- (1 - (exp(-1 * dens * pi)))
	writeRaster(pres, file=paste(brt,speclist[j],"_pres_2050mean.asc",sep=""), format="ascii", overwrite=TRUE)
	png(paste(brt,speclist[j],"_pres_2050mean.png",sep=""))
	plot(pres)
	dev.off()
	}
	
for (j in 1:length(speclist)) {
	try(dens <- raster(paste(brt,speclist[j],"_2020mean.asc",sep="")))
	pres <- (1 - (exp(-1 * dens * pi)))
	writeRaster(pres, file=paste(brt,speclist[j],"_pres_2020mean.asc",sep=""), format="ascii", overwrite=TRUE)
	png(paste(brt,speclist[j],"_pres_2020mean.png",sep=""))
	plot(pres)
	dev.off()
	}