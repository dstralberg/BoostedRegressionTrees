library(ggplot2)
library(reshape)
library(raster)
brt <- "L:/Boreal/DensityModelsV6/brt/cv/"
speclist <- read.csv("I:/BAM/BAMData/SpeciesClassesModv5.csv")
speclistB <- speclist[speclist$Boreal10==1,c(1,4)]
speclistB <- as.factor(as.character(speclistB[1:80,1]))

a <- raster(paste(brt,"_curtotdens.asc",sep=""))
i <- raster(paste(brt,speclistB[1],"_currmean.asc",sep=""))
p <- i/a
lnp <- log(p)
H <- p*lnp
for (j in 2:length(speclistB)) {
	i <- raster(paste(brt,speclistB[j],"_currmean.asc",sep=""))
	p <- i/a
	lnp <- log(p)
	H <- p*lnp
	H <- H + (p*lnp)
	}
Dcurr <- exp(-1*H)	
writeRaster(Dcurr, file=paste(brt,"_DcurrB10.asc",sep=""), format="ascii", overwrite=TRUE)
png(paste(brt,"_DcurrB10.png",sep=""))
plot(Dcurr)
dev.off()

a <- raster(paste(brt,"_2020totdens.asc",sep=""))
i <- raster(paste(brt,speclistB[1],"_2020mean.asc",sep=""))
p <- i/a
lnp <- log(p)
H <- p*lnp
for (j in 2:length(speclistB)) {
	i <- raster(paste(brt,speclistB[j],"_2020mean.asc",sep=""))
	p <- i/a
	lnp <- log(p)
	H <- p*lnp
	H <- H + (p*lnp)
	}
D2020 <- exp(-1*H)	
writeRaster(D2020, file=paste(brt,"_D2020B10.asc",sep=""), format="ascii", overwrite=TRUE)
png(paste(brt,"_D2020B10.png",sep=""))
plot(D2020)
dev.off()
	
a <- raster(paste(brt,"_2050totdens.asc",sep=""))
i <- raster(paste(brt,speclistB[1],"_2050mean.asc",sep=""))
p <- i/a
lnp <- log(p)
H <- p*lnp
for (j in 2:length(speclistB)) {
	i <- raster(paste(brt,speclistB[j],"_2050mean.asc",sep=""))
	p <- i/a
	lnp <- log(p)
	H <- p*lnp
	H <- H + (p*lnp)
	}
D2050 <- exp(-1*H)	
writeRaster(D2050, file=paste(brt,"_D2050B10.asc",sep=""), format="ascii", overwrite=TRUE)
png(paste(brt,"_D2050B10.png",sep=""))
plot(D2050)
dev.off()

a <- raster(paste(brt,"_2080totdens.asc",sep=""))
i <- raster(paste(brt,speclistB[1],"_2080mean.asc",sep=""))
p <- i/a
lnp <- log(p)
H <- p*lnp
for (j in 2:length(speclistB)) {
	i <- raster(paste(brt,speclistB[j],"_2080mean.asc",sep=""))
	p <- i/a
	lnp <- log(p)
	H <- p*lnp
	H <- H + (p*lnp)
	}
D2080 <- exp(-1*H)	
writeRaster(D2080, file=paste(brt,"_D2080B10.asc",sep=""), format="ascii", overwrite=TRUE)
png(paste(brt,"_D2080B10.png",sep=""))
plot(D2080)
dev.off()

#species richness
	
p <- raster(paste(brt,speclistB[1],"_pres_currmean.asc",sep=""))
SRcurr <- p
for (j in 2:length(speclistB)) {
	p <- raster(paste(brt,speclistB[j],"_pres_currmean.asc",sep=""))
	SRcurr <- SRcurr + p
	}
writeRaster(SRcurr, file=paste(brt,"_SRcurrB10.asc",sep=""), format="ascii", overwrite=TRUE)
png(paste(brt,"_SRcurrB10.png",sep=""))
plot(SRcurr)
dev.off()

p <- raster(paste(brt,speclistB[1],"_pres_2020mean.asc",sep=""))
SR2020 <- p
for (j in 2:length(speclistB)) {
	p <- raster(paste(brt,speclistB[j],"_pres_2020mean.asc",sep=""))
	SR2020 <- SR2020 + p
	}
writeRaster(SR2020, file=paste(brt,"_SR2020B10.asc",sep=""), format="ascii", overwrite=TRUE)
png(paste(brt,"_SR2020B10.png",sep=""))
plot(SR2020)
dev.off()
	
p <- raster(paste(brt,speclistB[1],"_pres_2050mean.asc",sep=""))
SR2050 <- p
for (j in 2:length(speclistB)) {
	p <- raster(paste(brt,speclistB[j],"_pres_2050mean.asc",sep=""))
	SR2050 <- SR2050 + p
	}
writeRaster(SR2050, file=paste(brt,"_SR2050B10.asc",sep=""), format="ascii", overwrite=TRUE)
png(paste(brt,"_SR2050.png",sep=""))
plot(SR2050)
dev.off()

p <- raster(paste(brt,speclistB[1],"_pres_2080mean.asc",sep=""))
SR2080 <- p
for (j in 2:length(speclistB)) {
	p <- raster(paste(brt,speclistB[j],"_pres_2080mean.asc",sep=""))
	SR2080 <- SR2080 + p
	}
writeRaster(SR2080, file=paste(brt,"_SR2080B10.asc",sep=""), format="ascii", overwrite=TRUE)
png(paste(brt,"_SR2080.png",sep=""))
plot(SR2080)
dev.off()	