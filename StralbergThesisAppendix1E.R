# R code to develop boosted regression tree models and generate current predictions
library(raster) #Reading, writing, manipulating, analyzing and modeling of gridded spatial data
library(dismo) #Species distribution modeling
library(gbm) #Generalized boosted regression models
library(sampling) #Functions for drawing and calibrating samples

###################
#Data Preparation

#Bird data compiled by Boreal Avian Modelling Project, availability subject to data agreements with individual providers (http://www.borealbirds.ca/user/contact.php) 
#speclist: species list (by four-letter code)
#PC: point-count dataframe
#xx: survey attributes dataframe
#RES: list of dataframes with offsets for each species (see QPAD estimation in 'detect' package, https://github.com/psolymos/QPAD) 
#XY: Point-count coordinate dataframe
 
#Field definitions
#ABUND = raw count for a given species and survey
#SPEC = 4-letter species code
#PKEY = unique survey ID
#SS = unique point location ID
#PCODE = unique project code
#SITE = unique site code (collection of points locations)
#YEAR = survey year
#A = estimated area surveyed (QPAD offset component)
#p = estimated singing rate (QPAD offset component)
#q = estimated probability of detection (QPAD offset component)

surveydate <- aggregate(PC$ABUND, by=list("PKEY"=PC$PKEY,"YEAR"=PC$YEAR,"SS"=PC$SS,"PCODE"=PC$PCODE,"SITE"=PC$SITE), FUN=sum)

#Load current climate raster layers (requires raster package)
setwd(curclimate) #set current climate directory
clim <- list.files(curclimate, pattern =".asc$")
curclim<-stack(raster(clim[1]), raster(clim[2]))
i<-3
while (i <= length(clim)) {
	curclim <- addLayer(curclim,raster(clim[i]))
	i<-i+1
	}
	
#Load landcover raster layers (requires raster package)
setwd(landcover) #set landcover directory
curlc <- list.files(landcover, pattern =".asc$")
lcstack <- stack(raster(curlc[1]), raster(curlc[2]))
i<-3
while (i <= length(curlc)) {
	lcstack <- addLayer(lcstack,raster(curlc[i]))
	i<-i+1
	}

#Load topoedaphic raster layers (requires raster package)
setwd(topo) #set topoedaphic directory
topoedaphic <- list.files(topo, pattern =".asc$")
topostack <- stack(raster(topoedaphic[1]), raster(topoedaphic[2]))
i<-3
while (i <= length(topoedaphic)) {
	topostack <- addLayer(topostack,raster(topoedaphic[i]))
	i<-i+1
	}	

#Combine climate, landcover, and topoedaphic raster layers in a single stack (requires raster package)
climstack <- curclim
climstack <- addLayer(climstack, topostack,lcstack)	

#Extract climate, landcover, and topography data by XY coordinates
sites <- aggregate(PC$ABUND, by = list("SITE" = PC$SITE, "PCODE"=PC$PCODE, "SS"=PC$SS), FUN = sum)
XY <- merge(XY, sites[,1:3], by="SS")
climxy <- cbind(XY,extract(climstack,as.matrix(cbind(XY[,2],XY[,3]))))
climxy<-cbind(climxy,extract(nalc,as.matrix(cbind(climxy[,2],climxy[,3])))) #Extract NALCMS landcover class for point filtering
names(climxy)[ncol(climxy)] <- "LCC"
climxy <- na.omit(climxy)
climxy$ID <- as.factor(climxy$ID)
###################

###################
#Bootstrap sampling of data locations, repeated 11 times to create 11 different resampled datasets
#Requires raster and sampling packages
#Additional field definitions
#LCC = landcover class code from NALCMS landcover dataset
#YearGuess = earliest possible year of mapped disturbance

#Define sampling strata as intersection between SITE and PCODE
survey <- xx[1:2]
survey$ID <- row.names(survey)
dat1 <- merge(survey[,1:2],climxy,by="SS")
dat1$group <- paste(dat1$PCODE, dat1$SITE, dat1$ID, sep="-") 

#Remove urban, agricultural, and open water points
dat1 <- dat1[(dat1$LCC %in% c(15,16,17)) == FALSE,] 

#Remove points with surveys conducted after disturbance event
disturb <- read.csv("Disturb.csv") #Intersection of X-Y coordinates with anthropogenic disturbance from Global Forest Watch data 
dat2 <- merge(dat1, disturb[,c(1,65)])
dat3 <- merge(dat2, surveydate[,1:5], by=c("SS","PKEY","PCODE","SITE"))
dat3$YearGuess <- ifelse(dat3$YearGuess == 0,9999,dat3$YearGuess)
dat3$keep <- ifelse(dat3$YearGuess < dat3$YEAR, 0, 1) 
dat4 <- dat3[dat3$keep == 1,]

#Count number of surveys within group (SITE x PCODE)
dat4$count <- 1
count <- aggregate(dat4$count, by=list("group" = dat4$group), FUN = sum) 
names(count)[2] <- "count"
dat4 <- merge(dat4[,1:ncol(dat4)-1],count)

set.seed(72189) #Set seed for repeatability (different in each iteration)

#Sample one point from each group with more than 10 surveys
datmany <- dat4[dat4$count>10,] 
datsamp <- stratified(datmany, 1, 1)

#Sample one point from a third of the groups with fewer than 10 surveys
datfew <- dat4[dat4$count<11,]
datfew1 <- aggregate(datfew$count, by=list("group"=datfew$group), FUN = sum)
datsamp1 <- as.data.frame(sample(datfew1$group, size=nrow(datfew1)/3, replace=FALSE))
names(datsamp1)[1] <- "group"
datsamp2 <- merge(datfew,datsamp1)
datsamp3 <- stratified(datsamp2, 1, 1)

#Assign weights to sampled points based on inverse of total number of surveys #within a 20 km x 20 km (5 pixel by 5 pixel) area
datsamp4 <- rbind(datsamp[,c(1:62,64:66,69)],datsamp3[,c(1:62,64:66,69)])
r2 <- raster(clim[1])
samprast <- rasterize(datsamp4[,6:7], r2, field=1)
sampsum25 <- focal(samprast, w=5, na.rm=TRUE)
datsamp5 <- cbind(datsamp4,extract(sampsum25,as.matrix(cbind(datsamp4[,6],datsamp4[,7]))))
names(datsamp5)[ncol(datsamp5)] <- "sampsum25"
datsamp5$wt <- 1/datsamp5$sampsum25

datsampx <- datsamp5 #where is x is the iteration
####################

####################
#Build and save models, predict and export rasters
#Requires dismo, gbm, and raster packages
setwd(w) #set working directory

for (j in 1:length(speclist)) {
	specdat <- PC[PC$SPECIES == as.character(speclist[j]),]
	dat1 <- merge(datsamp5,specdat[,1:6],by=c("SS","PKEY","SITE","PCODE"),all.x=TRUE)
	dat1$SPECIES <- as.character(speclist[j])
	dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND))
	off <- as.data.frame(cbind(xx[1],RES[spp==speclist[j]]))
	off$Species <- speclist[j]
	names(off) <- c("PKEY","A","p","q","SPECIES")
	off$offset <- off$A * off$p * off$q
	dat2 <- merge(dat1,off[,c(1,5:6)])
	dat2$logoffset <- log(dat2$offset)

#Build and predict climate-only models
bird.brt.step1 <- gbm.step(dat2, gbm.y = (ncol(dat2)-2), gbm.x = c(11,12,14,16,19,25,37), family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5, offset=dat2$logoffset, site.weights=dat2$wt)
save(bird.brt.step1, file=paste(w,speclist[j],"_brt_clim1.RData",sep=""))
pdf(paste(w,speclist[j],"_brtplotclim1.pdf",sep=""))
gbm.plot(bird.brt.step1)
	gbm.plot.fits(bird.brt.step1, v=1:7)
dev.off()
rast <- predict(climstack, bird.brt.step1, type="response", n.trees=bird.brt.step1$n.trees)
writeRaster(rast, filename=paste(w,speclist[j],"_brtpredclim_1.asc",sep=""), format="ascii",overwrite=TRUE)
	
#Build and predict climate + land-use + topography models
bird.brt.step2 <- gbm.step(dat2, gbm.y = (ncol(dat2)-2), gbm.x = c(11,12,14,16,19,25,37,41,53,54,56,57), family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5, offset=dat2$logoffset, site.weights=dat2$wt)
save(bird.brt.step2, file=paste(w,speclist[j],"_brt_climtop1.RData",sep=""))
pdf(paste(w,speclist[j],"_brtplotclimtop1.pdf",sep=""))
	gbm.plot(bird.brt.step2)
	gbm.plot.fits(bird.brt.step2, v=1:12)
dev.off()
rast <- predict(climstack, bird.brt.step2, type="response", n.trees=bird.brt.step2$n.trees)
writeRaster(rast, filename=paste(w,speclist[j],"_brtpredclimtop_1.asc",sep=""), format="ascii",overwrite=TRUE)
	}
#####################


#####################
#Generate future projections

#Load landcover raster layers (requires raster package)
setwd(landcover) #set landcover directory
curlc <- list.files(landcover, pattern =".asc$")
lcstack <- stack(raster(curlc[1]), raster(curlc[2]))
i<-3
while (i <= length(curlc)) {
	lcstack <- addLayer(lcstack,raster(curlc[i]))
	i<-i+1
	}

#Load topoedaphic raster layers (requires raster package)
setwd(topo) #set topoedaphic directory
topoedaphic <- list.files(topo, pattern =".asc$")
topostack <- stack(raster(topoedaphic[1]), raster(topoedaphic[2]))
i<-3
while (i <= length(topoedaphic)) {
	topostack <- addLayer(topostack,raster(topoedaphic[i]))
	i<-i+1
	}	
	
setwd(ft) #set future climate directory (time period / senario)
#Future climate rasters organized in sub-directories named according to GCM
#Very time consuming; only 4/20 GCMs + ensembles predicted
gcms <- read.csv("gcms.csv") #List of GCMs for A2 scenario
gcms <- as.factor(as.character(gcms[1:20,1]))

for (i in 1:length(gcms)) {
	fclimtop <- list.files(paste(ft,gcms[i],sep=""), pattern=".asc$")
	setwd(paste(ft,gcms[i],sep="")) #set working directory to GCM
	futclimtop <- stack(raster(fclimtop[1]), raster(fclimtop[2]))
	k<-3
	while (k <= length(fclimtop)) {
		futclimtop <- addLayer(futclimtop,raster(fclimtop[k]))
		k<-k+1
		}
	climstack <- futclimtop
	climstack <- addLayer(climstack, topostack,lcstack
	
	setwd(w) #Set working directory
	models <- list.files(w, pattern=".RData")
	try(rm(bird.brt.step1))
	try(rm(bird.brt.step2))
	for (j in 1:length(models)) {
		base <- gsub(".RData","",models[j])
		if(file.exists(paste(w,"gcms/",base,"_",gcms[i],".asc",sep="")) == FALSE) {
		load(paste(w,models[j],sep=""))
		try(rast <- predict(climstack, bird.brt.step2, type="response", n.trees=bird.brt.step2$n.trees))
		try(rast <- predict(climstack, bird.brt.step1, type="response", n.trees=bird.brt.step1$n.trees))
		writeRaster(rast, filename=paste(w,"gcms/",base,"_",gcms[i],".asc",sep=""), format="ascii",overwrite=TRUE)
		png(paste(w,"gcms/",base,"_",gcms[i],".png",sep=""))
		plot(rast, zlim=c(0,1))
		dev.off()
		try(rm(bird.brt.step1))
		try(rm(bird.brt.step2))
		}
		}
}
#####################
