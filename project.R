library(rgeos)
library(rgdal)
library(RColorBrewer)
library(GISTools)
library(raster)
library(gstat) 
library(nlme)
library(sp)
require(geosphere)
library(spdep)
library(spgwr)
library(GWmodel)
library(moments)
library(plotrix)
library(spatstat)
library(maptools)
library(plyr)


##### PREP POLLUTION DATA #####

#read pollution data
ozone = read.csv("OZONE_PICKDATA_2016-4-30.csv", header=T, sep=",")
partmat = read.csv("PM25HR_PICKDATA_2016-4-30.csv", header=T, sep=",")

#read EPA monitoring shapefile data
monitor <- readOGR(dsn="airmonitoringstations.shp", layer="airmonitoringstations")

#extract monitoring station points within the South coast
SC.monitor = monitor[monitor$AIRBASIN %in% c("South Coast"),]

#reproject the data to a suitable projection
#utm used because of the scale of the analysis
SC.monitor.t = spTransform(SC.monitor, CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 + towgs84=0,0,0"))

#read the california air basin spatial dataset
Ca.AirBasin = readOGR(dsn="CaAirBasin.shp",layer="CaAirBasin")

#extract the south coast air basin from the spatial dataset
SC.AirBasin = Ca.AirBasin[Ca.AirBasin$NAME %in% c("South Coast"),]

#reproject the south coast air basin spatial dataset to match the projection of the monitoring station dataset
SC.AirBasin.t = spTransform(SC.AirBasin, CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))

#### PREP CENSUS DATA ####

#read california census tracts
all.tracts = shapefile("cb_2015_06_tract_500k.shp")

#reproject the data to match the projection of the pollution datasets
all.tracts.t = spTransform(all.tracts, CRS("+proj=longlat +datum=NAD83 +no_defs + ellps=GRS80 +towgs84=0,0,0"))

#extract census tracts from that intersect with the sc air basin
SC.tracts.t <- intersect(all.tracts.t, SC.AirBasin.t)

#read in the census data and merge with census tracts
census.data = read.csv("ACS_14_5YR_B02001_with_ann.csv") #race
census2.data = read.csv("ACS_14_5YR_B19013_with_ann.csv") #income

#race data set up
race.data = census.data[,c("GEO.id", "HD01_VD01", "HD01_VD02","HD01_VD03", "HD01_VD04","HD01_VD05","HD01_VD06")]
colnames(race.data)[1:7] = c("AFFGEOID","Total","White","Black","AIndian","Asian","Hawaiian")
#delete unneccessary info
race.data <- race.data[c(-1),]

#function for getting diversity measure
getdiversity <- function(rowval) {
	n = rowval
	tempmax = max(as.numeric(as.character(race.data$White[n])),as.numeric(as.character(race.data$Black[n])),as.numeric(as.character(race.data$AIndian[n])),as.numeric(as.character(race.data$Asian[n])),as.numeric(as.character(race.data$Hawaiian[n])))
	diversity = tempmax / as.numeric(as.character(race.data$Total[n]))
	return(diversity)
} 
#add new measure onto dataset
for (i in 1:nrow(race.data)) {
	race.data[i] = cbind(race.data[i], diversity=getdiversity(i))
}

race.data = read.csv("race-edit.csv") #edited race data

#income data set up
income.data = census2.data[,c("GEO.id", "HD01_VD01")]
colnames(income.data)[1:2] = c("AFFGEOID","Income")

#change the name of GEO.id so it can be merged with shapefile
colnames(census.data)[1] = "AFFGEOID"

#merge the sc census tracts with the census data
SC.tracts.t = sp::merge(SC.tracts.t, race.data, by='AFFGEOID')
SC.tracts.t = sp::merge(SC.tracts.t, income.data, by='AFFGEOID')

#remove missing values
SC.tracts.t <- SC.tracts.t[complete.cases(SC.tracts.t$Diversity),]

#choropleth map of diversity values
shades = auto.shading(as.numeric(SC.tracts.t$Diversity),cols=brewer.pal(9,"Blues"))
choropleth(SC.tracts.t,SC.tracts.t$Diversity,shading=shades, border=NA)
plot(SC.monitor.t, add=T, col="red", cex=0.5, pch=16)
plot(SC.AirBasin.t, add=T, border="lightgrey")

#map the south coast air basin with census tracts and monitoring station locations
par(mfrow=c(1,1))
plot(SC.tracts.t, border="lightgrey", col=SC.tracts.t$Diversity)
plot(SC.AirBasin.t, add=T)
plot(SC.monitor.t, add=T, col="red", cex=0.5, pch=16)

#### PROCESS OZONE DATA ####

#calculate mean and max ozone level for each site for all readings
mean.partmat <- aggregate(value ~ site, partmat, mean)
max.partmat <- aggregate(value ~ site, partmat, max)

#join the mean and max ozone values to their respective monitoring stations
#and rename the first col of monitoring data to site for merging purposes
names(SC.monitor.t)[1] <- "site"

#merge SC.monitor with ozone mean and max data
mrg.tab.mean <- sp::merge(SC.monitor.t, mean.partmat, by="site", all.x=FALSE)
mrg.tab.max <- sp::merge(SC.monitor.t, max.partmat, by="site", all.x=FALSE)

#create a max and a mean spatialPointDataFrame from data
partmat.mean.spdf = na.omit(mrg.tab.mean)
partmat.max.spdf = na.omit(mrg.tab.max)

#add the merged moniotring stations to map (there are less stations now as a result of merge)
plot(partmat.mean.spdf,add=T,col="blue")

#### CORRELATION ANALYSIS ####

head(SC.tracts.t)

DV = as.numeric(SC.tracts.t$White)
INDV = as.numeric(SC.tracts.t$Black)
plot(DV,INDV)
cor.test(DV,INDV)

#### REGRESSION ANALYSIS ####

#first, plot data and fit linear line through it based on least squares
lmfit = lm(DV~INDV)
abline(lmfit)

#conduct regression analysis
linear.model = lm(DV~INDV)
summary(linear.model)
resids = residuals(linear.model)

#prepare data for plotting residuals from the regression analysis
#first convert census tracts to a data frame object (need to not use spdf for now)
tracts.df = as.data.frame(SC.tracts.t)

#combine census data with residuals from census data
tracts.df = cbind(tracts.df[1],resids)

#change col header name
colnames(tracts.df)[2] = "residuals"

#merge the residuals data set with the census data set
SC.tracts.t = merge(SC.tracts.t, tracts.df, by='AFFGEOID')

#set min and max lat/long for mapping residuals
minlat = min(SC.monitor$LATITUDE - 0.1)
maxlat = max(SC.monitor$LATITUDE + 0.1)
minlong = min(SC.monitor$LONGITUDE - 0.1)
maxlong = max(SC.monitor$LONGITUDE + 0.1)

#create a colouring scheme for choropleth map
shades = auto.shading(as.numeric(SC.tracts.t$residuals),cols=brewer.pal(9,"Greens"))

#create a choropleth map
choropleth(SC.tracts.t, SC.tracts.t$residuals, border=NA, shading=shades, xlim=c(minlong,maxlong), ylim=c(minlat,maxlat))

#conduct a global moran's i to see if residuals exhibit spatial autocorrelation
#create a weighted matrix file in order to perform moran's i analysis
SC.tracts.nb = poly2nb((SC.tracts.t), queen=FALSE)
SC.tracts.nb
#convert neighbour object to a listw object
SC.tract.lw = nb2listw(SC.tracts.nb, zero.policy = TRUE)
moran.test(resids,SC.tract.lw, zero.policy = TRUE)

#perform geographically weighted regression (only do this when spatially autocorrelated)
#(order of the dependent ~ independent terms is intentionally the opposite from the linear model)
gwr.res = gwr.basic(DV~INDV, data=SC.tracts.t, bw=20000, kernel='gaussian')
#examine the model summary, the results of the gloval regression, and the results of the geographically weighted regression
#note the gwr coefficeient estimates and the r squared value
gwr.res
head(gwr.res$SDF)

#function that will create a map of points
quick.map <- function(spdf,var,legend.title,main.title) {
	x <- spdf@data[,var]
	cut.vals <- pretty(x)
	x.cut <- cut(x,cut.vals)
	cut.levels <- levels(x.cut)
	cut.band <- match(x.cut, cut.levels)
	colors <- brewer.pal(length(cut.levels), 'Reds')
	par(mar=c(1,1,1,1))
	plot(SC.AirBasin,col='grey85')
	title(main.title)
	plot(spdf,add=TRUE,col=colors[cut.band],pch=16)
	legend('topleft',cut.levels,col=colors,pch=16,bty='n',title=legend.title)
}

#using the geographically weighted regresion data, plot the gw regression coefficients
quick.map(gwr.res$SDF,"DV","legend title","map title")
#map of geographically weighted regression residuals
quick.map(gwr.res$SDF, "residual","legend title","map title")


#### VORONOI TESSELATION ####
##(in order to get control points for kriging)

#prepare data
voronoipolygons = function(layer) {
	require(deldir)
	crds = layer@coords
	z = deldir(crds[,1], crds[,2])
	w = tile.list(z)
	polys = vector(mode='list', length=length(w))
	require(sp)
	for (i in seq(along=polys)) {
		pcrds = cbind(w[[i]]$x, w[[i]]$y)
		pcrds = rbind(pcrds, pcrds[1,])
		polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
	}
	SP = SpatialPolygons(polys)
	voronoi = SpatialPolygonsDataFrame(SP, data=data.frame(dummy = seq(length(SP)), row.names=sapply(slot(SP, 'polygons'), function(x) slot(x, 'ID'))))
}

#calculate voronoi tesselation
partmat.voro = voronoipolygons(partmat.mean.spdf)
proj4string(partmat.voro) = proj4string((partmat.mean.spdf)) #define the projection

#plot voronoi polygons on top of census data as a visual check
plot(SC.tracts.t, main="Voronoi for PM2.5", border="lightgrey", col=NA)
plot(partmat.voro, , border="blue", col=NA, add=T)

#create points for calculating the interpolated surface
s.grid = spsample(partmat.voro, type="regular", n=6000)

#define the same projection
proj4string(s.grid) = proj4string(partmat.mean.spdf)

#### PERFORM KRIGING ####

#estimating semivariogram
evgm <- variogram(value ~ 1, partmat.mean.spdf, boundaries=seq(0,200,l=16))
plot(evgm)
nugget = 2
sill = 35
range = 120
plot(evgm,model=vgm(sill,"Mat",range,nugget),main="Semivariogram for PM2.5 in Los Angeles", xlab="Distance", ylab="Semivariance")
fvmg <- fit.variogram(evgm, vgm(1.4e-05, "Mat", 50, 0))
fvmg <- fit.variogram(evgm, vgm(27, "Mat", 50, 0))
plot(evgm, fvmg)

#perform kriging
krig.est <- krige(partmat.mean.spdf$value~1, partmat.mean.spdf, newdata=s.grid, model=fvmg)

#convert kriging surface to a raster dataset
r = raster(extent(krig.est), ncols=115, nrows=52, crs=proj4string(krig.est))
# dim(krig.est) #"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
krig.rst = rasterize(krig.est, r, "var1.pred")


#plot raster interpolated surface with census tracts
plot(krig.rst, col=heat.colors(30, alpha=1),main="Kriging Interpolated Surface for PM2.5 in Los Angeles")
plot(SC.AirBasin.t, add=T)
plot(SC.tracts.t, main="Kriging", border="darkgrey", col=NA, add=T, lwd=0.4)
plot(centroid.spdf, col="black", pch=16, cex=0.4, add=T) #do after centroids

#### JOINING POLLUTION AND CENSUS DATA ####

#calculate centroids of census tracts polygons
centroid <- centroid(SC.tracts.t)
centroid.spdf <- as.data.frame(centroid)
coordinates(centroid.spdf) <- ~ V1 + V2

#extract numeric values from kriging raster
partmat.krig.c = raster::extract(krig.rst, centroid.spdf)

#convert numeric values to data.frame
partmat.krig.c.as.df = as.data.frame(partmat.krig.c)

#add data frame to SpatialPoints 
#this allows us to convert SpatialPoints (centroid.spdf) to SPDF
partmat.c.spdf = SpatialPointsDataFrame(centroid.spdf, partmat.krig.c.as.df)

#merge data extracted from kriging raster with census data for income by centroid data with attached information of income
census.pt = raster::intersect(partmat.c.spdf, SC.tracts.t) #use ozone.spdf instead of centroid scpd to assure spatial overlap between both datasets
dim(census.pt)
plot(census.pt)

#keep only informaion about diversity and pollution
fin.tab = data.frame(census.pt@data$partmat.krig.c, census.pt@data$Diversity)

#change names of columns to readable format
names(fin.tab) = c("mean.partmat", "div")

#check names and table
head(fin.tab)

#now data is ready to perform correlation or regression on pollution and socioeconomic variables


##### DESCRIPTIVE STATISTICS ####

all.data = read.csv("SCtracts-t.csv", header=T, sep=",")

partmat = read.csv("PM25HR_PICKDATA_2016-4-30.csv", header=T, sep=",")
head(partmat)
mean.partmat <- aggregate(value ~ site, partmat, mean)
max.partmat <- aggregate(value ~ site, partmat, max)

#mean of mean particulate matter 2.5 across whole of LA
mean(mean.partmat$value)
#same of max
mean(max.partmat$value)

median(mean.partmat$value)
names(sort(-table(mean.partmat$value)))[1] #mode
sd(mean.partmat$value)
skewness(mean.partmat$value)
kurtosis(mean.partmat$value)
sd(mean.partmat$value)/mean(mean.partmat$value) #coefficient of variation
min(mean.partmat$value)
max(mean.partmat$value)

race.data = read.csv("race-edit.csv") #edited race data
head(race.data)
mean(race.data$Diversity)
median(race.data$Diversity)
names(sort(-table(race.data$Diversity)))[1] #mode
sd(race.data$Diversity)
skewness(race.data$Diversity)
kurtosis(race.data$Diversity)
sd(race.data$Diversity)/mean(race.data$Diversity) #coefficient of variation
min(race.data$Diversity)
max(race.data$Diversity)

##using all.data##

#pm2.5 values
mean(all.data$mean.partmat)
median(all.data$mean.partmat)
names(sort(-table(all.data$mean.partmat)))[1] #mode
sd(all.data$mean.partmat)
skewness(all.data$mean.partmat)
kurtosis(all.data$mean.partmat)
sd(all.data$mean.partmat)/mean(all.data$mean.partmat) #coefficient of variation
min(all.data$mean.partmat)
max(all.data$mean.partmat)

#diversity values
mean(all.data$Diversity)
median(all.data$Diversity)
names(sort(-table(all.data$Diversity)))[1] #mode
sd(all.data$Diversity)
skewness(all.data$Diversity)
kurtosis(all.data$Diversity)
sd(all.data$Diversity)/mean(all.data$Diversity) #coefficient of variation
min(all.data$Diversity)
max(all.data$Diversity)

#histograms
hist(all.data$mean.partmat, main="Distribution of PM2.5", xlab="Particulate Matter 2.5", ylab="Frequency", col="lightblue")
hist(all.data$Diversity, main="Distribution of Diversity", xlab="Diversity", ylab="Frequency", col="lightblue")

#choropleth map of diversity values
shades = auto.shading(as.numeric(race.data$Diversity),cols=brewer.pal(5,"Blues"))
choropleth(SC.tracts.t,all.data$Diversity,shading=shades, border=NA, main="Diversity Levels in Los Angeles")
choro.legend(px="bottomright", sh=shades, bg="white", title="Diversity")
plot(SC.monitor.t, add=T, col="red", cex=0.5, pch=16)
plot(SC.AirBasin.t, add=T, border="lightgrey")

#choropleth map of particulate matter values
shades2 = auto.shading(as.numeric(all.data$mean.partmat),cols=brewer.pal(9,"Blues"))
choropleth(SC.tracts.t,all.data$mean.partmat,shading=shades2,border=NA, main="PM2.5 Levels in Los Angeles")
choro.legend(px="bottomright",sh=shades2,bg="white", title="PM2.5")
plot(SC.monitor.t, add=T, col="red", cex=0.5, pch=16)
plot(SC.AirBasin.t, add=T, border="lightgrey")



#### SPATIAL DESCRIPTIVE & SPATIAL INFERENTIAL MEASURES ####

#read EPA monitoring shapefile data
monitor <- readOGR(dsn="airmonitoringstations.shp", layer="airmonitoringstations")

#extract monitoring station points within the South coast
SC.monitor = monitor[monitor$AIRBASIN %in% c("South Coast"),]

#reproject the data to a suitable projection
#utm used because of the scale of the analysis
SC.monitor.t = spTransform(SC.monitor, CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 + towgs84=0,0,0"))

plot(SC.tracts.t, border="darkgrey", main="Monitoring Stations Across Los Angeles")
plot(all.tracts.t, border="lightgrey", add=T)
plot(SC.tracts.t,border="darkgrey", add=T)
plot(SC.AirBasin.t,border="darkgrey", add=T)
plot(SC.monitor.t,add=T,col="red",cex=0.5,pch=16)

#mean centre
meanlat = mean(SC.monitor.t$LATITUDE)
meanlong = mean(SC.monitor.t$LONGITUDE)
meancentre = c(meanlat,meanlong)

#weighted mean centre not neccessary because there is only 1 entry per location

#standard distance
tempdata = as.data.frame(SC.monitor.t)
tempdata = cbind(tempdata, xmean=tempdata$LATITUDE-meanlat)
tempdata = cbind(tempdata, ymean=tempdata$LONGITUDE-meanlong)
sumX = sum((tempdata$xmean)**2)
sumY = sum((tempdata$ymean)**2)
standdist = sqrt((sumX+sumY)/length(tempdata$LATITUDE))

#relative distance
longdist = max(tempdata$LONGITUDE) - min(tempdata$LONGITUDE)
latdist = max(tempdata$LATITUDE) - min(tempdata$LATITUDE)
sqarea = longdist*latdist
circleradius = sqrt(sqarea/pi)
reldist = standdist/circleradius

#map for displaying values
#plot(meancentre,add=T,col="blue",cex=0.5,pch=16)
plot(SC.monitor.t)
draw.circle(-117.7768,33.9788,0.009,col="green",border="green")
draw.circle(-117.7768,33.9788,standdist,col=rgb(0,0,1,0.4),border="blue")
legend(x="bottomleft",legend=c("Monitoring Stations","Mean Centre","Standard Distance"),col=c("red","green","blue"),bg="white",pch=16,border="white",cex=0.8)

#quadrat analysis

#create a window (boundary) object to define study site
window = as.owin(SC.tracts.t)
#prepare data for point pattern analysis
monitors.ppp = ppp(x=tempdata$LONGITUDE,y=tempdata$LATITUDE,window=window)
#define number of quadrats used in analysis
quads = 8
#use fxn quadratcount to calc number of points per quadrats
qcount = quadratcount(monitors.ppp,nx=quads,ny=quads)
#create map showing quadrats and points per quadrats
jpeg("Quadrat16.jpeg",2500,2000,res=300)
plot(monitors.ppp,pch="+",cex=0.5,border="lightgrey",main="Quadrat Analysis of Monitoring Stations")
plot(qcount,add=T,col="blue",lwd=0.5)
#now prepare data to be able to calculate the variables for a quadrat analysis
#first deifne the quadrat count dataset as a dataframe
qcount.df = as.data.frame((qcount))
#count the number of quadrats with distinct number of points
qcount.df = count(qcount.df,'Freq')
#change column names
colnames(qcount.df) = c("x","f")
#create new columns for total number of points and for fx^2
qcount.df = cbind(qcount.df, TotPoints = as.numeric(qcount.df$x)*as.numeric(qcount.df$f))
qcount.df = cbind(qcount.df, fx2 = (as.numeric(qcount.df$x)^2)*as.numeric(qcount.df$f))
#calculate the sum of each column which you will use as inputs into the formula for vmr
f.sum = sum(qcount.df$f)
TotPoints.sum = sum(qcount.df$TotPoints)
fx2.sum = sum(qcount.df$fx2)
#calculate variance, mean, vmr
variance = (fx2.sum * (TotPoints.sum**2 / quads**2)) / (quads**2 -1) ##CHECK THIS
mean = TotPoints.sum / f.sum
vmr = variance/mean
#perform the test statistic to test for existence of a random spatial pattern
chisquare = chisq.test(qcount.df)

#nearest neighbour analysis
nearestNeighbour = nndist(monitors.ppp)
nearestNeighbour = as.data.frame(as.numeric(nearestNeighbour))
colnames(nearestNeighbour) = "Distance"
nnd = sum(nearestNeighbour$Distance)/length(nearestNeighbour$Distance)
#to calculate density, its freq/distance


#### SPATIAL AUTOCORRELATION ####

#lagged mean plots
#neighbour lists extracted
SC.tracts.nb = poly2nb((SC.tracts.t), queen=FALSE)
#convert to a listw object
SC.tracts.lw = nb2listw(SC.tracts.nb, zero.policy=TRUE)
diversity.lagged.means <- lag.listw(SC.tracts.lw,all.data$Diversity)
choropleth(SC.tracts.t,smk.lagged.means,shades)

moran.plot(SC.tracts.t$Diversity,SC.tracts.lw,xlab="Diversity",ylab="Spatially Lagged Diversity",main="Moran Scatter Plot for Diversity in Los Angeles")

#global
#z - (global i - expectation) / sqrt(variance) = moran i stat stand dev
moran.test(SC.tracts.t$Diversity,SC.tracts.lw,zero.policy=TRUE)

#local morans i
div.lI <- localmoran(SC.tracts.t$Diversity,SC.tracts.lw, zero.policy=TRUE)
#create shading scheme
div.shades <- auto.shading(c(div.lI[,1],-div.lI[,1]),cols=brewer.pal(5,"Blues"))
choropleth(SC.tracts.t,div.lI[,1],shading=div.shades,border=NA,main="Local Moran's I for Diversity in Los Angeles")
choro.legend("bottomleft",sh=div.shades,fmt="%6.2f",xpd=TRUE,cex=0.7)
#for p value
pval.shades <- shading(c(0.01,0.05,0.1),cols=brewer.pal(5,"Blues"))
choropleth(SC.tracts.t,div.lI[,5],shading=pval.shades,border=NA,main="Local p-value for Diversity in Los Angeles")
choro.legend("bottomleft",sh=pval.shades,fmt="%6.2f",xpd=TRUE,cex=0.7)

#test for normality
#matrix of simulated local morans i
sim.I <- matrix(0,1000,length(SC.tracts.t$Diversity))
#evaluate simulated distribution
for(i in 1:1000) sim.I[i,] <- localmoran(sample(SC.tracts.t$Diversity),SC.tracts.lw)[,4]
qqnorm(sim.I[,1],main="QQ Plot for Diversity in Los Angeles")
qqline(sim.I[,1])



#### SPATIAL INTERPOLATION ####

#to create pollution surface from monitoring station data, voronoi then kriging, as not enough monitoring stations for all census tracts to be represented

#### CORRELATION ANALYSIS ####

IV = all.data$mean.partmat
DV = all.data$Diversity
plot(DV,IV, main="Correlation between Diversity and PM2.5 in Los Angeles", xlab="Particulate Matter 2.5", ylab="Diversity")
cor.test(DV,IV)


#### REGRESSION ANALYSIS ####

lmfit = lm(IV~DV)
abline(lmfit)
linear.model = lm(IV~DV)
summary(linear.model)
resids = residuals(linear.model)

write.csv(SC.tracts.t,file="for-resids.csv")
for.resids = read.csv("for-resids.csv", header=T, sep=",")

tracts.df = as.data.frame(for.resids)
tracts.df = cbind(tracts.df[1],resids)
colnames(tracts.df)[2] = "residuals"
merge.tracts = merge(for.resids,tracts.df,by='X')

minlat = min(SC.monitor$LATITUDE - 0.1)
maxlat = max(SC.monitor$LATITUDE + 0.1)
minlong = min(SC.monitor$LONGITUDE - 0.1)
maxlong = max(SC.monitor$LONGITUDE + 0.1)

shades = auto.shading(as.numeric(merge.tracts$residuals),cols=brewer.pal(9,"Blues"))
choropleth(merge.tracts, merge.tracts$residuals, border=NA, shading=shades, xlim=c(minlong,maxlong), ylim=c(minlat,maxlat))

#global moran's i to see if residuals exhibit spatial autocorrelation
#weighted matrix created to perform moran's i analysis
SC.tracts.nb = poly2nb((merge.tracts), queen=FALSE)
SC.tracts.nb
#convert neighbour object to a listw object
SC.tract.lw = nb2listw(SC.tracts.nb, zero.policy = TRUE)
moran.test(resids,SC.tract.lw, zero.policy = TRUE)

all.df = as.data.frame(all.data)
all.spdf = SpatialPointsDataFrame(all.data,all.df)
#perform geographically weighted regression (only if spatially autocorrelated)
gwr.res = gwr.basic(as.numeric(DV)~as.numeric(IV), data=SC.tracts.t, bw=20000, kernel='gaussian')

#note the gwr coefficeient estimates and the r squared value
gwr.res
head(gwr.res$SDF)

quick.map <- function(spdf,var,legend.title,main.title) {
	x <- spdf@data[,var]
	cut.vals <- pretty(x)
	x.cut <- cut(x,cut.vals)
	cut.levels <- levels(x.cut)
	cut.band <- match(x.cut, cut.levels)
	colors <- brewer.pal(length(cut.levels), 'Reds')
	par(mar=c(1,1,1,1))
	plot(SC.AirBasin,col='grey85')
	title(main.title)
	plot(spdf,add=TRUE,col=colors[cut.band],pch=16)
	legend('topleft',cut.levels,col=colors,pch=16,bty='n',title=legend.title)
}
#using the geographically weighted regresion data, plot the gw regression coefficients
quick.map(gwr.res$SDF,"Particulate Matter 2.5","","GW Regression Coefficients")
#map of geographically weighted regression residuals
quick.map(gwr.res$SDF, "Residual","","Geographically Weighted Regression Residuals")

