# exercises for ifgi - Esri development center workshop
# held Muenster, Mar 17-18 2011

###################################################
### chunk number 26: 
###################################################
#line 577 "spacetime.Rnw"

# here is an example where data in long format that come
# as a table in package plm, are combined (through USA state names)
# with state boundaries in package maps:

# first obtain the state boundaries, and their IDs (state names):
library(maps)
states.m = map('state', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(states.m$names, ":"), function(x) x[1])
IDs
    
# then, convert the states to a SpatialPolygons object:
library(maptools)
states = map2SpatialPolygons(states.m, IDs=IDs)
plot(states)

# next, import the Produc data, a table with space and time in
# columns 1 and 2:
library(plm)
data(Produc)
summary(Produc)

# then, construct the STFDF object:
library(spacetime)
yrs = 1970:1986
time = xts(1:17, as.POSIXct(paste(yrs, "-01-01", sep=""), tz = "GMT"))
# deselect District of Columbia, polygon 8, which is not present in Produc:
Produc.st = STFDF(states[-8], time, Produc[(order(Produc[2], Produc[1])),])

# and plot it:
stplot(Produc.st[,,"unemp"], yrs)


###################################################
### chunk number 28: 
###################################################
#line 638 "spacetime.Rnw"

# spatio-temporal interpolation: the Irish wind data set
# first: import the data:
library(gstat)
data(wind)
# this imports two data.frame's: wind.loc and wind:
summary(wind.loc)
summary(wind)
wind.loc$y = as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
wind.loc$x = as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))
coordinates(wind.loc) = ~x+y
proj4string(wind.loc) = "+proj=longlat +datum=WGS84"
# now wind.loc is converted into a georeferenced SpatialPointsDataFrame object:
summary(wind.loc)


###################################################
### chunk number 29: 
###################################################
#line 649 "spacetime.Rnw"

# load a map reference for these data, and plot the wind stations:
library(mapdata)
plot(wind.loc, xlim = c(-11,-5.4), ylim = c(51,55.5), axes=T, col="red",
	cex.axis =.7)
map("worldHires", add=T, col = grey(.5))
text(coordinates(wind.loc), pos=1, label=wind.loc$Station, cex=.7)


###################################################
### chunk number 30: 
###################################################
#line 666 "spacetime.Rnw"

# check the wind data: they are in space-wide form: stations form columns
wind[1:3,]


###################################################
### chunk number 31: 
###################################################
#line 671 "spacetime.Rnw"

# let us reference the time information into a time class:
wind$time = ISOdate(wind$year+1900, wind$month, wind$day)
summary(wind$time)

# the next few commands are some statistical magic; wind is sqrt-transformed
# and a smooth pattern is subtracted from the raw numbers:
wind$jday = as.numeric(format(wind$time, '%j'))
stations = 4:15
windsqrt = sqrt(0.5148 * wind[stations]) # knots -> m/s
Jday = 1:366
daymeans = apply(sapply(split(windsqrt - mean(windsqrt), wind$jday), mean), 2,  mean)
meanwind = lowess(daymeans ~ Jday, f = 0.1)$y[wind$jday]
velocities = apply(windsqrt, 2, function(x) { x - meanwind })

###################################################
### chunk number 32: 
###################################################
#line 685 "spacetime.Rnw"

# order locations to order of columns in wind;
# connect station names to location coordinates
wind.loc = wind.loc[match(names(wind[4:15]), wind.loc$Code),]
pts = coordinates(wind.loc[match(names(wind[4:15]), wind.loc$Code),])
rownames(pts) = wind.loc$Station
pts = SpatialPoints(pts)

# convert to utm zone 29, to be able to do interpolation in
# proper Euclidian (projected) space:
library(rgdal)
proj4string(pts) = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")
utm29 = CRS("+proj=utm +zone=29 +datum=WGS84 +ellps=WGS84")
t = xts(1:nrow(wind), wind$time)
pts = spTransform(pts, utm29)
# create the spacetime object:
# note the t() in:
w = STFDF(pts, t, data.frame(values = as.vector(t(velocities))))

# import, convert, and project the backdrop (country boundaries) data too:
library(maptools)
m = map2SpatialLines(
	map("worldHires", xlim = c(-11,-5.4), ylim = c(51,55.5), plot=F))
proj4string(m) = "+proj=longlat +datum=WGS84"
m = spTransform(m, utm29)

# setup a grid covering the region of m:
grd = SpatialPixels(SpatialPoints(makegrid(m, n = 300)),
	proj4string = proj4string(m))
# select data (arbitrarily) for april 1961:
w = w[, "1961-04"]
# 10 prediction time points, evenly spread over this month:
n = 10
tgrd = xts(1:n, seq(min(index(w)), max(index(w)), length=n))

# using a separable exponential covariance model, 
# with ranges 750 km and 1.5 day, do a space-time interpolation:
v = list(space = vgm(0.6, "Exp", 750000), time = vgm(1, "Exp", 1.5 * 3600 * 24))
pred = krigeST(sqrt(values)~1, w, STF(grd, tgrd), v)
wind.ST = STFDF(grd, tgrd, data.frame(sqrt_speed = pred))

###################################################
### chunk number 33: 
###################################################
#line 727 "spacetime.Rnw"

# plot the results:
layout = list(list("sp.lines", m, col='grey'),
	list("sp.points", pts, first=F, cex=.5))
stplot(wind.ST, col.regions=bpy.colors(),
	par.strip.text = list(cex=.5), sp.layout = layout)

# Exercises:
# 1. repeat the procedure but interpolate the daily 12:00 values in May 1962
# 2. plot, from the results obtained, only the first three days
# 3. repeat the procedure, but use a grid with 1000 grid cells instead of 300.
