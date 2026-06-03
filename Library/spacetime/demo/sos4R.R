library(gstat)

# move to an appropriate CRS:
RD = CRS(paste("+init=epsg:28992 +ellps=bessel",
 "+towgs84=565.237,50.0087,465.658,-0.406857,0.350733,-1.87035,4.0812"))

# in the current sos4R version, no2.spdf has reversed x and y. swap them:
x = as.data.frame(no2.spdf)
coordinates(x)=~lon+lat
no2.spdf=x
proj4string(no2.spdf) = proj4string(map.lines)

m = map2SpatialLines(map("worldHires", "Netherlands", plot=FALSE))
proj4string(m) = "+proj=longlat +datum=WGS84"

library(rgdal)
no2.spdf = spTransform(no2.spdf, RD)
m = spTransform(m, RD)

samplingTimes = sort(unique(no2.spdf$SamplingTime))
summary(samplingTimes)
no2.T1 = no2.spdf[no2.spdf$SamplingTime == samplingTimes[1],]

grd = SpatialPixels(SpatialPoints(makegrid(bbox(m), n = 10000)),
	proj4string = RD)

plot(grd)
plot(m, add=T)

names(no2.T1)[3] = "NO2"
NO2.idw =idw(NO2~1, no2.T1, grd)
lt = list(list("sp.lines", m, col = 'white'),
	list("sp.points", no2.T1, col = grey(.5)))
spplot(NO2.idw[1], col.regions = bpy.colors(), sp.layout=lt)

# exercise:
# 1. repeat the same procedure for the observation of the second sampling time
# 2. repeat the procedure for the observations of the last sampling time
