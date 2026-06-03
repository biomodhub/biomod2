# R script that downloads and analyzes some of the examples in
# Noel Cressie, Christopher K. Wikle, 2011: Statistics for
# spatio-temporal data; Wiley, NY.

# data obtained (checked Fri May 20 10:37:33 CEST 2011) from:
# ftp://ftp.wiley.com/public/sci_tech_med/spatio_temporal_data

getFILE = function(f) { 
	if(file.exists(f)) 
		f
	else paste("ftp://ftp.wiley.com/public/sci_tech_med/spatio_temporal_data/",
	  	f, sep="") 
}
# read locations:
ecd.ll = as.matrix(read.table(getFILE("ECDovelatlon.dat"), header=FALSE))
# space: convert to SpatialPoints, and set CRS:
library(sp)
ecd.ll = SpatialPoints(ecd.ll[,c(2,1)]) # lat-lon -> lon-lat
proj4string(ecd.ll) = "+proj=longlat +datum=WGS84 +ellps=WGS84"

# time:
ecd.years = 1986:2003
ecd.y = as.Date(paste(ecd.years, "-01-01", sep=""), "%Y-%m-%d")

# spacetime:
# read, and convert to matrix
ecd = as.matrix(read.table(getFILE("ECDoveBBS1986_2003.dat"), header=FALSE))
# set missing values
ecd[ecd == -1] = NA
library(spacetime)
ecd.st = STFDF(ecd.ll, ecd.y, data.frame(counts = as.vector(ecd)))
sel = (1:6) * 3
stplot(ecd.st[,sel], 
	ecd.years[sel], 
	col.regions=bpy.colors(),
	key.space = "right",
	main = "Figure 5.5"
)

ecd.xts = as(ecd.st, "xts")
ecd.yt = (apply(ecd.xts, 1, function(x)sum(x,na.rm=TRUE)))
plot(index(ecd.xts), ecd.yt, type = 'b',
	xlab = 'Year', ylab = 'Eurasian Collared Dove Count')
title("Figure 5.3")

sst = as.matrix(read.table(getFILE("SST011970_032003.dat"), header = FALSE))

# space:
sst.ll = as.matrix(read.table(getFILE("SSTlonlat.dat"), header = FALSE))
sst.ll = SpatialPoints(sst.ll)
proj4string(sst.ll) = "+proj=longlat +datum=WGS84"
gridded(sst.ll) = TRUE

# time:
sst.y = rep(1970:2003, each = 12, length.out = 399)
sst.m = rep(1:12, length.out = 399)
# dates:
sst.Date = as.Date(paste(sst.y, sst.m, "01", sep="-"), "%Y-%m-%d")
# as year/months:
library(zoo)
sst.ym = as.yearmon(sst.Date)

# landmask:
sea = (as.vector(read.table(getFILE("SSTlandmask.dat"), header=FALSE)) != 1)[,1]

# spacetime:
sst.st = STFDF(sst.ll, sst.ym, data.frame(sst = as.vector(sst)))

# make a large plot window for the following plot:
tsel = 1:144
stplot(sst.st[sea, tsel], 
	format(sst.Date, "%Y-%m")[tsel],
	layout = c(12,12))

# re-create figure 5.1, with time axis increasing upward: 
library(rgeos) # needed for "over" method, to select pixels on a line;
eq = SpatialLines(list(Lines(list(Line(cbind(c(0,360),c(0,0)))), "equator")),
	CRS(proj4string(sst.st)))
stplot(sst.st[eq, "1996::2003", drop=F], 
	mode = "xt",
	col.regions=bpy.colors(), cuts = 32,
	main="Figure 5.1",
	xlab="Longitude", ylab = "Time (Year)",
	scaleX=1
)
# see also hovmoller() at http://rastervis.r-forge.r-project.org/

# re-create Fig 5.4:
stplot(sst.st[sea,"1998-02::1999-01"],
	col.regions=bpy.colors(),
	# funny order; as.table=TRUE by default, so:
	index.cond=list(c(1,7,2,8,3,9,4,10,5,11,6,12)),
	layout = c(2, 6), scales = list(draw = TRUE),
	main = "Figure 5.4"
)

# Figure 5.17: EOF's
eof = eof(sst.st[sea,])
spplot(eof[1], col.regions = bpy.colors(), scales = list(draw=TRUE),
	main = "First EOF (Figure 5.17a)")
spplot(eof[2], col.regions = bpy.colors(), scales = list(draw=TRUE),
	main = "5.17c; Second EOF")

# 5.19: EOF summary stats
eof.summ = eof(sst.st[sea,], returnEOFs = FALSE)

library(xts)
eof.t = xts(predict(eof.summ), sst.ym)
plot(eof.t[,1], main = "5.17b (note that sign has flipped)")
plot(eof.t[,2], main = "5.17d")
# ... and so on.

v = eof.summ$sdev^2
plot(100*cumsum(v[1:100])/sum(v),
	ylim=c(30, 100), ylab = "Percent", xlab = "EOF", main = "Figure 5.19")

# canton data:
ca_nb = read.csv(getFILE("Canton_neighbor.csv"))
ca_ve = read.csv(getFILE("Canton_vertex.csv"))

xc = paste("X", 1:76, sep="")
yc = paste("Y", 1:76, sep="")
p = lapply(1:nrow(ca_ve), function(r) {
		pol = (na.omit(cbind(as.numeric(ca_ve[r,xc]), 
			as.numeric(ca_ve[r,yc]))))
		if (!identical(pol[1,],pol[nrow(pol),]))
			pol = rbind(pol, pol[1,])
		Polygon(pol)
   })
srl = lapply(1:268, function(i) {
		Polygons(p[i], as.character(ca_ve[i, "INSEE_ID"]))
   })
rownames(ca_nb) = ca_nb$INSEE_ID
can = SpatialPolygonsDataFrame(SpatialPolygons(srl), ca_nb[,1:9])
q = quantile(can$Z, c(0,.05,.20,.35,.65,.80,.95,1))
can$Zcl = factor(findInterval(can$Z, q, all.inside = TRUE),
	labels = c("< 5%", "5%-20%", "20%-35%", "35%-65%", "65%-80%", 
	"80%-95%", "> 95%"))
library(RColorBrewer)
# get missing polygons, paint them green:
srl.mv = lapply(269:300, function(i) {
		Polygons(p[i], as.character(ca_ve[i, "NO"]))
   })
mv = SpatialPolygons(srl.mv)
mv.col = "#80FF80"
lt = list("sp.polygons", mv, fill = mv.col)
# info of fig 4.11 (with some differences!) --
spplot(can["Zcl"], sp.layout = lt,
	col.regions= rev(brewer.pal(7, "RdBu")), cuts=6) 

# get neighbours:
nb = apply(ca_nb[paste("N", 1:12, sep="")], 1, function(x) x[x!=0])
# convert neighbour index to 1:268
nb = lapply(nb, function(x) which(row.names(can) %in% x))

plot(can)
plot(mv, col = mv.col, add = TRUE)
x = coordinates(can)
row.names(x) = row.names(can)
#for(i in 1:268) {
#	for (j in 1:length(nb[[i]]))
#		lines(rbind(x[i,], x[as.character(nb[[i]][j]),]),col=grey(.5))
#		lines(rbind(x[i,], x[nb[[i]][j],]),col=grey(.5))
#}
library(spdep)
class(nb) = "nb"
plot(nb, coordinates(can), col = grey(.5), add=TRUE)
lw = nb2listw(nb)
# need to figure out whether the output of this is nonsense:
spautolm(Z~X1+X2, as.data.frame(can), lw)

s = sst.st[,1:3]
stplot(s)
sp = sst.st[,1]
class(sp)
sp.x = spsample(sp, 630, "regular", offset = c(.5,.5))
stplot(s, sp.layout=list("sp.points", sp.x))
gridded(sp.x)=TRUE
sa = aggregate(s, sp.x, mean)
stplot(sa)
s1 = aggregate(s, sp.x[sample(630),], mean)
s2 = aggregate(s, sp.x[sample(630),], mean)
stplot(s1)
s1[1,1]
stplot(s2) # plot should be identical
# should be different, as pixel order changed:
s2[1,1]
