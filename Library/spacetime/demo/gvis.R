library(spacetime)
data(air)
#rural_PM10 = na.omit(as(rural[1:5,], "data.frame"))
#rural_PM10 = na.omit(as(na.locf(rural)[1:5,], "data.frame"))
#rural_PM10 = na.omit(as(rural[1:5,], "data.frame"))
rural_PM10 = as(rural[1:5,], "data.frame")
rural_PM10$PM10[is.na(rural_PM10$PM10)] = 0

library(googleVis)
# un-annotated timeLine:
TimeLine  <- gvisAnnotatedTimeLine(
	rural_PM10,
	datevar="time",
	numvar="PM10", 
	idvar="sp.ID",
	options=list(displayAnnotations=FALSE, width=900, height=600)
)
#plot(TimeLine)

rural_PM10$Annotation = rural_PM10$Title = as.character(NA)
row = which(rural_PM10$sp.ID == "DEBE056" & 
	rural_PM10$time == as.Date("2003-12-31"))
row
rural_PM10[row, "Title"] = "DEBE056"
rural_PM10[row, "Annotation"] = "Period with missing values drawn as line"
summary(rural_PM10)

AnnoTimeLine  <- gvisAnnotatedTimeLine(
	rural_PM10,
	datevar="time",
	numvar="PM10", 
	idvar="sp.ID",
	titlevar="Title", annotationvar="Annotation",
	options=list(displayAnnotations=TRUE, 
		zoomStartTime = as.Date("2003-07-01"), 
		zoomEndTime = as.Date("2004-07-01"),
		width=1200, height=600)
)
plot(AnnoTimeLine)

publish = FALSE
if (publish) {
	fname = paste(tempdir(), "/", AnnoTimeLine$chartid, ".html", sep="")
	target = "epebe_01@ifgifiles.uni-muenster.de:WWW/googleVis"
	scpcmd = paste("scp", fname, target)
	system(scpcmd)
}

# Questions / NOTES:
# does not deal well with NA values as response;
# how should gaps in the data (periods with NA values) be shown?

r = rural_PM10[1:100,]
r$PM10[20:80] = NA
TimeLine  <- gvisAnnotatedTimeLine(
	r,
	datevar="time",
	numvar="PM10", 
	idvar="sp.ID",
	options=list(displayAnnotations=FALSE, width=900, height=600)
)
plot(TimeLine)

# one single line chart with date:
LineChart = gvisLineChart(
	r, "time", "PM10"
)
plot(LineChart)

r2 = as.data.frame(as(rural[6:10,"2008"], "xts"))
r2$time = as.Date(rownames(r2))
LineChart = gvisLineChart(
	r2, "time", c("DEBE032", "DEHE046", "DEUB007", "DENW081", "DESH008"),
	options = list(width = 1200, focusTarget = "category",
	title = "PM10, 2008, for 5 German rural background stations",
	vAxis.logScale = TRUE) # vAxis.logScale doesn't seem to work...
)
plot(LineChart)

stopifnot(!is.projected(rural@sp))
sp = rural@sp
coord = coordinates(sp)
df = data.frame(cc = paste(coord[,2], coord[,1], sep=":"), 
	name = rownames(coord),
	stringsAsFactors = FALSE)
M2 <- gvisMap(df, "cc", "name",
              options=list(showTip = TRUE, mapType = 'normal',
              enableScrollWheel = TRUE, width = 1200, height = 700,
			  useMapTypeControl = TRUE))
plot(M2)


