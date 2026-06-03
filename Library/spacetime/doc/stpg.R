### R code from vignette source 'stpg.Rnw'

###################################################
### code chunk number 1: stpg.Rnw:47-48 (eval = FALSE)
###################################################
## library(RPostgreSQL)


###################################################
### code chunk number 2: stpg.Rnw:50-52
###################################################
library(sp)
library(spacetime)


###################################################
### code chunk number 3: stpg.Rnw:62-66
###################################################
dbname = "postgis"
user = "edzer"
password = "pw"
#password = ""


###################################################
### code chunk number 4: stpg.Rnw:70-73 (eval = FALSE)
###################################################
## drv <- dbDriver("PostgreSQL")
## con <- dbConnect(drv, dbname=dbname, user=user, password=password,
## host='localhost', port='5432')


###################################################
### code chunk number 5: stpg.Rnw:81-85 (eval = FALSE)
###################################################
## dbRemoveTable(con, "rural_attr")
## dbRemoveTable(con, "rural_space")
## dbRemoveTable(con, "rural_time")
## dbRemoveTable(con, "space_select")


###################################################
### code chunk number 6: stpg.Rnw:92-102 (eval = FALSE)
###################################################
## data(air)
## rural = STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
## rural = as(rural, "STSDF")
## p = rural@sp
## sp = SpatialPointsDataFrame(p, data.frame(geom_id=1:length(p)))
## library(rgdal)
## OGRstring = paste("PG:dbname=", dbname, " user=", user, 
## 	" password=", password, " host=localhost", sep = "")
## print(OGRstring)
## writeOGR(sp, OGRstring, "rural_space", driver = "PostgreSQL")


###################################################
### code chunk number 7: stpg.Rnw:107-108 (eval = FALSE)
###################################################
## subset(ogrDrivers(), name == "PostgreSQL")$write


###################################################
### code chunk number 8: stpg.Rnw:114-118 (eval = FALSE)
###################################################
## df = data.frame(time = index(rural@time), time_id = 1:nrow(rural@time))
## dbWriteTable(con, "rural_time", df)
## idx = "create index time_idx on rural_time (time);"
## dbSendQuery(con, idx)


###################################################
### code chunk number 9: stpg.Rnw:123-127 (eval = FALSE)
###################################################
## idx = rural@index
## names(rural@data) = "pm10" # lower case
## df = cbind(data.frame(geom_id = idx[,1], time_id = idx[,2]), rural@data)
## dbWriteTable(con, "rural_attr", df)


###################################################
### code chunk number 10: stpg.Rnw:136-143
###################################################
setClass("ST_PG", contains = "ST", 
	# slots = c(space_table = "character",
	representation(space_table = "character",
	time_table = "character",
	attr_table = "character",
	attr = "character",
	con = "PostgreSQLConnection"))


###################################################
### code chunk number 11: stpg.Rnw:146-154 (eval = FALSE)
###################################################
## rural_proxy = new("ST_PG", 
## 	#ST(rural@sp, rural@time, rural@endTime),
## 	as(rural, "ST"),
## 	space_table = "rural_space",
## 	time_table = "rural_time",
## 	attr_table = "rural_attr",
## 	attr = "pm10",
## 	con = con)


###################################################
### code chunk number 12: stpg.Rnw:161-181
###################################################
.SqlTime = function(x, j) {
	stopifnot(is.character(j))
	require(xts)
	t = .parseISO8601(j)
	t1 = paste("'", t$first.time, "'", sep = "")
	t2 = paste("'", t$last.time, "'", sep = "")
	what = paste("geom_id, time_id", paste(x@attr, collapse = ","), sep = ", ")
	paste("SELECT", what, "FROM", x@attr_table, "AS a JOIN", x@time_table,
		"AS b USING (time_id) WHERE b.time >= ", t1, "AND b.time <=", t2,";")
}
.SqlSpace = function(x, i) {
	stopifnot(is(i, "Spatial"))
	writeOGR(i, OGRstring, "space_select", driver = "PostgreSQL")
	what = paste("geom_id, time_id", paste(x@attr, collapse = ","), sep = ", ")
	paste("SELECT",  what, "FROM",  x@attr_table, 
		"AS a JOIN (SELECT p.wkb_geometry, p.geom_id FROM",
		x@space_table, " AS p, space_select AS q",
		"WHERE ST_Intersects(p.wkb_geometry, q.wkb_geometry))",
		"AS b USING (geom_id);")
}


###################################################
### code chunk number 13: stpg.Rnw:187-197
###################################################
setMethod("[", "ST_PG", function(x, i, j, ... , drop = TRUE) {
	stopifnot(missing(i) != missing(j)) # either of them present
	if (missing(j))
		sql = .SqlSpace(x,i)
	else
		sql = .SqlTime(x,j)
	print(sql)
	df = dbGetQuery(x@con, sql)
	STSDF(x@sp, x@time, df[x@attr], as.matrix(df[c("geom_id", "time_id")]))
})


###################################################
### code chunk number 14: stpg.Rnw:199-206 (eval = FALSE)
###################################################
## pm10_20050101 = rural_proxy[, "2005-01-01"]
## summary(pm10_20050101)
## summary(rural[,"2005-01-01"])
## 
## pm10_NRW = rural_proxy[DE_NUTS1[10,],]
## summary(pm10_NRW)
## summary(rural[DE_NUTS1[10,],])


###################################################
### code chunk number 15: stpg.Rnw:211-214 (eval = FALSE)
###################################################
## dim(pm10_NRW)
## pm10_NRW = pm10_NRW[T,]
## dim(pm10_NRW)


###################################################
### code chunk number 16: stpg.Rnw:217-220 (eval = FALSE)
###################################################
## object.size(rural)
## object.size(pm10_20050101)
## object.size(pm10_NRW)


###################################################
### code chunk number 17: stpg.Rnw:226-228 (eval = FALSE)
###################################################
## dbDisconnect(con)
## dbUnloadDriver(drv)


