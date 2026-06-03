### R code from vignette source 'stsubs.Rnw'

###################################################
### code chunk number 1: stsubs.Rnw:32-57
###################################################
library(sp)
library(xts)
library(spacetime)

# spatial component
spObj <- data.frame(x=c(0,2,1,0,2), y=c(2,2,1,0,0))
coordinates(spObj) <- ~x+y
row.names(spObj) <- c("A","B","C","D","E")

# temporal component
timeObj = as.POSIXct("2010-08-05")+3600*24*(0:1)

# data slot
m = 1:5*10 # means for each of the 5 point locations
mydata = data.frame(smpl1=rnorm(10, mean=rep(m, 2)),
                    smpl2=rnorm(10, mean=rep(m, 2)))
mydata[c(2, 4, 6, 7),] <- NA

STFDFobj <- STFDF(spObj, timeObj, mydata)
STSDFobj <- as(STFDFobj, "STSDF")
STIDFobj <- as(STFDFobj, "STIDF")

STFobj <- geometry(STFDFobj)
STSobj <- geometry(STSDFobj)
STIobj <- geometry(STIDFobj)


###################################################
### code chunk number 2: stsubs.Rnw:67-68
###################################################
STFobj[1:3,]


###################################################
### code chunk number 3: stsubs.Rnw:72-73
###################################################
STSobj[1:3,]


###################################################
### code chunk number 4: stsubs.Rnw:77-78
###################################################
STIobj[1:3,]


###################################################
### code chunk number 5: stsubs.Rnw:82-83
###################################################
STSobj[rep(3,3),]


###################################################
### code chunk number 6: stsubs.Rnw:87-88
###################################################
STFobj[3:1,]


###################################################
### code chunk number 7: stsubs.Rnw:95-96
###################################################
STFobj[c(TRUE,TRUE,FALSE,FALSE,TRUE),]


###################################################
### code chunk number 8: stsubs.Rnw:100-101
###################################################
STSobj[c(TRUE,TRUE,FALSE,FALSE,TRUE),]


###################################################
### code chunk number 9: stsubs.Rnw:105-106
###################################################
STSobj[c(TRUE,TRUE,FALSE,FALSE,TRUE,FALSE),]


###################################################
### code chunk number 10: stsubs.Rnw:115-116
###################################################
STFobj[c("A","B"),]


###################################################
### code chunk number 11: stsubs.Rnw:120-121
###################################################
STSobj[c("A","B"),]


###################################################
### code chunk number 12: stsubs.Rnw:125-126
###################################################
STIobj[c("C"),]


###################################################
### code chunk number 13: stsubs.Rnw:130-131
###################################################
STIobj[c("C","C"),]


###################################################
### code chunk number 14: stsubs.Rnw:135-136
###################################################
STFobj[c("E","A"),]


###################################################
### code chunk number 15: stsubs.Rnw:153-154
###################################################
STFobj[cbind(1:3,c(1,2,2)),]


###################################################
### code chunk number 16: stsubs.Rnw:158-159
###################################################
STSobj[cbind(1:3,c(1,1,1))]


###################################################
### code chunk number 17: stsubs.Rnw:163-164
###################################################
STIobj[cbind(1:5,5:1)]


