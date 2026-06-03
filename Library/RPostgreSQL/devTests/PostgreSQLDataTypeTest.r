#!/usr/bin/env r

## Create a database
tempdb <- "tempdbase"
system(paste("createdb", tempdb))

library(RPostgreSQL)
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname=tempdb)



## Test the numeric mapping
dbSendQuery(con, "create table testnumeric (intcolumn integer, floatcolumn float);")

i <- as.integer(10)
j <- as.numeric(56.6)

sql <- paste("insert into testnumeric ",
            "values (",i, "," ,j ,") ", sep="")
res <- dbSendQuery(con, sql)


dat <- dbReadTable(con, "testnumeric")
cat("Read Numeric values\n")

## now test the types of the colums we got
stopifnot( class(dat[,1]) == "integer" )
stopifnot( class(dat[,2]) == "numeric" )
cat("GOOD -- all numeric types are as expected\n")

## and test the values
stopifnot( identical( dat[1,1], i))
stopifnot( identical( dat[1,2], j))
cat("GOOD -- all numeric values are as expected\n")




## Test the logical mapping
dbSendQuery(con,"create table testlogical (col1 boolean, col2 boolean)")

i <- as.logical(TRUE)
j <- as.logical(FALSE)

sql <- paste("insert into testlogical ",
            "values (",i, "," ,j ,") ", sep="")
res <- dbSendQuery(con, sql);

dat <- dbReadTable(con, "testlogical")
cat("Read Logical values\n")

## now test the types of the colums we got
stopifnot( class(dat[,1]) == "logical" )
stopifnot( class(dat[,2]) == "logical" )
cat("GOOD -- all logical types are as expected\n")

## and test the values
stopifnot( identical( dat[1,1], i))
stopifnot( identical( dat[1,2], j))
cat("GOOD -- all logical values are as expected\n")




## Test the character mapping
dbSendQuery(con,"create table testchar (code char(3),city varchar(20),country text);")

i <- as.character("IN")
j <- as.character("Hyderabad")
k <- as.character("India")


sql <- paste("insert into testchar ",
            "values ('",i,"' , '",j ,"' , '",k,"') ", sep="")
res <- dbSendQuery(con, sql);

dat <- dbReadTable(con, "testchar")
cat("Read Character values\n")

## now test the types of the colums we got
stopifnot( class(dat[,1]) == "character" )
stopifnot( class(dat[,2]) == "character" )
stopifnot( class(dat[,3]) == "character" )
cat("GOOD -- all character types are as expected\n")


## and test the values
##stopifnot( identical( dat[1,1], i))
stopifnot( identical( dat[1,2], j))
stopifnot( identical( dat[1,3], k))
cat("GOOD -- all character values are as expected\n")


dbDisconnect(con)
dbUnloadDriver(drv)

system(paste("dropdb", tempdb))

cat("DONE\n") 
