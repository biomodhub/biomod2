#!/usr/bin/env r

## assign a basic time type
now <- Sys.time()

basicTypeTests <- function() {

    ## print invokes a conversion to char for 'display' even though the type is really POSIXct
    print(now)
    print(class(now))

    ## we can also convert to char() implicitly
    print(format(now))
    print(class(format(now)))

    ## but what is important is that 'now' is still a time type
    ## that we can 'compute'
    print(now)
    print(now + 60)  ## one minute later
    twomin <- now + 120
    print(as.numeric(now))

    ## and the time is even stored at millisecond granularity!!
    options("digits.secs"=7)		## need to tell R we want sub-second display up to 7 digits
    print(now)				## and now we do
    print(as.numeric(now), digits=16)	## same milli/microseconds here


    ## you can also go the other way and parse a datetime object from a char vector
    then <- strptime("2008-07-01 14:15:16", "%Y-%m-%d %H:%M:%S")
    print(then)
    print(class(then))
    ## and we can convert this from its default 'POSIXlt' ('long type) representation to 'POSIXct' ('compact type')
    then <- as.POSIXct(then)
    print(then)
    print(class(then))
}

dbTypeTests <- function(dateclass="timestamp without time zone") {
    cat("\n\n**** Trying with ", dateclass, "\n")
    tempdb <- "pgdatetime"
    system(paste("createdb", tempdb))   # create a temp database

    stopifnot(require(RPostgreSQL))
    drv <- dbDriver("PostgreSQL")
    con <- dbConnect(drv, dbname=tempdb)

    dbSendQuery(con, paste("create table foo (tt ", dateclass, ", zz integer);", sep=""))
    dbSendQuery(con, "insert into foo values('2008-07-01 14:15:16.123', 1);")

    dbSendQuery(con, paste("insert into foo values('", format(now), "', 2);", sep=""))

    res <- dbReadTable(con, "foo")  ## fails with 'RS-DBI driver warning: (unrecognized PostgreSQL field type 1184 in column 0)'
    print(res)

    ##res <- dbSendQuery(con, "select to_char(tt, 'YYYY-MM-DD HH24:MI:SS.US TZ') as character from foo;")
    res <- dbSendQuery(con, "select tt from foo;")
    data <- fetch(res, n=-1)
    #data <- dbGetQuery(con, "select tt from foo;")
    print(dbColumnInfo(res))
    ##times <- strptime(data[,1], "%Y-%m-%d %H:%M:%OS")
    times <- data[,1]
    print(times)
    print(class(times[1]))

    print(diff(times))	## yes we can compute on date times

    dbDisconnect(con)

    system(paste("dropdb", tempdb))   # create a temp database
}

dbTypeTests()
dbTypeTests("timestamp")
dbTypeTests("timestamp with time zone")
dbTypeTests("date")
#dbTypeTests("time with time zone")
#dbTypeTests("time without time zone")
