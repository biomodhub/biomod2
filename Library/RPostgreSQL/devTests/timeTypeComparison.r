#!/usr/bin/r

suppressMessages({
    library(RODBC)
    library(RPostgreSQL)
})

testTimeData <- function(tt = "date") {

    createTempTable(tt)

    channel <- odbcConnect("beancounter")
    odbcres <- sqlQuery(channel, "select * from timetest")
    close(channel)

    con <- dbConnect( dbDriver("PostgreSQL"), dbname="beancounter")
    pgres <- dbGetQuery(con, "select * from timetest")
    dbDisconnect(con)

    removeTempTable()

    cat("\nFor type='", tt, "'\nODBC returns ", format(odbcres[1,1]), " class ", class(odbcres[1,1]),
        "\nRPostgreSQL returns ", format(pgres[1,1]), " class ", class(pgres[1,1]), "\n", sep="")
}

createTempTable <- function(tt) {
    system(paste("psql beancounter -c \"create table timetest ( d1", tt, ")\" >/dev/null"))
    system("psql beancounter -c \"insert into timetest values( \'now\' )\" >/dev/null")
}

removeTempTable <- function() {
    system("psql beancounter -c \"drop table timetest\" >/dev/null")
}

testTimeData("date")
testTimeData("timestamp with time zone")
testTimeData("timestamp without time zone")
testTimeData("time with time zone")
testTimeData("time without time zone")
