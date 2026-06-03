#!/usr/bin/env r

usePG <- TRUE

if (usePG) {
    cat("Using Pg\n")
    tempdb <- "pgdatetime"
    system(paste("createdb", tempdb))   # create a temp database

    suppressMessages(library(RPostgreSQL))
    drv <- dbDriver("PostgreSQL")
    con <- dbConnect(drv, dbname=tempdb)

} else {
    cat("Using SQLite\n")
    tempdb <- "/tmp/tempdb.sqlite"          # assumed to not exist
    suppressMessages(library(RSQLite))
    drv <- dbDriver("SQLite")
    if (file.exists(tempdb))
        unlink(tempdb)
    con <- dbConnect(drv, tempdb)
}


sql <- "create table foo (i integer, r real, t text)"
res <- dbSendQuery(con, sql)
cat("Created table\n")

i <- as.integer(11)
r <- as.numeric(22.22)
txt <- as.character("blim blom")

sql <- paste("insert into foo ",
             "values (",
             i, ",",
             r, ", '",
             txt, "') ", sep="")
res <- dbSendQuery(con, sql)
cat("Wrote values\n")

df <- dbReadTable(con, "foo")
cat("Read values\n")

## now test the types of the colums we got
stopifnot( class(df[,1]) == "integer" )
stopifnot( class(df[,2]) == "numeric" )
stopifnot( class(df[,3]) == "character" )
cat("GOOD -- all types are as expected\n")

## and test the values
stopifnot( identical( df[1,1], i))
stopifnot( identical( df[1,2], r))
stopifnot( identical( df[1,3], txt))
cat("GOOD -- all values are as expected\n")

if (usePG) {
    dbDisconnect(con)
    dbUnloadDriver(drv)
    system(paste("dropdb", tempdb))
} else {
    dbDisconnect(con)
    dbUnloadDriver(drv)
    unlink(tempdb)
}

cat("DONE\n")
