#!/usr/bin/env r

cat("Testing the working of Transaction Management\n")

## Create a database
tempdb <- "tempdbase123"
system(paste("createdb", tempdb))

library(RPostgreSQL)
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname=tempdb)

## Test the numeric mapping
dbSendQuery(con, "create table book123list (intcolumn integer, floatcolumn float);")

## insert four rows into the table
dbSendQuery(con, "insert into book123list values(12,13.21);")
dbSendQuery(con, "insert into book123list values(50,11.21);")
dbSendQuery(con, "insert into book123list values(100,200.1);")
dbSendQuery(con, "insert into book123list values(5,3.56);")

cat("Table book123list contains the following records\n")
dbGetQuery(con, "select * from book123list")

cat("Test Run 1:\n")

dbBeginTransaction(con);
cat("Begin Transaction\n")
rs <- dbSendQuery(con, "DELETE from book123list WHERE intcolumn >= 50");
dbClearResult(rs);
cat("After Deletion\n");
dbGetQuery(con, "select * from book123list")
dbRollback(con)
cat("After Rolling back\n")

cat("Table book123list contains the following records\n")
dbGetQuery(con, "select * from book123list")

cat("Test Run 2:\n")

dbBeginTransaction(con);
cat("Begin Transaction\n")
dbGetQuery(con, "select * from book123list")[1, ];
rs <- dbSendQuery(con, "DELETE from book123list WHERE intcolumn >= 50");
dbClearResult(rs);
cat("After Deletion\n");
dbGetQuery(con, "select * from book123list")
dbCommit(con);
cat("After commiting the transaction\n")
dbGetQuery(con, "select * from book123list")

dbDisconnect(con)
dbUnloadDriver(drv)

system(paste("dropdb", tempdb))

cat("DONE\n") 
