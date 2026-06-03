## DEMO for illustrating all the commands in RPostgreSQL package.
## Run this file using 'cat demo.r | R --slave' at the command prompt.

## Create a database
tempdb <- "tempDEMO1dbase"
system(paste("dropdb", tempdb))
system(paste("createdb", tempdb))

cat("\nFor loading the library: library(RPostgreSQL)\n")
library(RPostgreSQL)

cat("\n## Demo for the working of methods in RPostgreSQL\n\n")

cat("##  1.dbDriver(\"driverName\") instantiates the driver object \nEg. drv <- dbDriver(\"PostgreSQL\") \n")
drv <- dbDriver("PostgreSQL")
drv

cat("\n## 2.dbConnect(drv,...) creates and opens a connection to the database implemented by the driver drv. Connection string should be specified with parameters like user, password, dbname, host, port, tty and options. For more details refer to the documentation \nEg. con <- dbConnect(drv,dbname=\"dbName\")\n")
con <- dbConnect(drv, dbname=tempdb)
con
con2<- dbConnect(drv, dbname=tempdb)
con2

cat("\n## 3.dbListConnection(drv, ...): List of connections handled by the driver \nEg. dbListConnections(drv)\n")
dbListConnections(drv)

cat("##  4. dbGetInfo(dbObject, ...) returns information about the dbObject like driver, connection or resultSet \n Eg. dbGetInfo(drv)\n")
dbGetInfo(drv)

cat("## 5. dbSendQuery(con, statement, ...) submits one statement to the database\n Eg. rs <- dbSendQuery(con,\"select * from TableName\")\n")
dbSendQuery(con, "create table book789listCosts (intcolumn integer, floatcolumn float);")
## insert four rows into the table
dbSendQuery(con, "insert into book789listCosts values(12,13.21);")
dbSendQuery(con, "insert into book789listCosts values(50,11.21);")
dbSendQuery(con, "insert into book789listCosts values(100,200.1);")
dbSendQuery(con, "insert into book789listCosts values(5,3.56);")
rs <- dbSendQuery(con, "select * from book789listCosts")

cat("\n## 6. fetch(rs,n, ...) fetches the next n elements from the result set. n=-1 means RETURN ALL ELEMENTS.\n Eg.fetch(rs,n=-1)\n")
fetch(rs,n=-1)

cat("\n## 7. dbGetQuery(con,statement, ...) submits, execute, and extract output in one operation. \nEg.dbGetQuery(con,\"select * from TableName\")\n")
dbGetQuery(con,"select * from book789listCosts")

cat("\n## 8. dbGetException(con, ...) returns the status of the last DBMS statement sent over the connection. \n Eg. dbGetException(con)\n")
dbGetException(con)

cat("\n## 9. dbListResults(con, ...) returns the resultsets active on the given connection. Please note that the current RPostgreSQL package can handle only one resultset per connection (which may change in the future).\n Eg. dbListResults(con)\n")
dbListResults(con)

cat("\n## 10. dbListTables(con, ...) returns the list of tables available on the connection. \n Eg. dbListTables(con)\n")
dbListTables(con)

cat("\n## 11. dbExistsTable(con, TableName, ...) checks whether a particular table exists on the given connection. Returns a logical.\n Eg. dbExistsTable(con,\"name\")\n")
dbExistsTable(con,"name")

cat("\n## 12. dbRemoveTable(con, TableName, ...) removes the specified table on the connection. Returns a logical indicating operation succeeded or not. \n Eg. dbRemoveTable(con,\"name\")\n")

cat("\n## 13. dbListFields(con, TableName, ...) returns the list of column names (fields) in the table.\n Eg. dbListFields(con,\"book789listCosts\")\n")

cat("\n## 14. dbColumnInfo(res, ...) produces a query that describes the output of the query.\n Eg. dbColumnInfo(rs)\n")
dbColumnInfo(rs)

cat("\n## 15. dbReadTable(conn, name, ...) imports the data stored remotely in the table name on connection conn. Use the ﬁeld row.names as the row.names attribute of the output data.frame. Returns a data.frame.\n Eg. dframe <- dbReadTable(con,\"book789listcosts\")\nAfter executing dbReadTable, the dframe contains:\n")
dframe <- dbReadTable(con,"book789listcosts")
dframe

cat("\n## 16. dbWriteTable(conn, name, value, ...) writes the contents of the dataframe \'value\' into the table name specified. Returns a logical indicating whether operation succeeded or not. \n Eg. dbWriteTable(con,\"newTable\",dframe)\n")
dbWriteTable(con,"newTable",dframe)
cat("Checking the contents of newTable: \n Eg. dbGetQuery(con,\"newTable\") \n")
dbGetQuery(con,"select * from newTable")

cat("\n## 17. dbGetStatement(res, ...) returns the DBMS statement associated with the result.\n Eg. dbGetStatement(rs)\n")
dbGetStatement(rs)

cat("\n## 18. dbGetRowsAffected(res, ...) returns the number of rows affected the executed statement. If no rows are affected, \"-1\" is returned. \n Eg. dbGetRowsAffected(rs)\n")
dbGetRowsAffected(rs)

cat("\n## 19. dbHasCompleted(res, ...) returns a logical to indicate whether an operation is completed or not\n Eg. dbHasCompleted(rs)\n")
dbHasCompleted(rs)

cat("\n## 20.dbGetRowCount(res, ...) returns number of rows fetched so far.\n Eg. dbGetRowCount(rs)\n")
dbGetRowCount(rs)

cat("\n## Commands for Transaction Management \n")

cat("\n## 21.dbBeginTransaction begins the PostgreSQL transaction. dbCommit commits the transaction while dbRollback rolls back the transaction. Returns a logical indicating whether the operation succeeded or not.\n")
cat("Eg1. dbBeginTransaction(con) \n")
dbBeginTransaction(con)
cat("Eg1. dbRemoveTable(con,\"newTable\") \n")
dbRemoveTable(con,"newTable")
cat("Eg1. dbExistsTable(con,\"newTable\") \n")
dbExistsTable(con,"newTable")
cat("Eg1. dbRollback(con)\n")
dbRollback(con)
cat("Eg1. dbExistsTable(con,\"newTable\") \n")
dbExistsTable(con,"newTable")
cat("\nEg2. dbBeginTransaction(con) \n")
dbBeginTransaction(con)
cat("Eg2. dbRemoveTable(con,\"newTable\") \n")
dbRemoveTable(con,"newTable")
cat("Eg2. dbExistsTable(con,\"newTable\") \n")
dbExistsTable(con,"newTable")
cat("Eg2. dbCommit(con)\n")
dbCommit(con)
cat("Eg2. dbExistsTable(con,\"newTable\") \n")
dbExistsTable(con,"newTable")


cat("\n## Commands for Freeing the Resources \n")

cat("\n## 22. dbClearResult(rs, ...) flushes any pending data and frees the resources used by result set.\n Eg.dbClearResult(resultSet)\n Cleared the ResultSet :")
dbClearResult(rs)

cat("\n## 23. dbDisconnect(con, ...) closes the connection. \nEg. dbDisconnect(con)")
cat("\nClosed the connections:\n")
dbDisconnect(con)
dbDisconnect(con2)

cat("\n## 24. dbUnloadDriver(drv,...) frees all the resources used by the driver\n Eg. dbUnloadDriver(drv)\n")
cat("Freed the Resources: ")
dbUnloadDriver(drv)

cat("\n")
system(paste("dropdb", tempdb))
cat("Demo Completed\n")
