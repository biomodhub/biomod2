#!/bin/bash

set -e

dbname="pgtemp"

createdb ${dbname} 
echo "Created ${dbname}"

cat <<EOF | R --slave
  library(RPostgreSQL)
  data(trees)
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, dbname="${dbname}")
  res <- dbWriteTable(con, "trees", trees)
  dbDisconnect(con)
EOF

psql ${dbname} -c "select * from trees;" 

dropdb ${dbname}
echo "Deleted ${dbname}"
