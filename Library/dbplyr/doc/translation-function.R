## -----------------------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)

## -----------------------------------------------------------------------------
library(dbplyr)
library(dplyr)

con <- simulate_dbi()

## -----------------------------------------------------------------------------
translate_sql((x + y) / 2, con = con)

## -----------------------------------------------------------------------------
translate_sql(x ^ 2L, con = con)
translate_sql(x ^ 2L, con = simulate_sqlite())
translate_sql(x ^ 2L, con = simulate_access())

## -----------------------------------------------------------------------------
# In SQLite variable names are escaped by double quotes:
translate_sql(x, con = con)
# And strings are escaped by single quotes
translate_sql("x", con = con)

## -----------------------------------------------------------------------------
translate_sql(substr(x, 5, 10), con = con)
translate_sql(log(x, 10), con = con)

## -----------------------------------------------------------------------------
translate_sql(1, con = con)
translate_sql(1L, con = con)

## -----------------------------------------------------------------------------
df <- tibble(
  x = c(10L, 10L, -10L, -10L), 
  y = c(3L, -3L, 3L, -3L)
)
mf <- tbl_memdb(df)

df %>% mutate(x %% y)
mf %>% mutate(x %% y)

## -----------------------------------------------------------------------------
translate_sql(mean(x), con = con)
translate_sql(mean(x, na.rm = TRUE), con = con)

## -----------------------------------------------------------------------------
translate_sql(mean(x, na.rm = TRUE), window = FALSE, con = con)

## -----------------------------------------------------------------------------
translate_sql(if (x > 5) "big" else "small", con = con)
translate_sql(switch(x, a = 1L, b = 2L, 3L), con = con)

## -----------------------------------------------------------------------------
translate_sql(foofify(x, y), con = con)

## -----------------------------------------------------------------------------
translate_sql(FOOFIFY(x, y), con = con)

## -----------------------------------------------------------------------------
translate_sql(x %LIKE% "%foo%", con = con)

## -----------------------------------------------------------------------------
translate_sql(x %||% y, con = con)

## -----------------------------------------------------------------------------
translate_sql(sql("x!"), con = con)
translate_sql(x == sql("ANY VALUES(1, 2, 3)"), con = con)

## -----------------------------------------------------------------------------
mf <- memdb_frame(x = 1, y = 2)

mf %>% 
  transmute(factorial = sql("x!")) %>% 
  show_query()

mf %>% 
  transmute(factorial = sql("CAST(x AS FLOAT)")) %>% 
  show_query()

## -----------------------------------------------------------------------------
try({
options(dplyr.strict_sql = TRUE)
translate_sql(glob(x, y), con = con)
})

## -----------------------------------------------------------------------------
knitr::include_graphics("windows.png", dpi = 300)

## -----------------------------------------------------------------------------
translate_sql(mean(G), con = con)
translate_sql(rank(G), con = con)
translate_sql(ntile(G, 2), con = con)
translate_sql(lag(G), con = con)

## -----------------------------------------------------------------------------
translate_sql(cummean(G), vars_order = "year", con = con)
translate_sql(rank(), vars_group = "ID", con = con)

## -----------------------------------------------------------------------------
# mutate(players,
#   min_rank(yearID),
#   order_by(yearID, cumsum(G)),
#   lead(G, order_by = yearID)
# )

