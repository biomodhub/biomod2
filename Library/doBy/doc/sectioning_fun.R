## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options("digits"=3)
library(doBy)
library(microbenchmark)
#devtools::load_all()

## -----------------------------------------------------------------------------
fun  <- function(a, b, c=4, d=9) {
    a + b + c + d
}

## -----------------------------------------------------------------------------
fun_def <- section_fun(fun, list(b=7, d=10))
fun_def
fun_body <- section_fun(fun, list(b=7, d=10), method="sub")
fun_body
fun_env <- section_fun(fun, list(b=7, d=10), method = "env")
fun_env

## -----------------------------------------------------------------------------
get_section(fun_env) 
## same as: attr(fun_env, "arg_env")$args 
get_fun(fun_env) 
## same as: environment(fun_env)$fun

## -----------------------------------------------------------------------------
fun(a=10, b=7, c=5, d=10)
fun_def(a=10, c=5)
fun_body(a=10, c=5)
fun_env(a=10, c=5)

## -----------------------------------------------------------------------------
inv_toep <- function(n) {
    solve(toeplitz(1:n))
}

## ----eval=F-------------------------------------------------------------------
# microbenchmark(
#     inv_toep(4), inv_toep(8), inv_toep(16),
#     times=3
# )

## -----------------------------------------------------------------------------
n.vec  <- c(4, 8, 16)
fun_list <- lapply(n.vec,
                   function(ni) {
                       section_fun(inv_toep, list(n=ni))
                   })
fun_list

## -----------------------------------------------------------------------------
fun_list[[1]]
fun_list[[1]]()

## -----------------------------------------------------------------------------
bquote_list <- function(fun_list) {
    lapply(fun_list, function(g){
        bquote(.(g)())
    })
}

## -----------------------------------------------------------------------------
bq_fun_list <- bquote_list(fun_list)
bq_fun_list
bq_fun_list[[1]]
eval(bq_fun_list[[1]])

## -----------------------------------------------------------------------------
microbenchmark(
  list = bq_fun_list,
  times = 5
)

## -----------------------------------------------------------------------------
n.vec  <- seq(20, 80, by=20)
fun_def <- lapply(n.vec,
                  function(n){
                      section_fun(inv_toep, list(n=n), method="def")
                  })
fun_body <- lapply(n.vec,
                  function(n){
                      section_fun(inv_toep, list(n=n), method="sub")
                  })
fun_env <- lapply(n.vec,
                  function(n){
                      section_fun(inv_toep, list(n=n), method="env")
                  })

names(fun_def)  <- paste0("def", n.vec)
names(fun_body) <- paste0("body", n.vec)
names(fun_env)  <- paste0("env", n.vec)

bq_fun_list <- bquote_list(c(fun_def, fun_body, fun_env))
bq_fun_list |> head()

mb <- microbenchmark(
  list  = bq_fun_list,
  times = 2
)
mb

