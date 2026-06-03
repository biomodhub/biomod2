## ----message=FALSE------------------------------------------------------------
library(htmlTable)
library(dplyr)
library(magrittr)
data("mtcars")

mtcars %<>%
  mutate(am = factor(am, levels = 0:1, labels = c("Automatic", "Manual")),
         vs = factor(vs, levels = 0:1, labels = c("V-shaped", "straight")))

mtcars %>% 
  head(3) %>% 
  select(Transmission = am, Gas = mpg, Weight = wt) %>% 
  htmlTable()

## -----------------------------------------------------------------------------
mtcars %>% 
  head(3) %>% 
  select(Transmission = am, Gas = mpg, Weight = wt) %>% 
  txtRound(digits = 1) %>% 
  htmlTable()

## -----------------------------------------------------------------------------
txtRound(c(1, 1.1034), digits = 2)

# Use a character to convert
txtRound("1.2333", digits = 2)

## -----------------------------------------------------------------------------
# Large numbers can be combined with the txtInt option
txtRound(12345.12, digits = 1, txtInt_args = TRUE)

txtRound(12345.12, digits = 1, txtInt_args = list(language = "se", html = FALSE))

## -----------------------------------------------------------------------------
mtcars %>% 
  head(3) %>% 
  select(mpg, wt) %>% 
  txtRound(mpg, wt_txt = wt, digits = 1)

## -----------------------------------------------------------------------------
mtcars %>% 
  head(3) %>% 
  select(mpg, qsec, wt) %>% 
  txtRound(digits = list(wt = 2, .default = 1))

## -----------------------------------------------------------------------------
mtcars_matrix <- mtcars %>% 
  select(mpg, qsec, wt) %>% 
  head(3) %>% 
  as.matrix()

mtcars_matrix %>% 
  txtRound(digits = 1)

## -----------------------------------------------------------------------------
mtcars_matrix %>% 
  txtRound(excl.cols = "^wt$",
           excl.rows = "^Mazda RX4$",
           digits = 1)

## -----------------------------------------------------------------------------
mtcars_matrix %>% 
  txtRound(digits = list(mpg = 0, wt = 2, .default = 1))

## -----------------------------------------------------------------------------
txtInt(1e7)

## -----------------------------------------------------------------------------
txtInt(1e7, language = "SI", html = FALSE)

txtInt(1e7, language = "SI", html = TRUE)

## -----------------------------------------------------------------------------
txtPval(c(0.1233213, 0.035, 0.001, 0.000001), html = FALSE)

# The < sign is less-than in html code '&lt;'
txtPval(c(0.05, 0.001, 0.000001), html = TRUE)

## -----------------------------------------------------------------------------
txtMergeLines("Line 1",
              "Line 2",
              "Line 3")

## -----------------------------------------------------------------------------
txtMergeLines("Line 1
               Line 2
               Line 3")

## -----------------------------------------------------------------------------
txtMergeLines("Line 1
               Line 2
               Line 3",
              html = FALSE)

