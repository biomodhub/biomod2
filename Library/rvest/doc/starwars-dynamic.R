## ----eval = FALSE, echo = FALSE-----------------------------------------------
# library(magrittr)
# crawl_html <- function(x) {
#   x |>
#     gsub("\r", "", .) |>
#     gsub("\n\n", "</p><p>", .) |>
#     gsub("\n", " ", .) |>
#     paste0("<p>", ., "</p>")
# }
# 
# fields <- c("episode_id", "title", "release_date", "director", "opening_crawl")
# repurrrsive::sw_films |>
#   lapply(\(film) film[fields]) |>
#   lapply(function(film) {
#     film$opening_crawl <- crawl_html(film$opening_crawl)
#     film
#   }) |>
#   jsonlite::write_json("vignettes/starwars.json", pretty = TRUE, auto_unbox = TRUE)

