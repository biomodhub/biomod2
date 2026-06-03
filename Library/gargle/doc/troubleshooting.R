## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(gargle)

## -----------------------------------------------------------------------------
gargle_verbosity()

## -----------------------------------------------------------------------------
# save current value
op <- options(gargle_verbosity = "debug")

gargle_verbosity()

# restore original value
options(op)

## -----------------------------------------------------------------------------
gargle_verbosity()

with_gargle_verbosity(
  "debug",
  gargle_verbosity()
)

gargle_verbosity()

f <- function() {
  local_gargle_verbosity("debug")
  gargle_verbosity()
}

f()

gargle_verbosity()

## -----------------------------------------------------------------------------
# gargle_oauth_sitrep()
# #' > 14 tokens found in this gargle OAuth cache:
# #' '~/Library/Caches/gargle'
# #'
# #' email                         app         scope                          hash...
# #' ----------------------------- ----------- ------------------------------ ----------
# #' abcdefghijklm@gmail.com       thingy      ...bigquery, ...cloud-platform 128f9cc...
# #' buzzy@example.org             gargle-demo                                15acf95...
# #' stella@example.org            gargle-demo ...drive                       4281945...
# #' abcdefghijklm@gmail.com       gargle-demo ...drive                       48e7e76...
# #' abcdefghijklm@gmail.com       tidyverse                                  69a7353...
# #' nopqr@ABCDEFG.com             tidyverse   ...spreadsheets.readonly       86a70b9...
# #' abcdefghijklm@gmail.com       tidyverse   ...drive                       d9443db...
# #' nopqr@HIJKLMN.com             tidyverse   ...drive                       d9443db...
# #' nopqr@ABCDEFG.com             tidyverse   ...drive                       d9443db...
# #' stuvwzyzabcd@gmail.com        tidyverse   ...drive                       d9443db...
# #' efghijklmnopqrtsuvw@gmail.com tidyverse   ...drive                       d9443db...
# #' abcdefghijklm@gmail.com       tidyverse   ...drive.readonly              ecd11fa...
# #' abcdefghijklm@gmail.com       tidyverse   ...bigquery, ...cloud-platform ece63f4...
# #' nopqr@ABCDEFG.com             tidyverse   ...spreadsheets                f178dd8...

## -----------------------------------------------------------------------------
# install.packages("googlesheets4")

## -----------------------------------------------------------------------------
knitr::include_graphics("deleted_client.png")

## -----------------------------------------------------------------------------
# only run the chunk below in settings that are known to be safe, i.e. where
# occasional, incidental failure is OK
can_decrypt <- gargle::secret_has_key("GARGLE_KEY")

## -----------------------------------------------------------------------------
library(gargle)

req <- request_build(
  method = "GET",
  path = "webfonts/v1/webfonts",
  params = list(
    sort = "popularity"
  ),
  key = gargle_api_key(),
  base_url = "https://www.googleapis.com"
)
resp <- request_make(req)
out <- response_process(resp)

lr <- gargle_last_response()
tmp <- tempfile("gargle-last-response-")
saveRDS(lr, tmp)
# you could share this .rds file with a colleague or the gargle maintainer

# how it would look to complete the round trip, i.e. load this on the other end
rt_lr <- readRDS(tmp)

all.equal(lr, rt_lr)

# clean up
unlink("tmp")

