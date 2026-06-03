## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# library(googlesheets4)
# 
# gs4_auth_configure(api_key = "YOUR_API_KEY_GOES_HERE")
# gs4_deauth()
# 
# # now you can read public resources, such as official example Sheets,
# # without any need for auth
# gs4_example("gapminder") |>
#   read_sheet()

## -----------------------------------------------------------------------------
# library(googledrive)
# 
# google_client <- gargle::gargle_oauth_client_from_json(
#   path = "/path/to/the/JSON/that/was/downloaded/from/gcp/console.json",
#   name = "acme-corp-google-client"
# )
# drive_auth_configure(client = google_client)
# 
# # now any new OAuth tokens are obtained with the configured client

## -----------------------------------------------------------------------------
# # googledrive
# drive_auth(path = "/path/to/your/service-account-token.json")

