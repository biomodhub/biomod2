## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(gargle)

## -----------------------------------------------------------------------------
(path_to_installed_client <- system.file(
  "extdata", "client_secret_installed.googleusercontent.com.json",
  package = "gargle"
))
jsonlite::prettify(scan(path_to_installed_client, what = character()))
(client <- gargle_oauth_client_from_json(path_to_installed_client))
class(client)

(path_to_web_client <- system.file(
  "extdata", "client_secret_web.googleusercontent.com.json",
  package = "gargle"
))
jsonlite::prettify(scan(path_to_web_client, what = character()))
(client <- gargle_oauth_client_from_json(path_to_web_client))
class(client)

## -----------------------------------------------------------------------------
# # BEFORE
# drive_auth_configure <- function(app, path, api_key) {
#   # not showing this code
#   .auth$set_app(app)
#   # more code we're not showing
# }
# 
# drive_oauth_app <- function() .auth$app
# 
# # AFTER
# drive_auth_configure <- function(client, path, api_key, app = deprecated()) {
#   if (lifecycle::is_present(app)) {
#     lifecycle::deprecate_warn(
#       "2.1.0",
#       "drive_auth_configure(app)",
#       "drive_auth_configure(client)"
#     )
#     drive_auth_configure(client = app, path = path, api_key = api_key)
#   }
# 
#   # not showing this code
#   .auth$set_client(client)
#   # more code we're not showing
# }
# 
# drive_oauth_client <- function() .auth$client
# 
# drive_oauth_app <- function() {
#   lifecycle::deprecate_warn(
#     "2.1.0", "drive_oauth_app()", "drive_oauth_client()"
#   )
#   drive_oauth_client()
# }

## -----------------------------------------------------------------------------
# usethis::use_lifecycle()

## -----------------------------------------------------------------------------
# drive_auth <- function(...) {
#   # code not shown
#   cred <- gargle::token_fetch(
#     scopes = scopes,
#     # app = drive_oauth_client() %||% <BUILT_IN_DEFAULT_CLIENT>,   # BEFORE
#     client = drive_oauth_client() %||% <BUILT_IN_DEFAULT_CLIENT>,  # AFTER
#     email = email,
#     path = path,
#     package = "googledrive",
#     cache = cache,
#     use_oob = use_oob,
#     token = token
#   )
#   # code not shown
# }

