## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# # googledrive::
# drive_auth <- function(email = gargle::gargle_oauth_email(),
#                        path = NULL,
#                        scopes = "https://www.googleapis.com/auth/drive",
#                        cache = gargle::gargle_oauth_cache(),
#                        use_oob = gargle::gargle_oob_default(),
#                        token = NULL) {
#   # this catches a common error, where the user passes JSON for an OAuth client
#   # to the `path` argument, which only expects a service account token
#   gargle::check_is_service_account(path, hint = "drive_auth_configure")
# 
#   cred <- gargle::token_fetch(
#     scopes = scopes,
#     client = drive_oauth_client() %||% <BUILT_IN_DEFAULT_CLIENT>,
#     email = email,
#     path = path,
#     package = "googledrive",
#     cache = cache,
#     use_oob = use_oob,
#     token = token
#   )
#   if (!inherits(cred, "Token2.0")) {
#     # throw an informative error here
#   }
#   .auth$set_cred(cred)
#   .auth$set_auth_active(TRUE)
# 
#   invisible()
# }

## -----------------------------------------------------------------------------
# .auth <- NULL

## -----------------------------------------------------------------------------
.onLoad <- function(libname, pkgname) {
  utils::assignInMyNamespace(
    ".auth",
    gargle::init_AuthState(package = "googledrive", auth_active = TRUE)
  )
  
  # other stuff
}

## -----------------------------------------------------------------------------
# library(googledrive)
# 
# # first: download the OAuth client as a JSON file
# drive_auth_configure(
#   path = "/path/to/the/JSON/that/was/downloaded/from/gcp/console.json"
# )
# 
# drive_oauth_client()
# #> <gargle_oauth_client>
# #> name: acme-corp-google-client
# #> id: 123456789.apps.googleusercontent.com
# #> secret: <REDACTED>
# #> type: installed
# #> redirect_uris: http://localhost

## -----------------------------------------------------------------------------
# library(googledrive)
# 
# drive_auth_configure(api_key = "123456789")
# 
# drive_api_key()
# #> "123456789"

## -----------------------------------------------------------------------------
# # googledrive::
# drive_auth(email = "janedoe_work@gmail.com")

## -----------------------------------------------------------------------------
# # googledrive::
# drive_auth <- function(email = gargle::gargle_oauth_email(),
#                        path = NULL,
#                        scopes = "https://www.googleapis.com/auth/drive",
#                        cache = gargle::gargle_oauth_cache(),
#                        use_oob = gargle::gargle_oob_default(),
#                        token = NULL) { ... }

## -----------------------------------------------------------------------------
# # googledrive::
# drive_auth(scopes = "https://www.googleapis.com/auth/drive.readonly")

## -----------------------------------------------------------------------------
# # googledrive::
# request_generate <- function(endpoint = character(),
#                              params = list(),
#                              key = NULL,
#                              token = drive_token()) {
#   ept <- drive_endpoint(endpoint)
#   if (is.null(ept)) {
#     # throw error about unrecognized endpoint
#   }
# 
#   ## modifications specific to googledrive package
#   params$key <- key %||% params$key %||%
#     drive_api_key() %||% <BUILT_IN_DEFAULT_API_KEY>
#   if (!is.null(ept$parameters$supportsAllDrives)) {
#     params$supportsAllDrives <- TRUE
#   }
# 
#   req <- gargle::request_develop(endpoint = ept, params = params)
#   gargle::request_build(
#     path = req$path,
#     method = req$method,
#     params = req$params,
#     body = req$body,
#     token = token
#   )
# }

## -----------------------------------------------------------------------------
# # googledrive::
# drive_token <- function() {
#   if (isFALSE(.auth$auth_active)) {
#     return(NULL)
#   }
#   if (!drive_has_token()) {
#     drive_auth()
#   }
#   httr::config(token = .auth$cred)
# }

## -----------------------------------------------------------------------------
# # googledrive::
# drive_has_token <- function() {
#   inherits(.auth$cred, "Token2.0")
# }

## -----------------------------------------------------------------------------
# library(googledrive)
# library(googlesheets4)
# 
# drive_auth(email = "jane_doe@example.com") # gets a suitably scoped token
#                                            # and stashes for googledrive use
# 
# gs4_auth(token = drive_token())            # registers token with googlesheets4
# 
# # now work with both packages freely ...

## -----------------------------------------------------------------------------
# library(googledrive)
# 
# drive_auth(email = "janedoe_work@gmail.com")
# # do stuff with Google Drive here, with Jane Doe's "work" account
# 
# drive_auth(email = "janedoe_personal@gmail.com")
# # do other stuff with Google Drive here, with Jane Doe's "personal" account
# 
# drive_auth(path = "/path/to/a/service-account.json")
# # do other stuff with Google Drive here, using a service account

