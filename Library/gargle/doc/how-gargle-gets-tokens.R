## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(gargle)

## -----------------------------------------------------------------------------
# token_fetch(scopes, ...)

## -----------------------------------------------------------------------------
writeLines(names(cred_funs_list()))

## -----------------------------------------------------------------------------
# token_fetch(token = <TOKEN2.0>)
# 
# credentials_byo_oauth2(
#   token = <TOKEN2.0>
# )

## -----------------------------------------------------------------------------
# token_fetch(scopes = <SCOPES>, path = "/path/to/your/service-account.json")
# 
# # credentials_byo_oauth2() fails because no `token`,
# # which leads to this call:
# credentials_service_account(
#   scopes = <SCOPES>,
#   path = "/path/to/your/service-account.json"
# )

## -----------------------------------------------------------------------------
# token_fetch(scopes = <SCOPES>, path = "/path/to/your/external-account.json")
# 
# # credentials_byo_oauth2() fails because no `token`,
# # credentials_service_account() fails because the JSON provided via
# #   `path` is not of type "service_account",
# # which leads to this call:
# credentials_external_account(
#   scopes = <SCOPES>,
#   path = "/path/to/your/external-account.json"
# )

## -----------------------------------------------------------------------------
# token_fetch(scopes = <SCOPES>)
# 
# # credentials_byo_oauth2() fails because no `token`,
# # credentials_service_account() fails because no `path`,
# # credentials_external_account() fails because no `path`,
# # which leads to this call:
# credentials_app_default(
#   scopes = <SCOPES>
# )

## -----------------------------------------------------------------------------
# ${GOOGLE_APPLICATION_CREDENTIALS}
# ${CLOUDSDK_CONFIG}/application_default_credentials.json
# 
# # on Windows:
# %APPDATA%\gcloud\application_default_credentials.json
# %SystemDrive%\gcloud\application_default_credentials.json
# C:\gcloud\application_default_credentials.json
# 
# # on not-Windows:
# ~/.config/gcloud/application_default_credentials.json

## -----------------------------------------------------------------------------
# token_fetch(scopes = <SCOPES>)
# # or perhaps
# token_fetch(scopes = <SCOPES>, service_account = <SERVICE_ACCOUNT>)
# 
# # credentials_byo_oauth2() fails because no `token`,
# # credentials_service_account() fails because no `path`,
# # credentials_external_account() fails because no `path`,
# # credentials_app_default() fails because no ADC found,
# # which leads to one of these calls:
# credentials_gce(
#   scopes = <SCOPES>,
#   service_account = "default"
# )
# # or
# credentials_gce(
#   scopes = <SCOPES>,
#   service_account = <SERVICE_ACCOUNT>
# )

## -----------------------------------------------------------------------------
# token_fetch(scopes = <SCOPES>)
# 
# # credentials_byo_oauth2() fails because no `token`,
# # credentials_service_account() fails because no `path`,
# # credentials_external_account() fails because no `path`,
# # credentials_app_default() fails because no ADC found,
# # credentials_gce() fails because not on GCE,
# # which leads to this call:
# credentials_user_oauth2(
#   scopes = <SCOPES>,
#   app = <OAUTH_APP>,
#   package = "<PACKAGE>"
# )

## -----------------------------------------------------------------------------
# # user initiates auth or does something that triggers it indirectly
# THINGY_auth()
# 
# # which then calls
# gargle::token_fetch(
#   scopes  = <SCOPES_NEEDED_FOR_THE_THINGY_API>,
#   app     = thingy_app(),
#   package = "thingyr"
# )
# 
# # which leads to this call:
# credentials_user_oauth2(
#   scopes  = <SCOPES_NEEDED_FOR_THE_THINGY_API>,
#   app     = thingy_app(),
#   package = "thingyr"
# )

## -----------------------------------------------------------------------------
# gargle2.0_token(
#   email   = gargle_oauth_email(),
#   app     = thingy_app(),
#   package = "thingyr",
#   scope   = <SCOPES_NEEDED_FOR_THE_THINGY_API>,
#   cache   = gargle_oauth_cache()
# )

## -----------------------------------------------------------------------------
# The thingyr package is requesting access to your Google account.
# Enter '1' to start a new auth process or select a pre-authorized account.
# 1: Send me to the browser for a new auth process.
# 2: janedoe_personal@gmail.com
# 3: janedoe@example.com
# 4: janedoe_work@gmail.com
# Selection:

## -----------------------------------------------------------------------------
# thingy_auth(email = "janedoe_work@gmail.com")

## -----------------------------------------------------------------------------
# gargle_oauth_sitrep()
# #> 14 tokens found in this gargle OAuth cache:
# #> ~/Library/Caches/gargle
# #'
# #' email                         app         scopes                         hash...
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
writeLines(names(cred_funs_list()))

## -----------------------------------------------------------------------------
# gargle::cred_funs_add(credentials_gce = NULL)

