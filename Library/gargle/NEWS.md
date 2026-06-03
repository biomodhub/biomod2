# gargle 1.6.1

* `gargle_last_response()` and `gargle_last_content()` are newly exported.

* `response_process(error_class =)` is a new optional argument to add a subclass to any error thrown while processing a response.

# gargle 1.6.0

* When retrying a request, the messaging reveals more detail about the failed
  request and counts down until the retry happens (part of #64).
* Fixed some bugs around surfacing the details of a request failure:
  - Better handling when the error details include a localized message (#293)
  - Defensive escaping of `{..}` in Google-provided error messages, to prevent
    `cli::cli_abort()` from trying (and failing) to do string interpolation
    (https://github.com/tidyverse/googlesheets4/issues/319)
* gargle is better able to detect when it's running on Posit Workbench, but not
  necessarily in RStudio, such as in Positron or VS Code (#291).

# gargle 1.5.2

* Fixed a bug in an internal helper that validates input specifying a service
  account. The helper targets a common mistake where the JSON for an OAuth
  client is provided to an argument that is meant for a service account (#270).

# gargle 1.5.1

* Completed some overlooked, unfinished work around the OAuth "app" to "client"
  transition that affected out-of-band auth (#263, #264).

* The `secret_*()` functions are more discoverable via documentation.

# gargle 1.5.0

* gargle's existing unexported `secret_*()` functions are deprecated, in favor
  of new, exported `secret_*()` functions that are built on or inlined from
  httr2. The `vignette("managing-tokens-securely")` is updated to reflect the
  new, recommended strategy for encrypting secrets.
  - `secret_encrypt_json()` / `secret_decrypt_json()` are new gargle-specific
    functions.
  - `secret_write_rds()` / `secret_read_rds()`, `secret_make_key()`, and
    `secret_had_key()` are basically copies of their httr2 counterparts.
  - Legacy functions to move away from: `secret_pw_name()`, `secret_pw_gen()`,
    `secret_pw_exists()`, `secret_pw_get()`, `secret_can_decrypt()`,
    `secret_read()`, `secret_write()`.
  - The new approach makes it much easier to use gargle functions to encrypt and
    decrypt credentials in a project that is *not* necessarily an R package.

* The transition from OAuth "app" to OAuth "client" is fully enacted now. This
  process tarted in v1.3.0, when the `"gargle_oauth_client"` class was
  introduced, to support the new pseudo-OOB auth flow. The deprecations are
  implemented to preserve backwards compatibility for some time. In this
  release, function, argument, and field names are all updated to the "client"
  terminology:

  - `init_AuthState(client =)` instead of `init_AuthState(app =)`
  - `AuthState$client` instead of `AuthState$app`
  - `AuthState$set_client()` instead of `AuthState$set_app()`
  - `gargle2.0_token(client =)` instead of `gargle2.0_token(app =)`
  - `credentials_user_oauth2(client =)` instead of
    `credentials_user_oauth2(app =)`

  A new `vignette("oauth-client-not-app")` explains how a wrapper package should
  adapt.

* When the `"gargle_verbosity"` option is set to `"debug"`, there are more debugging messages around user credentials. Specifically, more information is available on the email, OAuth client, and scopes, with the goal of better understanding why a cached token is (or is not) being used.

* `check_is_service_account()` is a new function for use in wrapper packages to throw a more informative error when a user provides JSON for an OAuth client to an argument that is expecting JSON for a service account.

* `response_process()` has improved handling of responses that represent an HTTP error with HTML content (as opposed to the expected and preferred JSON) (#254).

* `response_process(call = caller_env())` is a new argument that is passed along to various helpers, which can improve error reporting for user-facing functions that call `response_process()` (#255).

# gargle 1.4.0

## Google Compute Engine

* `credentials_gce(scopes = NULL)` is now equivalent to `credentials_gce(scopes = "https://www.googleapis.com/auth/cloud-platform")`, i.e. there's an even stronger current towards the recommended "cloud-platform" scope.

* `credentials_gce(scopes =)` now includes those `scopes` in its request to the metadata server for an access token (#216). Note that the scopes for a GCE access token are generally pre-determined for the instance and its associated service account at creation/launch time and these requested `scopes` will have no effect. But this seems to do no harm and it is possible that there are contexts where this is useful.

* `credentials_gce()` now emits considerably more information when the `"gargle_verbosity"` option is set to `"debug"`. For example, it reports mismatches between requested scopes and instance scopes and between requested scopes and the access token's actual scopes.

* `credentials_gce()` stores the actual scopes of the received access token, which can differ from the requested scopes. This is also noted when the `"gargle_verbosity"` option is set to `"debug"`.

* The `GceToken` R6 class gains a better `$print()` method that is more similar to gargle's treatment of tokens obtained with other flows.

## Behaviour in a cloud/server context

* gargle is better able to detect when it's running on Posit Workbench or RStudio Server, e.g., in a subprocess.

* `gargle_oauth_client_type()` is a new function that returns either "installed"
or "web".
It returns the value of the new global option by the same name (`"gargle_oauth_client_type"`), if defined.
If the option is not defined, returns "web" on RStudio Server, Posit Workbench, Posit Cloud, or Google Colaboratory and "installed" otherwise.
In the context of out-of-band (OOB) auth, an "installed" client type leads to the conventional OOB flow (only available for GCP projects in testing mode) and a "web" client leads to the new pseudo-OOB flow.
The option and accessor have been added to cover contexts other than those mentioned above where it is helpful to request a "web" client.

* `credentials_user_oauth2()` now works in Google Colaboratory (#140).

## Everything else

* gargle now elicits user input via `readline()`, instead of via `utils::menu()`, which is favorable for interacting with the user in a Jupyter notebook (#242).

* The roxygen templating functions that wrapper packages can use to generate standardized documentation around auth have been updated to reflect gargle's pivot from OAuth "app" to "client".
Changes of note:

  - `PREFIX_auth_configure_description()` crosslinks to `PREFIX_oauth_client()`
    now, not `PREFIX_oauth_app()`. So this assumes the package has indeed
    introduced the `PREFIX_oauth_client()` function (and, presumably, has
    deprecated `PREFIX_oauth_app()`).
  - `PREFIX_auth_configure_params()` gains `client` argument.
  - `PREFIX_auth_configure_params()` deprecates the `app` argument and uses a
    lifecycle badge. It is assumed that the badge SVG is present, which can be
    achieved with `usethis::use_lifecycle()`.
  - `PREFIX_auth_configure_params()` crosslinks to
    `gargle::gargle_oauth_client_from_json()`. The wrapper package therefore
    needs to state a minimum version for gargle, e.g. `gargle (>= 1.3.0)` (or
    higher).

* `credentials_byo_oauth2()` works now for (variations of) service account tokens, as intended, not just for user tokens (#250). It also emits more information about scopes when the `"gargle_verbosity"` option is set to `"debug"`.

# gargle 1.3.0

## (Partial) deprecation out-of-band (OOB) auth flow

On February 16, 2022, Google announced the gradual deprecation of the out-of-band (OOB) OAuth flow.
OOB **still works** if the OAuth client is associated with a GCP project that is in testing mode and this is not going away.
But OOB is no longer supported for projects in production mode.
To be more accurate, some production-mode projects have gotten an extension to permit the use of OOB auth for a bit longer, but that's just a temporary reprieve.

The typical user who will (eventually) be impacted is:

* Using R via RStudio Server, Posit Workbench, or Posit Cloud.
* Using tidyverse packages such as googledrive, googlesheets4, or bigrquery.
* Relying on the built-in OAuth client. Importantly, this client is associated
  with a GCP project that is in production mode.

The phased deprecation of OOB is nearly complete and we expect conventional OOB to stop working with the built-in tidyverse OAuth client on February 1, 2023, at the latest.

**In preparation for this, gargle has gained support for a new flow, which we call pseudo-OOB (in contrast to conventional OOB)**.
The pseudo-OOB flow is triggered when `use_oob = TRUE` (an existing convention in gargle and gargle-using packages) and the configured OAuth client is of "Web application" type.
The gargle/googledrive/googlesheets4/bigrquery packages should now default to a "Web application" client on RStudio Server, Posit Workbench and Posit Cloud, leading the user through the pseudo-OOB flow.
Other than needing to re-auth once, affected users should still find that things "just work".

Read the `vignette("auth-from-web")` for more.

## gargle-specific notion of OAuth client

`gargle_oauth_client()` is a new constructor for an S3 class by the same name.
There are two motivations:

  - To adjust to Google's deprecation of conventional OOB and to support
    gargle's new pseudo-OOB flow, it is helpful for gargle to know whether an
    OAuth client ID is of type "Web application" or "Desktop app". That means we
    need a Google- and gargle-specific notion of an OAuth client, so we can
    introduce a `type` field.
  - A transition from httr to httr2 is on the horizon, so it makes sense to
    look more toward `httr2:oauth_client()` than to `httr::oauth_app()`.
    gargle's vocabulary is generally shifting towards "client" and away from
    "app".

`oauth_app_from_json()` has therefore been (soft) deprecated, in favor of a new function `gargle_oauth_client_from_json()`, which is the preferred way to instantiate an OAuth client, since the downloaded JSON conveys the client type and redirect URI(s).
As a bridging measure, `gargle_oauth_client` currently inherits from httr's `oauth_app`, but this probably won't be true in the long-term.

`gargle_client(type =)` replaces `gargle_app()`.

## Google Compute Engine and Google Kubernetes Engine

`credentials_gce()` no longer asks the user about initiating an OAuth cache, which is not relevant to that flow (#221).

`gce_instance_service_accounts()` is a newly exported utility that exposes the service accounts available from the metadata server for the current instance (#234).

The global option `"gargle.gce.timeout"` is newly documented in `credentials_gce()`.
This controls the timeout, in seconds, for requests to the metadata server.
The default value (or strategy) for setting this should often suffice, but the option exists for those with an empirical need to increase the timeout (#186, #195).

`vignette("non-interactive-auth")` has a new section "Workload Identity on Google Kubernetes Engine (GKE)" that explains how gargle supports the use of workload identity for applications running on GKE. This is the recommended method of auth in R code running on GKE that needs to access other Google Cloud services, such as the BigQuery API (#197, #223, @MarkEdmondson1234).

## Credential function registry

It's gotten a bit easier to work with the credential registry.
The primary motivation is that, for example, on Google Compute Engine, you might
actually want to suppress auth with the default service account and auth as a
normal user instead.
This is especially likely to come up with gmailr / the Gmail API.

* The credential-fetcher `credentials_byo_oauth2()` has been moved to the very
  beginning of the default registry. The logic is that a user who has specified
  a non-`NULL` value of `token` must mean business and does not want automagic
  auth methods like ADC or GCE to be tried before using their `token`
  (#187, #225).

* The `...` in `cred_funs_all()` are now
  [dynamic dots](https://rlang.r-lib.org/reference/dyn-dots.html) (#224).

* Every registered credential function must have a unique name now.
  This is newly enforced by `cred_funs_add()` and `cred_funs_set()` (#224).

* `cred_funs_list_default()` is a new function that returns gargle's default
  list of credential functions (#226).

* `cred_funs_add(cred_fun = NULL)` is now available to remove a credential
  function from the registry (#224).

* `with_cred_funs()` and `local_cred_funs()` are new helpers for making narrowly
  scoped changes to the registry (#226).

* The `ls` argument of `cred_funs_set()` has been renamed to `funs` (#226).

* In general, credential registry functions now return the current registry,
  invisibly (#224).

# gargle 1.2.1

* Help files below `man/` have been re-generated, so that they give rise to valid HTML5. (This is the impetus for this release, to keep the package safely on CRAN.)

* We have switched to newer oauth2.googleapis.com-based OAuth2 URIs, moving away from the accounts.google.com and googleapis.com/oauth2 equivalents.

* `credentials_gce()` no longer validates the requested scopes against instance scopes.
In practice, it's easy for this check to be more of a nuisance than a help (#161, #185 @craigcitro).

* `request_retry()` retries for an expanded set of HTTP codes: 408, 429, 500, 502, 503. Previously, retries were limited to 429 (#169).

## Dependency changes

* The minimum versions of rlang and testthat have been bumped. The motivation is to exploit and adapt to the changes to the display of error messages.

# gargle 1.2.0

## Workload identity federation

`credentials_external_account()` is a new function that implements "workload identity federation", a new (as of April 2021) keyless authentication mechanism.
This allows applications running on a non-Google Cloud platform, such as AWS, to access Google Cloud resources without using a conventional service account token, eliminating the security problem posed by long-lived, powerful service account credential files.

`credentials_external_account()` has been inserted into the default registry of credential-fetchers tried by `token_fetch()`, which makes it automatically available in certain wrapper packages, such as bigrquery.
`credentials_app_default()` recognizes the JSON configuration for an external account and passes such a call along to `credentials_external_account()`.

This new feature is still experimental and currently only supports AWS.
This [blog post](https://cloud.google.com/blog/products/identity-security/enable-keyless-access-to-gcp-with-workload-identity-federation) provides a good high-level introduction to workload identity federation.

## Other changes

The `email` argument of `credentials_user_oauth2()` accepts domain-only email specification via a glob pattern.
The goal is to make it possible for code like `PKG_auth(email = "*@example.com")` to identify a suitable cached token when executed on the machine of either `alice@example.com` or `bob@example.com`.

gargle now throws errors via `cli::cli_abort()`, which means error messages now have the same styling as informational messages.

## Dependency changes

aws.ec2metadata and aws.signature are new in Suggests.

# gargle 1.1.0

## OAuth token cache

Two changes affect stored user OAuth tokens:

* The default cache location has moved, to better align with general
  conventions around where to cache user data. Here's how that looks for a
  typical user:
  - Typical before, macOS: `~/.R/gargle/gargle-oauth`
  - Typical after, macOS: `~/Library/Caches/gargle`
  - Typical before, Windows: `C:/Users/jane/.R/gargle/gargle-oauth`
  - Typical after, Windows: `C:/Users/jane/AppData/Local/gargle/gargle/Cache`
* Tokens created with one of the built-in OAuth apps provided by the tidyverse
  packages are checked for validity. Tokens made with an old app are deleted.
  Note that we introduced a new OAuth app in gargle v1.0.0 and the previous
  app could be disabled at any time.
  - Nickname of previous tidyverse OAuth app: `tidyverse-calliope`
  - Nickname of tidyverse OAuth app as of gargle v1.0.0: `tidyverse-clio`

For users who accept all default behaviour around OAuth, these changes just mean you will see some messages about cleaning and moving the token cache.
These users can also expect to go through interactive auth (approximately once per package / API), to obtain fresh tokens made with the current tidyverse OAuth app.

If the rolling of the tidyverse OAuth app is highly disruptive to your workflow, this is a good wake-up call that you should be using your own OAuth app or, perhaps, an entirely different auth method, such as using a service account token in non-interactive settings.
As always, these articles explain how to take more control of auth:
 * <https://gargle.r-lib.org/articles/get-api-credentials.html>
 * <https://gargle.r-lib.org/articles/non-interactive-auth.html>

## User interface

The user interface has gotten more stylish, thanks to the cli package (<https://cli.r-lib.org>).

All errors thrown by gargle route through `rlang::abort()`, providing better access to the backtrace and, potentially, error data.
These errors have, at the very least, the `gargle_error` class and may also have additional subclasses.

`gargle_verbosity()` replaces `gargle_quiet()`.
Each such function is (or was) a convenience wrapper to query the option with that name.
Therefore, the option named "gargle_verbosity" now replaces "gargle_quiet".
If "gargle_verbosity" is unset, the old "gargle_quiet" is still consulted, but the user is advised to update their usage.

The new "gargle_verbosity" option is more expressive and has three levels:

* "debug", equivalent to the previous `gargle_quiet = FALSE`. Use for debugging
  and troubleshooting.
* "info" (the default), basically equivalent to the previous
  `gargle_quiet = TRUE`. Since gargle is not a user-facing package, it has very
  little to say and only emits messages that end users really need to see.
* "silent", no previous equivalent and of little practical significance. But it
  can be used to suppress all gargle messages.

The helpers `with_gargle_verbosity()` and `local_gargle_verbosity()` make it easy to temporarily modify the verbosity level, in the spirit of the [withr package](https://withr.r-lib.org).

## Other changes

There is special error handling when OAuth token refresh fails, due to deletion of the associated OAuth app.
This should help users who are relying on the default app provided by a package and, presumably, they need to update that package (#168).

`gargle_oob_default()` returns `TRUE` unconditionally when running in RStudio Server.

`response_process()` gains a `remember` argument.
When `TRUE` (the default), gargle stores the most recent response internally (with auth tokens redacted).
Unexported functions `gargle:::gargle_last_response()` and `gargle:::gargle_last_content()` facilitate *post mortem* analysis of, e.g., a failed request (#152).

`google.rpc.ErrorInfo` errors are explicitly handled now, resulting in a more informative error message.

`request_retry()` is better able to detect when the per-user quota has been exhausted (vs. the per-project quota), resulting in a more informed choice of backoff.

## Dependency changes

cli is new in Imports.

rstudioapi is new in Imports.

rappdirs is new in Imports.

httpuv is new in Suggests. We encourage its installation in interactive sessions, if we're about to initiate OAuth flow, unless it's clear that out-of-band auth is inevitable.

gargle now relies on testthat >= 3.0.0 and, specifically, uses third edition features.

mockr is new in Suggests, since `testthat::use_mock()` is superseded.

# gargle 1.0.0

* Better handling of `BadRequest` errors, i.e. more specifics are revealed.

* `oauth_app_from_json` now supports JSON files from the "Web application"
  client type (#155).

* `request_retry()` is a drop-in substitute for `request_make()` that uses (modified) exponential backoff to retry requests that fail with error `429 RESOURCE_EXHAUSTED` (#63).

* Credentials used in selected client packages have been rolled. Users of bigrquery, googledrive, and googlesheets4 can expect a prompt to re-authorize the "Tidyverse API Packages" when using an OAuth user token. This has no impact on users who use their own OAuth app (i.e. client ID and secret) or those who use service account tokens.

# gargle 0.5.0

* [Troubleshooting gargle auth](https://gargle.r-lib.org/articles/troubleshooting.html)
  is a new vignette.

* All user-facing messaging routes through `rlang::inform()`, which (as of
  rlang 0.4.2) prints to standard output in interactive sessions and to
  standard error in non-interactive sessions (#133). Messaging remains under
  the control of the `"gargle_quiet"` option, which defaults to `TRUE`.

* The `Gargle2.0` class gains its own `$refresh()` method, which removes a
  token from gargle's cache when it cannot be refreshed (#79).

* `credentials_service_account()` and `credentials_app_default()` gain an
  optional `subject` argument, which can be used to pass a subject claim along
  to `httr::oauth_service_token()` (#131, @samterfa).

* `request_make()` defaults to `encode = "json"`, which differs from the httr
  default, but aligns better with Google APIs (#124).

* `field_mask()` is a utility function for constructing a
  Protocol-Buffers-style, JSON-encoded field mask from a named R list.

* All R6 classes use the new documentation capabilities that appeared in
  roxygen2 7.0.0.

* OAuth2 flow can only be initiated when `rlang::is_interactive()` is `TRUE`. If
  a new token is needed in a non-interactive session, gargle now throws an
  error (#113).

* The application default credentials path is fixed on non-Windows platforms
  (#115, @acroz).

* `request_develop()` can accept a parameter that appears in both the path and
  the body (#123).

* `response_process()` explicitly declares the UTF-8 encoding of the content in
  Google API responses [tidyverse/googlesheets4#26](https://github.com/tidyverse/googlesheets4/issues/26).

* `response_process()` is able to expose details for a wider set of errors.

# gargle 0.4.0

* Eliminated uninformative failure when OAuth tokens cached on R <= 3.5 are re-loaded on R >= 3.6. The change to the default serialization version (2 vs. 3) creates an apparent mismatch between a token's hash and its key. Instead of inexplicably failing, now we attempt to repair the cache and carry on (#109, [tidyverse/googledrive#274](https://github.com/tidyverse/googledrive/issues/274).

* In a non-interactive context, gargle will use a cached OAuth token, if it discovers (at least) one, even if the user has not given explicit instructions. We emit a recommendation that the user make their intent unambiguous and link to the vignette on non-interactive auth (#92).

* gargle consults the option `"httr_oob_default"`, if the option `"gargle_oob_default"` is unset. This is part of an effort to automatically detect the need for out-of-band auth in more situations (#102).

* `credentials_service_account()` checks explicitly that `type` is `"service_account"`. This makes it easier to detect a common mistake, where the JSON for an OAuth client is provided instead of the JSON representing a service account (#93).

* `credentials_gce()` gains `cloud-platform` as a default scope, assuming that the typical user wants to "View and manage your data across Google Cloud Platform services" (#110, @MarkEdmondson1234).

# gargle 0.3.1

* [Non-interactive auth](https://gargle.r-lib.org/articles/non-interactive-auth.html) is a new vignette that serves as a guide for any client packages that use gargle for auth.

* `credentials_gce()` might actually work now (#97, @wlongabaugh).

* `credentials_app_default()` got a small bug fix relating to putting the token in the header (r-dbi/bigrquery#336)

* `token_fetch()` silently catches warnings, in addition to errors, as it falls
  through the registry of credential-fetching methods (#89).

* The yes/no asking if it's OK to cache OAuth tokens prints fully now
  (r-dbi/bigrquery#333).

# gargle 0.3.0

* The unexported functions available for generating standardized docs for
  `PKG_auth` functions in client packages have been updated.

* `token_userinfo()`, `token_email()`, and `token_tokeninfo()` are newly
  exported helpers that retrieve information for a token.

* `AuthState$set_app()` and `AuthState$set_api_key()` now allow setting a value
  of `NULL`, i.e. these fields are easier to clear.

* `credentials_byo_oauth2()` gains the ability to ingest a token from an object of class `httr::request`, i.e. to retrieve the `auth_token` component that holds an object of class `httr::Token2.0` that has been processed with `httr::config()`.

# gargle 0.2.0

* All built-in API credentials have been rotated and are stored internally in a way that reinforces appropriate use. There is a new [Privacy policy](https://tidyverse.org/google_privacy_policy/) as well as a [policy for authors of packages or other applications](https://tidyverse.org/google_privacy_policy/#policies-for-authors-of-packages-or-other-applications). This is related to a process to get the gargle project verified, which affects the OAuth2 capabilities and the consent screen.

* New vignette on "How to get your own API credentials", to help other package authors or users obtain their own API key or OAuth client ID and secret.

* `credentials_byo_oauth2()` is a new credential function. It is included in the
default registry consulted by `token_fetch()` and is tried just before
`credentials_user_oauth2()`.

# gargle 0.1.3

* Initial CRAN release
