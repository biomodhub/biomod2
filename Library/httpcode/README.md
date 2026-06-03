httpcode
========



[![cran checks](https://cranchecks.info/badges/worst/httpcode)](https://cranchecks.info/pkgs/httpcode)
[![Build Status](https://travis-ci.org/sckott/httpcode.svg)](https://travis-ci.org/sckott/httpcode)
[![codecov](https://codecov.io/gh/sckott/httpcode/branch/master/graph/badge.svg)](https://codecov.io/gh/sckott/httpcode)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/httpcode)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/httpcode)](https://cran.r-project.org/package=httpcode)

`httpcode` is a tiny R package to search for and show http code messages and description. It's a port of the Python httpcode (https://github.com/rspivak/httpcode) library.

`httpcode` has no dependencies.

Follows RFC 2616 (https://www.ietf.org/rfc/rfc2616.txt), and for additional codes 
following RFC 6585 (https://tools.ietf.org/html/rfc6585).

Structure of information for each status code:

* `status_code` - the status code
* `message` - very brief message explaining the code
* `explanation` - more verbose explanation, but still short
* `explanation_verbose` - the complete explanation

## Installation

Stable version


```r
install.packages("httpcode")
```

Development version


```r
remotes::install_github("sckott/httpcode")
```


```r
library("httpcode")
```

## Search by http code


```r
http_code(100)
#> <Status code: 100>
#>   Message: Continue
#>   Explanation: Request received, please continue
```


```r
http_code(400)
#> <Status code: 400>
#>   Message: Bad Request
#>   Explanation: Bad request syntax or unsupported method
```


```r
http_code(503)
#> <Status code: 503>
#>   Message: Service Unavailable
#>   Explanation: The server cannot process the request due to a high load
```


```r
http_code(999)
#> Error: No description found for code: 999
```

## Get verbose status code description


```r
http_code(100, verbose = TRUE)
#> <Status code: 100>
#>   Message: Continue
#>   Explanation: Request received, please continue
#>   Verbose Explanation: The client SHOULD continue with its request. This
#>     interim response is used to inform the client that the initial part of the
#>     request has been received and has not yet been rejected by the server. The
#>     client SHOULD continue by sending the remainder of the request or, if the
#>     request has already been completed, ignore this response. The server MUST
#>     send a final response after the request has been completed. See section
#>     8.2.3 for detailed discussion of the use and handling of this status code.
```


```r
http_code(400, verbose = TRUE)
#> <Status code: 400>
#>   Message: Bad Request
#>   Explanation: Bad request syntax or unsupported method
#>   Verbose Explanation: The request could not be understood by the server due to
#>     malformed syntax. The client SHOULD NOT repeat the request without
#>     modifications.
```

# Fuzzy code search


```r
http_code('1xx')
#> [[1]]
#> <Status code: 100>
#>   Message: Continue
#>   Explanation: Request received, please continue
#> 
#> [[2]]
#> <Status code: 101>
#>   Message: Switching Protocols
#>   Explanation: Switching to new protocol; obey Upgrade header
#> 
#> [[3]]
#> <Status code: 102>
#>   Message: Processing
#>   Explanation: This code indicates that the server has received and is processing the request, but no response is available yet (WebDAV; RFC 2518)
```


```r
http_code('3xx')
#> [[1]]
#> <Status code: 300>
#>   Message: Multiple Choices
#>   Explanation: Object has several resources -- see URI list
#> 
#> [[2]]
#> <Status code: 301>
#>   Message: Moved Permanently
#>   Explanation: Object moved permanently -- see URI list
#> 
...
```


```r
http_code('30[12]')
#> [[1]]
#> <Status code: 301>
#>   Message: Moved Permanently
#>   Explanation: Object moved permanently -- see URI list
#> 
#> [[2]]
#> <Status code: 302>
#>   Message: Found
#>   Explanation: Object moved temporarily -- see URI list
```


```r
http_code('30[34]')
#> [[1]]
#> <Status code: 303>
#>   Message: See Other
#>   Explanation: Object moved -- see Method and URL list
#> 
#> [[2]]
#> <Status code: 304>
#>   Message: Not Modified
#>   Explanation: Document has not changed since given time
```

## Search by message


```r
http_search("request")
#> [[1]]
#> <Status code: 100>
#>   Message: Continue
#>   Explanation: Request received, please continue
#> 
#> [[2]]
#> <Status code: 101>
#>   Message: Switching Protocols
#>   Explanation: Switching to new protocol; obey Upgrade header
#> 
...
```


```r
http_search("forbidden")
#> [[1]]
#> <Status code: 403>
#>   Message: Forbidden
#>   Explanation: Request forbidden -- authorization will not help
```


```r
http_search("too")
#> [[1]]
#> <Status code: 400>
#>   Message: Bad Request
#>   Explanation: Bad request syntax or unsupported method
#> 
#> [[2]]
#> <Status code: 403>
#>   Message: Forbidden
#>   Explanation: Request forbidden -- authorization will not help
#> 
#> [[3]]
#> <Status code: 413>
#>   Message: Request Entity Too Large
#>   Explanation: Entity is too large.
#> 
#> [[4]]
#> <Status code: 414>
#>   Message: Request-URI Too Long
#>   Explanation: URI is too long.
#> 
#> [[5]]
#> <Status code: 425>
#>   Message: Too Early
#>   Explanation: Indicates that the server is unwilling to risk processing a request that might be replayed.
#> 
#> [[6]]
#> <Status code: 429>
#>   Message: Too Many Requests
#>   Explanation: The user has sent too many requests in a given amount of time ("rate limiting") (RFC 6585)
#> 
#> [[7]]
#> <Status code: 431>
#>   Message: Request Header Fields Too Large
#>   Explanation: The server is unwilling to process the request because its header fields are too large. The request may be resubmitted after reducing the size of the request header fields (RFC 6585)
#> 
#> [[8]]
#> <Status code: 494>
#>   Message: Request header too large (Nginx)
#>   Explanation: Client sent too large request or too long header line.
```


```r
http_search("birds")
#> Error: No status code found for search: birds
```


## Bugs/features?

See [issues](https://github.com/sckott/httpcode/issues)

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[coc]: https://github.com/sckott/httpcode/blob/master/CODE_OF_CONDUCT.md
