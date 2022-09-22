#' Dummy function to clean working directory after package checks
#'
#' This is the last function that will be check and will remove residual files 
#' from other examples checking such as 'GuloGulo'. This is required to pass
#' CRAN checks.
#' @return
#' nothing returned
#' @keywords internal
#' @export
#'
#' @examples
#' 
#' zzz_bm()

zzz_bm <- function(){
  unlink('GuloGulo', recursive = TRUE)
}