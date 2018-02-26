`.functionkeep` <-
function(object, AIC)
{
    list(df.resid=object$df.resid, deviance=object$deviance, term=as.character(object$formula)[3], AIC=AIC)
}

