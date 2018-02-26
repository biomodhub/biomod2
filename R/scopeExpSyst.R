`.scopeExpSyst` <-
function(enviroTrain, mod)
{
    i <- 1
    junk2 <- c()
    while(i <= dim(enviroTrain)[2]) {
        
	      vname <- names(enviroTrain)[i]
	      
        if(mod=="NNET" | mod=="FDA" | mod=="GLMs" | mod=="CTA" | mod=="GBM") junk <- vname
        if(mod == "GLMq") {
            if(is.numeric(enviroTrain[,i]))      junk <- paste(vname, "+I(", vname, "^2)+I(",vname, "^3)", sep="")
            else if(is.factor(enviroTrain[,i]))  junk <- vname
        }
        if(mod == "GLMp") {
            if(is.numeric(enviroTrain[,i]))     junk <- paste(vname, "+I(", vname, "^2)+I(",vname, "^3)+", "poly(", vname, ",2) + poly(", vname, ",3)", sep="")
            else if(is.factor(enviroTrain[,i])) junk <- vname
        }
        junk2 <- c(junk2, junk)
        i <- i + 1
    }

    junk2 <- eval(parse(text=paste("~", paste(junk2, collapse="+"))))
    return(junk2)
}

