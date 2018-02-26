`.scope` <-
function(enviroTrain, Smoother, degree)
{
    XXX <- enviroTrain
    deg <- degree
    vnames <- names(XXX[])
    step.list <- as.list(vnames)
    names(step.list) <- vnames
    NbVar <- dim(enviroTrain)[2]
    i <- 1
    while(i <= NbVar) {
        vname <- names(XXX)[i]
        # loops through independent variable names
        junk <- c(paste("1 + ",vname, sep=""))
        # minimum scope
        if(is.numeric(XXX[,i])) {
            junk <- c(junk, paste(Smoother, "(", vname, ",", deg, ")", sep=""))
            junk <- eval(parse(text=paste("~", paste(junk, collapse="+"))))
        }
        else if(is.factor(XXX[,i])) {
            junk <- c(junk, paste(vname, sep=""))
            junk <- eval(parse(text=paste("~", paste(junk, collapse="+"))))
        }
        step.list[[vname]] <- junk
        i <- i + 1
    }
    
    return(step.list)
}

`.scope2` <-
function(enviroTrain, formula, Smoother, degree)
{
  # 0. args checking
  if(is.character(formula)) formula <- as.formula(formula)
  if(!inherits(formula,"formula")) stop("formula must be a formula object")
  
  if(is.matrix(enviroTrain)) enviroTrain <- as.data.frame(enviroTrain)
  
  # 1. detect factoriel variables
  factVar <- as.list(names(enviroTrain))
  factVar <- lapply(factVar, is.factor)
  names(factVar) <- names(enviroTrain)
  
  # 2. create the output squeletom
  step.list <- as.list(attr(terms(formula),"term.labels"))
  
  # 3. filling the output obj
  step.list <- lapply(step.list, function(x){
    junk <- c(paste("~1 + ",x, sep=""))
    if(length(factVar[[x]])){ # x is a simple variable
      if(!factVar[[x]]){ # x is not a factor
        junk <- paste(junk, " + ", Smoother, "(", x, ",", degree, ")", sep="")
      }
    } else{
      junk <- paste(junk, " + ", Smoother, "(", x, ",", degree, ")", sep="")
    }
    return(formula(junk))
  })
  
  names(step.list) <- attr(terms(formula),"term.labels")
  
  return(step.list)
}
