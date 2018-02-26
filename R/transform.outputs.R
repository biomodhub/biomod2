setGeneric( ".transform.outputs", 
            def = function(modOut, out = 'evaluation',...){
                    standardGeneric( ".transform.outputs" )
                    } )

setMethod('.transform.outputs', signature(modOut='array'), 
  function(modOut, out = 'evaluation'){
    # check out attr
    if(!(out %in% c('evaluation', 'prediction', 'var.import', 'calib.failure', 'models.run', 'prediction.eval' ) )){
      stop(paste("out argument must be one of ", toString(c('evaluation', 'prediction', 'var.import',
                                                            'calib.failure', 'models.run', 'prediction.eval' ))))
    }
    
    # check dim of input list
    if(length(dim(modOut)) != 4 ){
      cat('\n',dim(modOut),'\n')
      print(dimnames(modOut))
      warning("Not computed .transform.outputs because of an imcompatible input list dimention", immediate=T)
      return(NULL)
    }

    if(dim(modOut)[4] == 1 & length(unlist(strsplit(unlist(dimnames(modOut)[4]),'_'))) == 1 ){
      dataset.names <- 'AllData'
    } else{
      if(length(dimnames(modOut)[[4]]) > 0){
        dataset.names <- unlist(sapply(unlist(dimnames(modOut)[4]), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
      } else {
        dataset.names <- paste('PA', 1:dim(modOut)[4])
      }
    }

    run.eval.names <- sub('_','',unlist(dimnames(modOut)[3]))
    mod.names <- unlist(dimnames(modOut)[2])
    
    if (out=='evaluation'){
      if( is.null(modOut['evaluation',1,1,1])){ return(NULL) }
      eval.meth.names <- rownames(as.data.frame(modOut['evaluation',1,1,1]))
      eval.col.names <- colnames(as.data.frame(modOut['evaluation',1,1,1]))
  
      eval.out <- array(data = unlist(modOut['evaluation',,,]),
                        dim = c(length(eval.meth.names),
                                length(eval.col.names),
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(eval.meth.names,
                                         eval.col.names,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(eval.out) 
    }
    
    if (out=='prediction'){
      if( is.null(modOut['pred',1,1,1])){ return(NULL) }
      nb.pts.pred <- length(as.numeric(unlist(modOut['pred',1,1,1])))
      pred.out <- array(data = unlist(modOut['pred',,,]),
                        dim = c(nb.pts.pred,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(NULL,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(pred.out) 
    }
    
    if (out=='prediction.eval'){
      if( is.null(modOut['pred.eval',1,1,1])){ return(NULL) }
      nb.pts.pred.eval <- length(as.numeric(unlist(modOut['pred.eval',1,1,1])))
      pred.eval.out <- array(data = unlist(modOut['pred.eval',,,]),
                        dim = c(nb.pts.pred.eval,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(NULL,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(pred.eval.out) 
    }
    
    if (out=='var.import'){
      if( is.null(unlist(modOut['var.import',1,1,1]))){ return(NULL) }
      nb.var <- length(as.numeric(unlist(modOut['var.import',1,1,1])))
  
      vi.out <- array(data = unlist(modOut['var.import',,,]),
                        dim = c(nb.var,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(paste('Var',1:nb.var,sep=''), # to change
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(vi.out) 
    }
    
    if (out == 'calib.failure'){
      cf.out <- unlist(modOut['calib.failure',,,])
      return(cf.out[!is.null(cf.out)])
    }
    
    if (out == 'models.run'){
      mod.run.out <- unlist(modOut['ModelName',,,])
      return(mod.run.out[!is.null(mod.run.out)])
    }
    
  })

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
setMethod('.transform.outputs', signature(modOut='list'), 
  function(modOut, out = 'evaluation', dim.names = NULL){

    # check out attr
    if(!(out %in% c('evaluation', 'prediction', 'prediction.eval', 'var.import', 'calib.failure',
                    'models.run', 'EF.prediction', 'EF.PCA.median', 'EF.evaluation' ) )){
      stop(paste("out argument must be one of ", toString(c('evaluation', 'prediction', 'prediction.eval', 'var.import',
                                                            'calib.failure', 'models.run', 'EF.prediction',
                                                            'EF.PCA.median', 'EF.evaluation'))))
    }
    
    if(length(modOut) == 1 & length(unlist(strsplit(unlist(names(modOut)),'_'))) == 1 ){
      dataset.names <- 'AllData'
    } else{
      if(is.null(dim.names)){
        dataset.names <- unlist(sapply(unlist(names(modOut)), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
      } else{
        dataset.names <- unlist(dim.names[1])
      }
    }

    if(is.null(dim.names)){
      run.eval.names <- sub('_','',unlist(names(modOut[[1]]))) # may be good here to test that all names are identics
    
      mod.names <- unlist(names(modOut[[1]][[1]]))
    } else{
      run.eval.names <- unlist(dim.names[2])
      mod.names <- unlist(dim.names[3])
    }
    
    if (out=='evaluation'){
      
      eval.tab <- NULL
      nb_pa <- length(modOut)
      nb_run <- length(modOut[[1]])
      nb_mod <- length(modOut[[1]][[1]])
      
      for(i in 1:nb_pa){
        for(j in 1:nb_run){
          for(k in 1:nb_mod){
            eval.tab <- modOut[[i]][[j]][[k]][['evaluation']]
            if(!is.null(eval.tab)){ break }
          }
          if(!is.null(eval.tab)){ break }  
        }
        if(!is.null(eval.tab)){ break }
      }
      
      if( is.null(eval.tab)){ return(NULL) }
  
      eval.meth.names <- rownames(as.data.frame(eval.tab))
      eval.col.names <- colnames(as.data.frame(eval.tab))
      
      eval.out <- lapply(names(modOut),function(d1){ # data set
                    lapply(names(modOut[[d1]]), function(d2){ # run eval
                      lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
                        if(is.null(modOut[[d1]][[d2]][[d3]][['calib.failure']])){
                          return(data.frame(modOut[[d1]][[d2]][[d3]][['evaluation']]))
                        } else { matrix(NA, ncol=length(eval.col.names), nrow=length(eval.meth.names), dimnames=list(eval.meth.names,eval.col.names))}
                      })
                    })
                  })
  
      eval.out <- array(data = unlist(eval.out),
                        dim = c(length(eval.meth.names),
                                length(eval.col.names),
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(eval.meth.names,
                                         eval.col.names,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(eval.out) 
    }
    
    if (out=='prediction'){

      pred.tab <- NULL
      nb_pa <- length(modOut)
      nb_run <- length(modOut[[1]])
      nb_mod <- length(modOut[[1]][[1]])
      
      for(i in 1:nb_pa){
        for(j in 1:nb_run){
          for(k in 1:nb_mod){
            pred.tab <- modOut[[i]][[j]][[k]][['pred']]
            if(!is.null(pred.tab)){ break }
          }
          if(!is.null(pred.tab)){ break }  
        }
        if(!is.null(pred.tab)){ break }
      }
      
      if( is.null(pred.tab)){ return(NULL) }
      

      nb.pts.pred <- length(as.numeric(pred.tab))
      
      pred.out <- lapply(names(modOut),function(d1){ # data set
                    lapply(names(modOut[[d1]]), function(d2){ # run eval
                      lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
                        if(is.null(modOut[[d1]][[d2]][[d3]][['pred']])){
                          return(rep(NA,nb.pts.pred))
                        } else{
                          return(as.numeric(modOut[[d1]][[d2]][[d3]][['pred']]))
                        }
                      })
                    })
                  })
      
      pred.out <- array(data = unlist(pred.out),
                        dim = c(nb.pts.pred,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(NULL,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(pred.out) 
    }
    
    if (out=='prediction.eval'){
      pred.eval.tab <- NULL
      nb_pa <- length(modOut)
      nb_run <- length(modOut[[1]])
      nb_mod <- length(modOut[[1]][[1]])
      
      for(i in 1:nb_pa){
        for(j in 1:nb_run){
          for(k in 1:nb_mod){
            pred.eval.tab <- modOut[[i]][[j]][[k]][['pred.eval']]
            if(!is.null(pred.eval.tab)){ break }
          }
          if(!is.null(pred.eval.tab)){ break }  
        }
        if(!is.null(pred.eval.tab)){ break }
      }
      
      if( is.null(pred.eval.tab)){ return(NULL) }
      

      nb.pts.pred.eval <- length(as.numeric(pred.eval.tab))
     
      pred.eval.out <- lapply(names(modOut),function(d1){ # data set
                    lapply(names(modOut[[d1]]), function(d2){ # run eval
                      lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
                        if(is.null(modOut[[d1]][[d2]][[d3]][['pred.eval']])){
                          return(rep(NA,nb.pts.pred.eval))
                        } else{
                          return(as.numeric(modOut[[d1]][[d2]][[d3]][['pred.eval']]))
                        }
                      })
                    })
                  })
      
      pred.eval.out <- array(data = unlist(pred.eval.out),
                        dim = c(nb.pts.pred.eval,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(NULL,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(pred.eval.out) 
    }
    
    if (out=='var.import'){
      vi.tab <- NULL
      nb_pa <- length(modOut)
      nb_run <- length(modOut[[1]])
      nb_mod <- length(modOut[[1]][[1]])
      
      for(i in 1:nb_pa){
        for(j in 1:nb_run){
          for(k in 1:nb_mod){
            vi.tab <- modOut[[i]][[j]][[k]][['var.import']]
            if(!is.null(vi.tab)){ break }
          }
          if(!is.null(vi.tab)){ break }  
        }
        if(!is.null(vi.tab)){ break }
      }
      
      if( is.null(vi.tab)){ return(NULL) }                             
      
      nb.var <- length(as.numeric(unlist(vi.tab)))
      
      ef.mod <- grep(pattern="EF.",mod.names) # EF models
      if(length(ef.mod)>0){
        kept.mod <- mod.names[-ef.mod]
      } else{
        kept.mod <- mod.names
      }
                                        
      vi.out <- lapply(names(modOut),function(d1){ # data set
                  lapply(names(modOut[[d1]]), function(d2){ # run eval
                    lapply(kept.mod, function(d3){ # models without EF ones
                      if(is.null(modOut[[d1]][[d2]][[d3]][['var.import']])){
                        return(rep(NA,nb.var))
                      } else{                      
                        return(as.matrix(modOut[[d1]][[d2]][[d3]][['var.import']]))
                      }
                    })
                  })
                })                                 
  
      vi.out <- array(data = unlist(vi.out),
                        dim = c(nb.var,
                                length(kept.mod),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(names(modOut[[1]][[1]][[1]][['var.import']]), # to change
                                         kept.mod,
                                         run.eval.names,
                                         dataset.names))
      
     return(vi.out) 
    }
    
    if (out == 'calib.failure'){
      cf.out <- lapply(names(modOut),function(d1){ # data set
                  lapply(names(modOut[[d1]]), function(d2){ # run eval
                    lapply(names(modOut[[d1]][[d2]]), function(d3){ # models 
                      return(modOut[[d1]][[d2]][[d3]][['calib.failure']])
                    })
                  })
                })
      cf.out <- unlist(cf.out)
      if(length(cf.out)) cf.out <- na.omit(cf.out)
      if(length(cf.out)) cf.out <- cf.out[!is.null(cf.out)]
      if(!length(cf.out)) cf.out <- 'none'
      return(cf.out)
    }
    
    if (out == 'models.run'){
      mod.run.out <- lapply(names(modOut),function(d1){ # data set
                  lapply(names(modOut[[d1]]), function(d2){ # run eval
                    lapply(names(modOut[[d1]][[d2]]), function(d3){ # models 
                      return(as.character(modOut[[d1]][[d2]][[d3]][['ModelName']]))
                    })
                  })
                })
      mod.run.out <- unlist(mod.run.out)
      if(length(mod.run.out)) mod.run.out <- na.omit(mod.run.out)
      if(length(mod.run.out)) mod.run.out <- mod.run.out[!is.null(mod.run.out)]
      if(!length(mod.run.out)) mod.run.out <- 'none'
      return(mod.run.out)
    }
    

    if (out == 'EF.prediction'){
      if( is.null(modOut[[1]][[1]][[1]][['EM']])){ return(NULL) }

      nb.pts.ef.pred <- length(as.numeric(unlist(modOut[[1]][[1]][[1]][['EM']])))
      
      ef.pred.out <- lapply(1:length(modOut),function(d1){ # data set
                        lapply(1:length(modOut[[d1]]), function(d2){ # run eval
                          lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
                            return(as.numeric(modOut[[d1]][[d2]][[d3]][['EM']]))
                          })
                        })
                      })
      
      ef.pred.out <- array( data = unlist(ef.pred.out),
                            dim = c(nb.pts.ef.pred,
                                    length(modOut[[1]][[1]]),
                                    length(modOut[[1]]),
                                    length(modOut)),
                            dimnames = list(NULL,
                                             mod.names,
                                             run.eval.names,
                                             dataset.names))
      
     return(ef.pred.out) 
    }

    if (out == 'EF.PCA.median'){
      if( is.null(modOut[[1]][[1]][[1]][['PCA.median']])){ return(NULL) }

      ef.pca.out <- lapply(1:length(modOut),function(d1){ # data set
                        lapply(1:length(modOut[[d1]]), function(d2){ # run eval
                          lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
                            return(as.character(modOut[[d1]][[d2]][[d3]][['PCA.median']]))
                          })
                        })
                      })
      
      ef.pca.out <- array( data = unlist(ef.pca.out),
                            dim = c(1,
                                    length(modOut[[1]][[1]]),
                                    length(modOut[[1]]),
                                    length(modOut)),
                            dimnames = list(NULL,
                                             mod.names,
                                             run.eval.names,
                                             dataset.names))
      
     return(ef.pca.out) 
    }

    if (out == 'EF.evaluation'){
      if( is.null(modOut[[1]][[1]][[1]][['EM.eval']])){ return(NULL) }
      eval.meth.names <- rownames(as.data.frame(modOut[[1]][[1]][[1]][['EM.eval']]))
      eval.col.names <- colnames(as.data.frame(modOut[[1]][[1]][[1]][['EM.eval']]))
      
      ef.eval.out <- lapply(1:length(modOut),function(d1){ # data set
                    lapply(1:length(modOut[[d1]]), function(d2){ # run eval
                      lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
                        return(data.frame(modOut[[d1]][[d2]][[d3]][['EM.eval']]))
                      })
                    })
                  })
  
      ef.eval.out <- array(data = unlist(ef.eval.out),
                        dim = c(length(eval.meth.names),
                                length(eval.col.names),
                                length(modOut[[1]][[1]]),
                                length(modOut[[1]]),
                                length(modOut)),
                        dimnames = list(eval.meth.names,
                                         eval.col.names,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(ef.eval.out)
    }
    
  })

DF_to_ARRAY <- function(df){
#   cat("\n*** class(df) = ", class(df))
#   cat("\n*** colnames(df) = ", colnames(df))
  if(!is.data.frame(df) & !is.matrix(df)){
    if(is.list(df)){
      df.names <- names(df)
      df <- as.data.frame(df)
      names(df) <- df.names
    } else{
      stop("You have to give a data.frame")
    }
  }
  
  a <- sapply(strsplit(colnames(df), '_'), tail, n=3)
  b <- lapply(1:3, function(id) return(unique(a[id,])))
  array.dim.names <- c(list(character(0)),rev(b))
#   array.dim.names <- c(list(c(NULL)),rev(apply(sapply(strsplit(colnames(df), '_'), tail, n=3),1,unique)))
  
  array.dim <- c(nrow(df),sapply(array.dim.names[-1],length))
  array.out <- array(data=NA, dim=array.dim, dimnames=array.dim.names)

  for(x in colnames(df)){
    dimTmp <- rev(tail(unlist(strsplit(x, '_')), n=3))
    array.out[,dimTmp[1],dimTmp[2],dimTmp[3]] <- df[,x]
  }
  return(array.out)
}

LIST_to_ARRAY <- function(ll){
  test <- sapply(ll, is.array)
  if(!all(test)) stop("list elements should be arrays")
  test <- sapply(ll,dim)
  test <- apply(test,1,function(x){length(unique(x))==1})
  if(!all(test)) stop("list elements differ in dimension")
  
  formal.dim.names <- dimnames(ll[[1]])
  new.dim.names <- rev(apply(sapply(strsplit(names(ll), '_'), tail, n=3),1,unique))
  array.dim.names <- c(formal.dim.names,new.dim.names)
  array.dim <- sapply(array.dim.names,length)
  
  array.out <- array(data=NA, dim=array.dim, dimnames=array.dim.names)
  
  for(x in names(ll)){
    dimTmp <- rev(tail(unlist(strsplit(x, '_')), n=3))
    dimTmp <- paste( paste(rep(",",length(formal.dim.names)),collapse="") , paste("'",dimTmp,"'",sep="",collapse=","),collapse="")
    
    eval(parse(text=paste("array.out[",dimTmp,"] <-  ll[[x]]",sep="")))
  }
  return(array.out)
}
