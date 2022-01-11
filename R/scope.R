# .scope <- function(enviroTrain, Smoother, degree)
# {
#   XXX <- enviroTrain
#   deg <- degree
#   vnames <- names(XXX[])
#   step.list <- as.list(vnames)
#   names(step.list) <- vnames
#   NbVar <- dim(enviroTrain)[2]
#   i <- 1
#   while (i <= NbVar)
#   {
#     vname <- names(XXX)[i]
#     # loops through independent variable names
#     junk <- paste0("1 + ", vname)
#     # minimum scope
#     if (is.numeric(XXX[, i])) {
#       junk <- c(junk, paste0(Smoother, "(", vname, ",", deg, ")"))
#       junk <- eval(parse(text = paste("~", paste(junk, collapse = "+"))))
#     } else if (is.factor(XXX[, i])) {
#       junk <- c(junk, vname)
#       junk <- eval(parse(text = paste("~", paste(junk, collapse = "+"))))
#     }
#     step.list[[vname]] <- junk
#     i <- i + 1
#   }
#   
#   return(step.list)
# }
# 
# # .scope2 <- function(enviroTrain, formula, Smoother, degree)
# # {
# #   # 0. args checking
# #   if (is.character(formula)) { formula <- as.formula(formula) }
# #   if (!inherits(formula, "formula")) { stop("formula must be a formula object") }
# #   if (is.matrix(enviroTrain)) { enviroTrain <- as.data.frame(enviroTrain) }
# #   
# #   # 1. detect factoriel variables
# #   factVar <- as.list(names(enviroTrain))
# #   factVar <- lapply(factVar, is.factor)
# #   names(factVar) <- names(enviroTrain)
# #   
# #   # 2. create the output squeletom
# #   step.list <- as.list(attr(terms(formula), "term.labels"))
# #   
# #   # 3. filling the output obj
# #   step.list <- lapply(step.list, function(x)
# #   {
# #     junk <- paste0("~1 + ", x)
# #     if (length(factVar[[x]])) { # x is a simple variable
# #       if (!factVar[[x]]) { # x is not a factor
# #         junk <- paste0(junk, " + ", Smoother, "(", x, ",", degree, ")")
# #       }
# #     } else {
# #       junk <- paste0(junk, " + ", Smoother, "(", x, ",", degree, ")")
# #     }
# #     return(formula(junk))
# #   })
# #   names(step.list) <- attr(terms(formula), "term.labels")
# #   
# #   return(step.list)
# # }
# 
# ###################################################################################################
# 
# .scopeExpSyst <- function(enviroTrain, mod)
# {
#   i <- 1
#   junk2 <- c()
#   while (i <= dim(enviroTrain)[2])
#   {
#     vname <- names(enviroTrain)[i]
#     
#     if (mod %in% c("NNET", "FDA", "GLMs", "CTA", "GBM")) {
#       junk <- vname
#     } else if (mod == "GLMq") {
#       if (is.numeric(enviroTrain[, i])) {
#         junk <- paste0(vname, "+I(", vname, "^2)+I(", vname, "^3)")
#       } else if (is.factor(enviroTrain[, i])) {
#         junk <- vname
#       }
#     } else if (mod == "GLMp") {
#       if (is.numeric(enviroTrain[, i])) {
#         junk <- paste0(vname, "+I(", vname, "^2)+I(", vname, "^3)+", "poly(", vname, ",2) + poly(", vname, ",3)")
#       } else if(is.factor(enviroTrain[, i])) {
#         junk <- vname
#       }
#     }
#     junk2 <- c(junk2, junk)
#     i <- i + 1
#   }
#   
#   junk2 <- eval(parse(text = paste("~", paste(junk2, collapse = "+"))))
#   return(junk2)
# }
