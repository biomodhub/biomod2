## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(fig.align = "center", 
               out.width = "90%",
               fig.width = 6, fig.height = 5,
               dev.args = list(pointsize=10),
               par = TRUE, # needed for setting hook 
               collapse = TRUE, # collapse input & output code in chunks
               warning = FALSE)

knit_hooks$set(par = function(before, options, envir)
  { if(before && options$fig.show != "none") 
       par(family = "sans", mar=c(4.1,4.1,1.1,1.1), mgp=c(3,1,0), tcl=-0.5)
})
set.seed(1) # for exact reproducibility

## ----message = FALSE, echo=-2-------------------------------------------------
library(mclust)
cat(mclust:::mclustStartupMessage(), sep="")

## -----------------------------------------------------------------------------
data(diabetes)
#-- data corrections (see help("diabetes", package = "mclust")) --#
names(diabetes) <- c("class", "fpg", "glucose", "insulin")
diabetes$glucose[104] <- 455
#-----------------------------------------------------------------#
class <- diabetes$class
table(class)
X <- diabetes[,-1]
head(X)
clPairs(X, class)

BIC <- mclustBIC(X)
plot(BIC)
summary(BIC)

mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)

plot(mod1, what = "classification")
table(class, mod1$classification)

plot(mod1, what = "uncertainty")

ICL <- mclustICL(X)
summary(ICL)
plot(ICL)

LRT <- mclustBootstrapLRT(X, modelName = "VVV")
LRT

## -----------------------------------------------------------------------------
(hc1 <- hc(X, modelName = "VVV", use = "SVD"))
BIC1 <- mclustBIC(X, initialization = list(hcPairs = hc1)) # default 
summary(BIC1)

(hc2 <- hc(X, modelName = "VVV", use = "VARS"))
BIC2 <- mclustBIC(X, initialization = list(hcPairs = hc2))
summary(BIC2)

(hc3 <- hc(X, modelName = "EEE", use = "SVD"))
BIC3 <- mclustBIC(X, initialization = list(hcPairs = hc3))
summary(BIC3)

## -----------------------------------------------------------------------------
BIC <- mclustBICupdate(BIC1, BIC2, BIC3)
summary(BIC)
plot(BIC)

## ----echo=-1------------------------------------------------------------------
set.seed(20181116)
data(galaxies, package = "MASS") 
galaxies <- galaxies / 1000
BIC <- NULL
for(j in 1:20)
{
  rBIC <- mclustBIC(galaxies, verbose = FALSE,
                    initialization = list(hcPairs = hcRandomPairs(galaxies)))
  BIC <- mclustBICupdate(BIC, rBIC)
}
summary(BIC)
plot(BIC)
mod <- Mclust(galaxies, x = BIC)
summary(mod)

## -----------------------------------------------------------------------------
data(iris)
class <- iris$Species
table(class)
X <- iris[,1:4]
head(X)
mod2 <- MclustDA(X, class, modelType = "EDDA")
summary(mod2)
plot(mod2, what = "scatterplot")
plot(mod2, what = "classification")

## -----------------------------------------------------------------------------
data(banknote)
class <- banknote$Status
table(class)
X <- banknote[,-1]
head(X)
mod3 <- MclustDA(X, class)
summary(mod3)
plot(mod3, what = "scatterplot")
plot(mod3, what = "classification")

## -----------------------------------------------------------------------------
cv <- cvMclustDA(mod2, nfold = 10)
str(cv)
unlist(cv[3:6])
cv <- cvMclustDA(mod3, nfold = 10)
str(cv)
unlist(cv[3:6])

## -----------------------------------------------------------------------------
data(acidity)
mod4 <- densityMclust(acidity)
summary(mod4)
plot(mod4, what = "BIC")
plot(mod4, what = "density", data = acidity, breaks = 15)
plot(mod4, what = "diagnostic", type = "cdf")
plot(mod4, what = "diagnostic", type = "qq")

## -----------------------------------------------------------------------------
data(faithful)
mod5 <- densityMclust(faithful)
summary(mod5)
plot(mod5, what = "BIC")
plot(mod5, what = "density", type = "hdr", data = faithful, points.cex = 0.5)
plot(mod5, what = "density", type = "persp")

## -----------------------------------------------------------------------------
boot1 <- MclustBootstrap(mod1, nboot = 999, type = "bs")
summary(boot1, what = "se")
summary(boot1, what = "ci")

## ----echo=-1, fig.width=6, fig.height=7---------------------------------------
par(mfrow=c(4,3))
plot(boot1, what = "pro")
plot(boot1, what = "mean")

## -----------------------------------------------------------------------------
boot4 <- MclustBootstrap(mod4, nboot = 999, type = "bs")
summary(boot4, what = "se")
summary(boot4, what = "ci")

## ----echo=-1------------------------------------------------------------------
par(mfrow=c(2,2))
plot(boot4, what = "pro")
plot(boot4, what = "mean")

## -----------------------------------------------------------------------------
mod1dr <- MclustDR(mod1)
summary(mod1dr)
plot(mod1dr, what = "pairs")
plot(mod1dr, what = "boundaries", ngrid = 200)

mod1dr <- MclustDR(mod1, lambda = 1)
summary(mod1dr)
plot(mod1dr, what = "scatterplot")
plot(mod1dr, what = "boundaries", ngrid = 200)

## -----------------------------------------------------------------------------
mod2dr <- MclustDR(mod2)
summary(mod2dr)
plot(mod2dr, what = "scatterplot")
plot(mod2dr, what = "boundaries", ngrid = 200)

mod3dr <- MclustDR(mod3)
summary(mod3dr)
plot(mod3dr, what = "scatterplot")
plot(mod3dr, what = "boundaries", ngrid = 200)

## -----------------------------------------------------------------------------
mclust.options("bicPlotColors")
mclust.options("classPlotColors")

## ----eval=FALSE---------------------------------------------------------------
# palette.colors(palette = "Okabe-Ito")

## -----------------------------------------------------------------------------
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
"#D55E00", "#CC79A7", "#999999")

## -----------------------------------------------------------------------------
bicPlotColors <- mclust.options("bicPlotColors")
bicPlotColors[1:14] <- c(cbPalette, cbPalette[1:5])
mclust.options("bicPlotColors" = bicPlotColors)
mclust.options("classPlotColors" = cbPalette[-1])

clPairs(iris[,-5], iris$Species)
mod <- Mclust(iris[,-5])
plot(mod, what = "BIC")
plot(mod, what = "classification")

## -----------------------------------------------------------------------------
sessionInfo()

