## ----load_library-------------------------------------------------------------
library(ecospat)
citation("ecospat")

## -----------------------------------------------------------------------------
data(ecospat.testData)
names(ecospat.testData)

## -----------------------------------------------------------------------------
data(ecospat.testNiche.inv)
names(ecospat.testNiche.inv)

## -----------------------------------------------------------------------------
data(ecospat.testNiche.nat)
names(ecospat.testNiche.nat)

## -----------------------------------------------------------------------------
if(requireNamespace("ape")){
  fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
  tree<-ape::read.tree(fpath)
  tree$tip.label
  plot(tree, cex=0.6)
}

## ----mantel_cor---------------------------------------------------------------
ecospat.mantel.correlogram(dfvar=ecospat.testData[c(2:16)],colxy=1:2, n=100, 
                           colvar=3:7, max=1000, nclass=10, nperm=100)

## -----------------------------------------------------------------------------
colvar <- ecospat.testData[c(4:8)]
x <- cor(colvar, method="pearson")
ecospat.npred (x, th=0.75)

## -----------------------------------------------------------------------------
x <- cor(colvar, method="spearman")
ecospat.npred (x, th=0.75)

## -----------------------------------------------------------------------------
x <- ecospat.testData[c(4:8)]
p<- x[1:90,] #A projection dataset.
ref<- x[91:300,] # A reference dataset

## -----------------------------------------------------------------------------
ecospat.climan(ref,p)

## -----------------------------------------------------------------------------
x <- ecospat.testData[c(2,3,4:8)]
proj<- x[1:90,] #A projection dataset.
cal<- x[91:300,] #A calibration dataset
mess.object<-ecospat.mess (proj, cal, w="default")
ecospat.plot.mess (mess.object, cex=1, pch=15)

## -----------------------------------------------------------------------------
if(requireNamespace("ape")){
  fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
  tree <- ape::read.tree(fpath)
  data <- ecospat.testData[9:52]
  pd<- ecospat.calculate.pd(tree, data, method = "spanning", type = "species", root = TRUE, average = FALSE, verbose = FALSE )
  plot(pd)
}

## -----------------------------------------------------------------------------
library(ade4)
inv <- ecospat.testNiche.inv
nat <- ecospat.testNiche.nat
pca.env <- ade4::dudi.pca(rbind(nat,inv)[,3:10],scannf=F,nf=2) 
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

## -----------------------------------------------------------------------------
# PCA scores for the whole study area
scores.globclim <- pca.env$li

# PCA scores for the species native distribution
scores.sp.nat <- ade4::suprow(pca.env,nat[which(nat[,11]==1),3:10])$li

# PCA scores for the species invasive distribution
scores.sp.inv <- ade4::suprow(pca.env,inv[which(inv[,11]==1),3:10])$li

# PCA scores for the whole native study area
scores.clim.nat <- ade4::suprow(pca.env,nat[,3:10])$li

# PCA scores for the whole invaded study area
scores.clim.inv <- ade4::suprow(pca.env,inv[,3:10])$li

## -----------------------------------------------------------------------------
# gridding the native niche
grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim, 
                                       glob1=scores.clim.nat,
                                       sp=scores.sp.nat, R=100,
                                       th.sp=0) 

## -----------------------------------------------------------------------------
# gridding the invasive niche
grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.inv,
                                       sp=scores.sp.inv, R=100,
                                       th.sp=0) 

## -----------------------------------------------------------------------------
# Compute Schoener's D, index of niche overlap
D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor = TRUE)$D 
D.overlap

## -----------------------------------------------------------------------------
eq.test <- ecospat.niche.equivalency.test(grid.clim.nat, grid.clim.inv,rep=10,
                                          intersection = 0.1,
                                          overlap.alternative =  "higher",
                                          expansion.alternative = "lower",
                                          stability.alternative = "higher",
                                          unfilling.alternative = "lower")

## -----------------------------------------------------------------------------
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")

## -----------------------------------------------------------------------------
sim.test <- ecospat.niche.similarity.test(grid.clim.nat, grid.clim.inv,rep=10,
                                          overlap.alternative =  "higher",
                                          expansion.alternative = "lower",
                                          stability.alternative = "higher",
                                          unfilling.alternative = "lower",
                                          intersection = 0.1,
                                          rand.type=1) 

## -----------------------------------------------------------------------------
ecospat.plot.overlap.test(sim.test, "D", "Similarity")

## ----niche.dyn----------------------------------------------------------------
niche.dyn <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv)

## -----------------------------------------------------------------------------
ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv, quant=0.25, interest=2,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")

ecospat.shift.centroids(scores.sp.nat, scores.sp.inv, scores.clim.nat, scores.clim.inv)

## -----------------------------------------------------------------------------
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity")
ecospat.plot.overlap.test(sim.test, "stability", "Similarity")
ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")

## -----------------------------------------------------------------------------
# gridding the native niche
grid.clim.t.nat <- ecospat.grid.clim.dyn(glob=as.data.frame(rbind(nat,inv)[,10]), 
                                         glob1=as.data.frame(nat[,10]),
                                    sp=as.data.frame(nat[which(nat[,11]==1),10]),
                                    R=1000, th.sp=0) 
# gridding the invaded niche
grid.clim.t.inv <- ecospat.grid.clim.dyn(glob=as.data.frame(rbind(nat,inv)[,10]),
                                         glob1=as.data.frame(inv[,10]),
                                         sp=as.data.frame(inv[which(inv[,11]==1),10]), 
                                         R=1000, th.sp=0) 
t.dyn<-ecospat.niche.dyn.index (grid.clim.t.nat, grid.clim.t.inv)
ecospat.plot.niche.dyn(grid.clim.t.nat, grid.clim.t.inv, quant=0, 
                       interest=2, title= "Niche Overlap", 
                       name.axis1="Average temperature")

## ----co_occ-------------------------------------------------------------------
data <- ecospat.testData[c(9:16,54:57)]

## -----------------------------------------------------------------------------
ecospat.co_occurrences (data)

## ----Cscore-------------------------------------------------------------------
data<- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
nperm <- 100
outpath <- getwd()
ecospat.Cscore(data, nperm, outpath)

## ----cor-plot-----------------------------------------------------------------
data <- ecospat.testData[,4:8]
ecospat.cor.plot(data)

## -----------------------------------------------------------------------------
data <- ecospat.testData
caleval <- ecospat.caleval (data = ecospat.testData[53], xy = data[2:3],
                            row.num = 1:nrow(data), nrep = 2, ratio = 0.7,
                            disaggregate = 0.2, pseudoabs = 100, npres = 10,
                            replace = FALSE)
head(caleval)

## -----------------------------------------------------------------------------
fit <- ecospat.testData$glm_Saxifraga_oppositifolia

## -----------------------------------------------------------------------------
obs<-ecospat.testData$glm_Saxifraga_oppositifolia[which(ecospat.testData$Saxifraga_oppositifolia==1)]


## ----boyce--------------------------------------------------------------------
ecospat.boyce (fit, obs, nclass = 0, window.w = "default", res = 100, 
               PEplot = TRUE)$cor

## -----------------------------------------------------------------------------
eval<-ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
pred<-ecospat.testData[c(73:92)]


## -----------------------------------------------------------------------------
CommunityEval<-ecospat.CommunityEval (eval, pred, proba = TRUE, ntir=5,verbose = T)


## -----------------------------------------------------------------------------
library(biomod2)

## -----------------------------------------------------------------------------
# species
# occurrences
xy <- inv[,1:2]
head(xy)
sp_occ <- inv[11]

# env
current <- inv[3:7]
head(current)

## BIOMOD
t1 <- Sys.time()
sp<-1

## -----------------------------------------------------------------------------
### Formating the data with the BIOMOD_FormatingData() function form the package biomod2

myBiomodData <- biomod2::BIOMOD_FormatingData( resp.var = as.numeric(sp_occ[,sp]),
                                      expl.var = current,
                                      resp.xy = xy,
                                      resp.name = colnames(sp_occ)[sp])

## ----ESM.Modeling-------------------------------------------------------------
### Calibration of simple bivariate models

# remove insivible(capture.output)) to print output in the console
# this is just to keep the vignette short
invisible(capture.output(my.ESM <- ecospat.ESM.Modeling( data=myBiomodData,
                                models=c('GLM'),
                                NbRunEval=2,
                                DataSplit=70,
                                weighting.score=c("AUC"),
                                parallel=F)
        )
)

## ----ESM.EnsembleModeling-----------------------------------------------------
### Evaluation and average of simple bivariate models to ESMs
my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM,weighting.score=c("SomersD"),threshold=0)

## ----ESM.Projection-----------------------------------------------------------
### Projection of simple bivariate models into new space 
my.ESM_proj_current <- ecospat.ESM.Projection(ESM.modeling.output=my.ESM,
                                              new.env=current)

## ----ESM.EnsembleProjection---------------------------------------------------
### Projection of calibrated ESMs into new space 
my.ESM_EFproj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=my.ESM_proj_current,
                                                        ESM.EnsembleModeling.output=my.ESM_EF)

## -----------------------------------------------------------------------------
proba <- ecospat.testData[,73:92]

## -----------------------------------------------------------------------------
sr <- as.data.frame(rowSums(proba))

## ----SESAM--------------------------------------------------------------------
prr<-ecospat.SESAM.prr(proba, sr)
head(prr)[,1:4]

## -----------------------------------------------------------------------------
presence<-ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
pred<-ecospat.testData[c(73:92)]

## -----------------------------------------------------------------------------
nbpermut <- 100
outpath <- getwd()
ecospat.cons_Cscore(presence, pred, nbpermut, outpath)

