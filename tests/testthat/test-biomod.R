library(readr)
library(dplyr)
library(rlang)
library(raster)
# library(biomod2)
devtools::load_all(".")

# test with binary respone variable?
binaryResp <- FALSE

# species occurrences
species.dat <-
  read_csv(
    system.file("external/species/mammals_table.csv", package = "biomod2")
  )

head(species.dat)

# the name of studied species
resp.name <- 'GuloGulo'

# the presence/absences data for our species
resp.var <- species.dat %>% pull(resp.name)

if (isFALSE(binaryResp)) {
  resp.var[resp.var == 1] <- round(rnorm(sum(resp.var == 1), mean = 0.2, sd = 0.1), 2)
  resp.var[resp.var > 1] <- 1
  resp.var[resp.var < 0] <- 0

  ## convert to something like abundance
  resp.var <- resp.var*10000
}

# the XY coordinates of species data
resp.xy <- species.dat %>% select_at(c('X_WGS84', 'Y_WGS84'))

# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
expl.name <- paste0('bio', c(3, 4, 7, 11, 12), '.grd')
expl.var <-
  purrr::map(
    expl.name,
    ~ raster(
      system.file(file.path('external/bioclim/current', .x), package="biomod2"),
      rat = FALSE)
  ) %>%
  stack()

# 1. Formatting Data
bm.formdat <-
  BIOMOD_FormatingData(
    binaryResp = binaryResp,
    resp.var = resp.var,
    expl.var = expl.var,
    resp.xy = resp.xy,
    resp.name = resp.name
  )

# 2. Defining Models Options using default options.
bm.opt <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation
if (binaryResp) {
  bm.opt <- BIOMOD_ModelingOptions()
  models <- c('SRE','RF', "CTA", 'MAXENT.Phillips.2')
  metric.eval <- c('TSS','ROC')
} else {
  metric.eval <- c("R2", "RMSE")
  models <- c(
    "GLM",
    "RF",
    "GBM",
    "MARS",
    "GAM"
  )
  bm.opt <- BIOMOD_ModelingOptions(
                                   GLM = list(family = Gamma(link = log)),
                                   RF = list(do.classif = FALSE),
                                   MARS = list(glm = list(family = binomial), interaction.level = 1),
                                   GAM = list(k = 3)
                                   )
}
set.seed(123)
bm.mod <-
  BIOMOD_Modeling(
    bm.formdat,
    modeling.id = "test",
    models = models,
    bm.options = bm.opt,
    nb.rep = 2,
    data.split.perc = 80,
    var.import = 0,
    metric.eval = metric.eval,
    do.full.models = FALSE,
  )

## print a summary of modeling stuff
bm.mod

bm.maxent.mod.list <-
  BIOMOD_LoadModels(bm.mod, models = 'RF')

bm.maxent.mod.1 <- get(bm.maxent.mod.list[1])

bm_VariablesImportance(bm.model = bm.maxent.mod.1,
                       expl.var = as.data.frame(expl.var),
                       method = "full_rand",
                       nb.rep = 3)

# 4.1 Projection on current environmental conditions
if (binaryResp) {
  metric.binary <- 'TSS'
} else {
  metric.binary <- NULL
}
bm.proj <-
  BIOMOD_Projection(
    bm.mod = bm.mod,
    new.env = expl.var,
    proj.name = 'current',
    models.chosen = 'all',
    metric.binary = metric.binary,
    compress = FALSE,
    build.clamping.mask = FALSE,
    on_0_1000 = FALSE
  )

if (!binaryResp) {
  metric.eval <- "R2"
}


metric.select.thresh <- rep(0.3, length(metric.eval))

set.seed(123)
bm.ensemb <- BIOMOD_EnsembleModeling(bm.mod,
                                     em.by = "all",
                                     metric.select = metric.eval,
                                     metric.select.thresh = metric.select.thresh,
                                     metric.eval = metric.eval)

bm.ensembProj <- BIOMOD_EnsembleForecasting(bm.em = bm.ensemb,
                                            bm.proj = bm.proj,
                                            on_0_1000 = FALSE)
plot(bm.ensembProj)

# bm.maxent.mod.list <-
#   BIOMOD_LoadModels(bm.mod, models='MAXENT.Phillips.2')
#
# bm.maxent.mod.1 <- get(bm.maxent.mod.list[1])
#
# bm.maxent.proj.1a <- predict(bm.maxent.mod.1, newdata = expl.var)

## TEST ALL MODELS --------------------------------------------------
# Load species occurrences (6 species available)
myFile <- system.file('external/species/mammals_table.csv', package = 'biomod2')
DataSpecies <- read.csv(myFile, row.names = 1)
head(DataSpecies)

# Select the name of the studied species
myRespName <- 'GuloGulo'

# Get corresponding presence/absence data
myResp <- as.numeric(DataSpecies[, myRespName])

# Get corresponding XY coordinates
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]

# Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myFiles <- paste0('external/bioclim/current/bio', c(3, 4, 7, 11, 12), '.grd')
myExpl <- stack(system.file(myFiles, package = 'biomod2'))


# ---------------------------------------------------------------
# Format Data with true absences
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

# Create default modeling options
myBiomodOptions <- BIOMOD_ModelingOptions()


# ---------------------------------------------------------------
# Model single models
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = 'AllModels',
                                    bm.options = myBiomodOptions,
                                    nb.rep = 2,
                                    data.split.perc = 80,
                                    metric.eval = c('TSS','ROC'),
                                    var.import = 3,
                                    do.full.models = FALSE)

myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.7),
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 3,
                                      prob.mean = TRUE,
                                      prob.median = TRUE,
                                      prob.cv = TRUE,
                                      prob.ci = TRUE,
                                      prob.ci.alpha = 0.05,
                                      committee.averaging = TRUE,
                                      prob.mean.weight = TRUE,
                                      prob.mean.weight.decay = 'proportional')
## ---------------------------------------------------------------
## Evaluate models with Boyce index and MPA
myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut,
                                  bm.em = myBiomodEM)
myBiomodPO$Cutoff

## Evaluate models with Boyce index and MPA (using background data)
myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut,
                                  bm.em = myBiomodEM,
                                  bg.env = myExpl[])

## ----------------------------------------------------------------
## Future projections

## project current first
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  proj.name = 'Current',
                                  new.env = myExpl,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)

## ---------------------------------------------------------------
## Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myFiles = paste0('external/bioclim/future/bio', c(3, 4, 7, 11, 12), '.grd')
myExplFuture = stack(system.file(myFiles, package = 'biomod2'))

## Project onto future conditions
myBiomodProjectionFuture <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                              proj.name = 'Future',
                                              new.env = myExplFuture,
                                              models.chosen = 'all',
                                              metric.binary = 'TSS',
                                              build.clamping.mask = TRUE)

## Load current and future binary projections
proj.current <- stack("GuloGulo/proj_Current/proj_Current_GuloGulo_TSSbin.grd")
proj.future <- stack("GuloGulo/proj_Future/proj_Future_GuloGulo_TSSbin.grd")

## Compute differences
myBiomodRangeSize <- BIOMOD_RangeSize(proj.current = proj.current, proj.future = proj.future)

myBiomodRangeSize$Compt.By.Models
plot(myBiomodRangeSize$Diff.By.Pixel)
