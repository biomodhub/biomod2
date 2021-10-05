library(readr)
library(dplyr)
library(rlang)
# library(biomod2)
devtools::load_all(".")

# test with binary respone variable?
binaryResp <- FALSE

# species occurrences
species.dat <-
  read_csv(
    system.file("external/species/mammals_table.csv", package="biomod2")
  )

head(species.dat)

# the name of studied species
resp.name <- 'GuloGulo'

# the presence/absences data for our species
resp.var <- species.dat %>% pull(resp.name)

if (!binaryResp) {
  resp.var[resp.var == 1] <- round(rnorm(sum(resp.var == 1), mean = 0.2, sd = 0.1), 2)
  resp.var[resp.var > 1] <- 1
  resp.var[resp.var < 0] <- 0
}

# the XY coordinates of species data
resp.xy <- species.dat %>% select_at(c('X_WGS84', 'Y_WGS84'))


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
expl.name <- paste0('bio', c(3, 4, 7, 11, 12), '.grd')
expl.var <-
  purrr::map(
    expl.name,
    ~ raster::raster(
      system.file(file.path('external/bioclim/current', .x), package="biomod2"),
      rat = FALSE)
  ) %>%
  raster::stack()

# 1. Formatting Data
bm.formdat <-
  BIOMOD_FormatingData(binaryResp = binaryResp,
                       resp.var = resp.var,
                       expl.var = expl.var,
                       resp.xy = resp.xy,
                       resp.name = resp.name
  )

## binaryResp = FALSE works with:
## GLM, RF, GAM, MARS

# 2. Defining Models Options using default options.

# 3. Doing Modelisation
if (binaryResp) {
  bm.opt <- BIOMOD_ModelingOptions()
  models <- c('SRE','RF', 'MAXENT.Phillips.2')
  models.eval.meth <- c('TSS','ROC')
} else {
  models.eval.meth <- c("R2", "RMSE")
  models <- c("GLM", "RF", "MARS", "GAM")
  bm.opt <- BIOMOD_ModelingOptions(RF = list(do.classif = FALSE),
                                   GAM = list(k = 3),
                                   MARS = list(glm = list(family = binomial), interaction.level = 1))
}


bm.mod <-
  BIOMOD_Modeling(
    bm.formdat,
    models = models,
    models.options = bm.opt,
    NbRunEval = 2,
    DataSplit = 80,
    VarImport = 0,
    models.eval.meth = models.eval.meth,
    do.full.models = FALSE,
    modeling.id = "test"
  )

## print a summary of modeling stuff
bm.mod

# 4.1 Projection on current environmental conditions
if (binaryResp) {
  binary.meth <- 'TSS'
} else {
  binary.meth <- NULL
}
bm.proj <-
  BIOMOD_Projection(
    modeling.output = bm.mod,
    new.env = expl.var,
    proj.name = 'current',
    selected.models = 'all',
    binary.meth = binary.meth,
    compress = FALSE,
    build.clamping.mask = FALSE
  )

if (binaryResp) {
  eval.metric.quality.threshold <- NULL
} else {
  eval.metric.quality.threshold <- 0.3
}

bm.ensemb <- BIOMOD_EnsembleModeling(bm.mod, em.by = "all",
                                     eval.metric = models.eval.meth,
                                     eval.metric.quality.threshold = eval.metric.quality.threshold,
                                     models.eval.meth = models.eval.meth)


bm.ensembProj <- BIOMOD_EnsembleForecasting(EM.output = bm.ensemb, projection.output = bm.proj)
plot(bm.ensembProj)

# bm.maxent.mod.list <-
#   BIOMOD_LoadModels(bm.mod, models='MAXENT.Phillips.2')
#
# bm.maxent.mod.1 <- get(bm.maxent.mod.list[1])
#
# bm.maxent.proj.1a <- predict(bm.maxent.mod.1, newdata = expl.var)

