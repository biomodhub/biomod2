library(readr)
library(dplyr)
library(rlang)
# library(biomod2)
devtools::load_all(".")

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
  BIOMOD_FormatingData(
    resp.var = resp.var,
    expl.var = expl.var,
    resp.xy = resp.xy,
    resp.name = resp.name
  )

# 2. Defining Models Options using default options.
bm.opt <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation
bm.mod <-
  BIOMOD_Modeling(
    bm.formdat,
    modeling.id = "test",
    models = c('SRE','RF', 'MAXENT.Phillips.2'),
    bm.options = bm.opt,
    nb.rep = 2,
    data.split.perc = 80,
    var.import = 0,
    metric.eval = c('TSS','ROC'),
    do.full.models = FALSE,
  )

## print a summary of modeling stuff
bm.mod

# 4.1 Projection on current environmental conditions

bm.proj <-
  BIOMOD_Projection(
    bm.mod = bm.mod,
    new.env = expl.var,
    proj.name = 'current',
    models.chosen = 'all',
    metric.binary = 'TSS',
    compress = FALSE,
    build.clamping.mask = FALSE
  )


# bm.maxent.mod.list <-
#   BIOMOD_LoadModels(bm.mod, models='MAXENT.Phillips.2')
#
# bm.maxent.mod.1 <- get(bm.maxent.mod.list[1])
#
# bm.maxent.proj.1a <- predict(bm.maxent.mod.1, newdata = expl.var)

