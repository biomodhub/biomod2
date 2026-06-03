## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(
  fig.width = 7, fig.height = 5, fig.align = "center",
  warning = FALSE, message = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
# install.packages("Bchron")

## -----------------------------------------------------------------------------
library(Bchron)

## ----fig.align='center',fig.width=6,fig.height=5------------------------------
ages1 <- BchronCalibrate(
  ages = 11553,
  ageSds = 230,
  calCurves = "intcal20",
  ids = "Ox-123456"
)
summary(ages1)
plot(ages1)

## ----results='hide'-----------------------------------------------------------
ages2 <- BchronCalibrate(
  ages = c(3445, 11553, 7456),
  ageSds = c(50, 230, 110),
  calCurves = c("intcal20", "intcal20", "shcal20")
)
summary(ages2)
plot(ages2)

## ----eval=FALSE---------------------------------------------------------------
# plot(ages2, date = "Date3")
# plot(ages2, date = 1:2)

## ----fig.align='center',fig.width=6,fig.height=5------------------------------
ages3 <- BchronCalibrate(
  ages = c(3445, 11000),
  ageSds = c(50, 200),
  positions = c(100, 150),
  calCurves = c("intcal20", "normal")
)
summary(ages3)
plot(ages3)

## -----------------------------------------------------------------------------
library(ggplot2)
plot(ages3) +
  labs(
    x = "Age (years BP)",
    y = "Depth (cm)",
    title = "Two dates at different depths"
  )

## -----------------------------------------------------------------------------
plot(ages3, ageScale = "bc", scaleReverse = FALSE) +
  labs(
    x = "Age (years BC/AD)",
    y = "Depth (cm)",
    title = "Two dates at different depths"
  )

## -----------------------------------------------------------------------------
plot(ages1, includeCal = TRUE)

## -----------------------------------------------------------------------------
ages4 <- BchronCalibrate(
  ages = c(3445, 11553, 7456),
  ageSds = c(50, 230, 110),
  calCurves = rep("intcal20", 3)
)
plot(ages4, includeCal = TRUE, fillCol = "orange")

## -----------------------------------------------------------------------------
# First create age samples for each date
age_samples <- sampleAges(ages3)
# Now summarise them with quantile - this gives a 95% credible interval
apply(age_samples, 2, quantile, prob = c(0.025, 0.975))

## -----------------------------------------------------------------------------
apply(age_samples, 2, quantile, prob = c(0.5))

## -----------------------------------------------------------------------------
age_diff <- age_samples[, 2] - age_samples[, 1]
qplot(age_diff,
  geom = "histogram", bins = 30,
  main = "Age difference between 3445+/-50 and 11553+/-230",
  ylab = "Frequency", xlab = "Age difference"
)

## -----------------------------------------------------------------------------
data(Glendalough)
print(Glendalough)

## ----results='hide'-----------------------------------------------------------
GlenOut <- with(
  Glendalough,
  Bchronology(
    ages = ages,
    ageSds = ageSds,
    calCurves = calCurves,
    positions = position,
    positionThicknesses = thickness,
    ids = id,
    predictPositions = seq(0, 1500, by = 10)
  )
)

## ----eval=FALSE---------------------------------------------------------------
# help(Bchronology)

## ----eval=FALSE---------------------------------------------------------------
# summary(GlenOut)

## -----------------------------------------------------------------------------
summary(GlenOut, type = "convergence")
summary(GlenOut, type = "outliers")

## ----fig.align='center',fig.width=6,fig.height=5------------------------------
plot(GlenOut)

## ----results='hide'-----------------------------------------------------------
predictAges <- predict(GlenOut,
  newPositions = c(150, 725, 1500),
  newPositionThicknesses = c(5, 0, 20)
)
predictAges <- predict(GlenOut,
  newPositions = seq(0, 1500, by = 10)
)

## ----results ='hide'----------------------------------------------------------
acc_rate <- summary(GlenOut,
  type = "acc_rate",
  probs = c(0.25, 0.5, 0.75)
)

## ----eval=FALSE---------------------------------------------------------------
# acc_rate_2 <- acc_rate %>%
#   gather(
#     key = Percentile,
#     value = value,
#     -age_grid
#   )
# ggplot(acc_rate_2, aes(
#   x = age_grid,
#   y = value,
#   colour = Percentile
# )) +
#   scale_colour_brewer(palette = 1) +
#   geom_line() +
#   theme_bw() +
#   scale_x_reverse() +
#   labs(
#     y = "cm per year",
#     x = "Age (cal years BP)"
#   )

## ----eval=FALSE---------------------------------------------------------------
# sed_rate <- summary(GlenOut,
#   type = "sed_rate", useExisting = FALSE,
#   probs = c(0.25, 0.5, 0.75)
# )

## ----eval=FALSE---------------------------------------------------------------
# sed_rate_2 <- sed_rate %>%
#   gather(
#     key = Percentile,
#     value = value,
#     -position_grid
#   )
# ggplot(sed_rate_2, aes(
#   x = position_grid,
#   y = value,
#   colour = Percentile
# )) +
#   scale_colour_brewer(palette = 1) +
#   geom_line() +
#   theme_bw() +
#   scale_x_reverse() +
#   labs(
#     y = "Years per cm",
#     x = "Depth (cm)"
#   )

## -----------------------------------------------------------------------------
summary(GlenOut, type = "max_var")

## -----------------------------------------------------------------------------
max_var <- summary(GlenOut, type = "max_var", numPos = 1)
plot(GlenOut) +
  labs(
    title = "Glendalough",
    x = "Age (cal years BP)",
    y = "Depth (cm)"
  ) +
  geom_hline(yintercept = max_var)

## ----eval=FALSE---------------------------------------------------------------
# dateInfluence(GlenOut,
#   whichDate = "Beta-100901",
#   measure = "absMedianDiff"
# )

## ----eval=FALSE---------------------------------------------------------------
# dateInfluence(GlenOut,
#   whichDate = "all",
#   measure = "absMedianDiff"
# )

## ----eval=FALSE---------------------------------------------------------------
# dateInfluence(GlenOut,
#   whichDate = "all",
#   measure = "KL"
# )

## -----------------------------------------------------------------------------
data(TestChronData)
data(TestRSLData)

## ----messages=FALSE, results='hide', eval=FALSE-------------------------------
# RSLchron <- with(
#   TestChronData,
#   Bchronology(
#     ages = ages,
#     ageSds = ageSds,
#     positions = position,
#     positionThicknesses = thickness,
#     ids = id,
#     calCurves = calCurves,
#     predictPositions = TestRSLData$Depth
#   )
# )
# RSLrun <- with(
#   TestRSLData,
#   BchronRSL(RSLchron,
#     RSLmean = RSL,
#     RSLsd = Sigma,
#     degree = 3
#   )
# )

## ----eval=FALSE---------------------------------------------------------------
# summary(RSLrun, type = "RSL", age_grid = seq(0, 2000, by = 250))
# plot(RSLrun, type = "RSL") + ggtitle("Relative sea level plot")
# plot(RSLrun, type = "rate") + ggtitle("Rate of RSL change") +
#   ylab("Rate (mm per year)")

## ----results='hide', eval=FALSE-----------------------------------------------
# data(Sluggan)
# SlugDens <- BchronDensity(
#   ages = Sluggan$ages,
#   ageSds = Sluggan$ageSds,
#   calCurves = Sluggan$calCurves
# )

## ----eval=FALSE---------------------------------------------------------------
# summary(SlugDens, prob = 0.95)

## ----eval=FALSE---------------------------------------------------------------
# plot(SlugDens, xlab = "Age (cal years BP)", xaxp = c(0, 16000, 16))

## ----eval=FALSE---------------------------------------------------------------
# SlugDensFast <- BchronDensityFast(
#   ages = Sluggan$ages,
#   ageSds = Sluggan$ageSds,
#   calCurves = Sluggan$calCurves
# )

## ----eval = FALSE-------------------------------------------------------------
# # Load in the calibration curve with:
# intcal09 <- read.table(system.file("extdata/intcal09.14c",
#   package = "Bchron"
# ), sep = ",")
# 
# # Run createCalCurve
# createCalCurve(
#   name = "intcal09", calAges = intcal09[, 1],
#   uncalAges = intcal09[, 2], oneSigma = intcal09[, 3]
# )

## ----eval = FALSE-------------------------------------------------------------
# file.copy("intcal09.rda", system.file("data", package = "Bchron"))

## ----eval = FALSE-------------------------------------------------------------
# age_09 <- BchronCalibrate(
#   age = 15500, ageSds = 150, calCurves = "intcal09",
#   ids = "IntCal09"
# )
# age_20 <- BchronCalibrate(
#   age = 15500, ageSds = 150, calCurves = "intcal20",
#   ids = "Intcal20"
# )
# library(ggplot2)
# plot(age_09) +
#   geom_line(
#     data = as.data.frame(age_20$Intcal20),
#     aes(x = ageGrid, y = densities), col = "red"
#   ) +
#   ggtitle("Intcal09 vs Intcal20")

## ----message=FALSE------------------------------------------------------------
age_rc <- BchronCalibrate(
  age = 3000 - 80,
  ageSds = sqrt(50^2 + 30^2),
  calCurves = "marine20",
  ids = "Res_corr"
)
summary(age_rc)

## ----eval = FALSE-------------------------------------------------------------
# myChron <- with(
#   myDataSet,
#   Bchronology(
#     ages = ages + offsets,
#     ageSds = sqrt(ageSds + offsetSds),
#     calCurves = calCurves,
#     positions = positions,
#     positionThicknesses = thicknesses,
#     ids = ids
#   )
# )

## ----results = 'hide'---------------------------------------------------------
unCal1 <- unCalibrate(2350, type = "ages")

## -----------------------------------------------------------------------------
print(unCal1)

## ----results = 'hide'---------------------------------------------------------
unCal2 <- unCalibrate(
  calAge = c(2350, 4750, 11440),
  type = "ages",
  calCurve = "intcal20"
)

## -----------------------------------------------------------------------------
print(unCal2)

## ----results = 'hide'---------------------------------------------------------
ageRange <- seq(8000, 9000, by = 5)
c14Ages <- unCalibrate(ageRange,
  type = "ages"
)
load(system.file("data/intcal20.rda", package = "Bchron"))
ggplot(intcal20, aes(x = V1, y = V2)) +
  geom_line() +
  theme_bw() +
  scale_x_continuous(limits = range(ageRange)) +
  scale_y_continuous(limits = range(c14Ages)) +
  geom_rug(
    data = data.frame(c14Ages), aes(y = c14Ages),
    sides = "l", inherit.aes = F
  ) +
  geom_rug(
    data = data.frame(ageRange), aes(x = ageRange),
    sides = "b", inherit.aes = F
  ) +
  labs(x = "Calendar years BP", y = "14C years BP")

## -----------------------------------------------------------------------------
calAge <- BchronCalibrate(
  ages = 11255,
  ageSds = 25,
  calCurves = "intcal20"
)
calSampleAges <- sampleAges(calAge)

## ----results = 'hide'---------------------------------------------------------
unCal <- unCalibrate(calSampleAges,
  type = "samples"
)

## -----------------------------------------------------------------------------
print(unCal)

## ----results = 'hide'---------------------------------------------------------
data(Glendalough)
GlenOut <- with(
  Glendalough,
  Bchronology(
    ages = ages,
    ageSds = ageSds,
    calCurves = calCurves,
    positions = position,
    positionThicknesses = thickness,
    ids = id,
    predictPositions = seq(0, 1500, by = 10)
  )
)

## -----------------------------------------------------------------------------
library(ggplot2)
plot(GlenOut) +
  labs(
    title = "Glendalough",
    x = "Age (cal years BP)",
    y = "Depth (cm)"
  )

## -----------------------------------------------------------------------------
new_cal <- BchronCalibrate(
  ages = 7000,
  ageSds = 40,
  calCurve = "intcal20"
)

## -----------------------------------------------------------------------------
library(ggridges)
plot(GlenOut) +
  geom_ridgeline(
    data = as.data.frame(new_cal$Date1),
    aes(
      x = ageGrid,
      y = 600,
      height = densities * 10000, # Note the 10000 came from trial and error
      group = "New date",
    ),
    fill = "grey",
    colour = "black"
  ) +
  annotate("text", x = 9000, y = 570, label = "New date")

## ----results = 'hide'---------------------------------------------------------
GlenOut2 <- with(
  Glendalough[-2,],
  Bchronology(
    ages = ages,
    ageSds = ageSds,
    calCurves = calCurves,
    positions = position,
    positionThicknesses = thickness,
    ids = id,
    predictPositions = seq(0, 1500, by = 10)
  )
)

## -----------------------------------------------------------------------------
alpha <- 0.95
chronRange <- data.frame(
      chronLow = apply(GlenOut2$thetaPredict, 2, "quantile", probs = (1 - alpha) / 2),
      chronMed = apply(GlenOut2$thetaPredict, 2, "quantile", probs = 0.5),
      chronHigh = apply(GlenOut2$thetaPredict, 2, "quantile", probs = 1 - (1 - alpha) / 2),
      positions = GlenOut2$predictPositions
    )

## -----------------------------------------------------------------------------
ageGrid <- with(chronRange, seq(min(chronLow), max(chronHigh),
      length = nrow(chronRange)
    ))
chronRangeSwap <- data.frame(
      Age = ageGrid,
      positionLow = with(chronRange, approx(chronLow, positions,
        xout = ageGrid,
        rule = 2
      )$y),
      Position = with(chronRange, approx(chronMed, positions,
        xout = ageGrid,
        rule = 2
      )$y),
      positionHigh = with(chronRange, approx(chronHigh, positions,
        xout = ageGrid,
        rule = 2
      )$y),
      Date = "Bchron",
      densities = NA,
      height = NA
    )

## -----------------------------------------------------------------------------
plot(GlenOut, dateLabels = FALSE, chronTransparency = 0.3) + 
  geom_ribbon(
        data = chronRangeSwap,
        aes(
          x = Age,
          ymin = positionLow,
          ymax = positionHigh
        ),
        colour = "red",
        fill = "red",
        alpha = 0.3
      )

