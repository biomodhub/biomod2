## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(hypervolume)
library(palmerpenguins)
library(ggplot2)
library(gridExtra)
library(raster)

## ----results = "hide"---------------------------------------------------------
data(penguins)
data(quercus)

## ----results = "hide", eval = FALSE-------------------------------------------
# data("quercus")
# data_alba = subset(quercus, Species=="Quercus alba")[,c("Longitude","Latitude")]
# data_rubra = subset(quercus, Species=="Quercus rubra")[,c("Longitude","Latitude")]
# climatelayers <- getData('worldclim', var='bio', res=10, path=tempdir())
# 
# # z-transform climate layers to make axes comparable
# climatelayers_ss = climatelayers[[c(1,4,12,15)]]
# for (i in 1:nlayers(climatelayers_ss))
# {
#   climatelayers_ss[[i]] <- (climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / cellStats(climatelayers_ss[[i]], 'sd')
# }
# climatelayers_ss_cropped = crop(climatelayers_ss, extent(-150,-50,15,60))
# 
# # extract transformed climate values
# climate_alba = extract(climatelayers_ss_cropped, data_alba)
# climate_rubra = extract(climatelayers_ss_cropped, data_rubra)
# combined_sample = data.frame(rbind(climate_rubra, climate_alba))
# combined_sample["Species"] = quercus$Species
# 
# # Create hypervolumes
# hv_alba = hypervolume(climate_alba)
# hv_rubra = hypervolume(climate_rubra)

## ----eval=FALSE---------------------------------------------------------------
# alba_seq = hypervolume_resample("alba_seq", hv_alba, "bootstrap seq", n = 20, seq = seq(100, 1700, 400), cores = 32)
# rubra_seq = hypervolume_resample("rubra_seq", hv_rubra, "bootstrap seq", n = 20, seq = seq(100, 2100, 400), cores = 32)
# 
# # Funnel Plots
# alba_plot = hypervolume_funnel(alba_seq) +
#   geom_point(aes(y = upperq)) +
#   geom_point(aes(y = lowerq)) +
#   geom_point(aes(y = sample_mean), col = "blue") +
#   theme_bw() +
#   labs(title = "a)", subtitle = NULL) +
#   ylab("Volume") +
#   ylim(0, 0.8) +
#   xlim(0, 2100)
# rubra_plot = hypervolume_funnel(rubra_seq) +
#   geom_point(aes(y = upperq)) +
#   geom_point(aes(y = lowerq)) +
#   geom_point(aes(y = sample_mean), col = "blue") +
#   theme_bw() +
#   labs(title = "b)", subtitle = NULL) +
#   ylab("Volume") +
#   ylim(0, 0.8) +
#   xlim(0, 2100)
# grid.arrange(alba_plot, rubra_plot, nrow = 1)

## ----out.width="100%", echo=FALSE---------------------------------------------
url = "funnel_plot_comparisons.png"
knitr::include_graphics(url)

## ----eval = FALSE-------------------------------------------------------------
# data("quercus")
# data_alba = subset(quercus, Species=="Quercus alba")[,c("Longitude","Latitude")]
# climatelayers <- getData('worldclim', var='bio', res=10, path=tempdir())
# 
# # z-transform climate layers to make axes comparable
# climatelayers_ss = climatelayers[[c(1,4,12,15)]]
# for (i in 1:nlayers(climatelayers_ss))
# {
#   climatelayers_ss[[i]] <- (climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / cellStats(climatelayers_ss[[i]], 'sd')
# }
# climatelayers_ss_cropped = crop(climatelayers_ss, extent(-150,-50,15,60))
# 
# # extract transformed climate values
# climate_alba = extract(climatelayers_ss_cropped, data_alba)
# 
# # Generate Hypervolumes
# hv_alba = hypervolume(climate_alba,name='alba',samples.per.point=10)
# 
# # Give all samples to the east of -75 longitude double the weight of all other observations to compensate for undersampling
# weights = ifelse(data_alba$Longitude > -75, 2, 1)
# 
# # Create new hypervolume using weighted bootstrapping
# weighted_hv_alba = hypervolume_resample(hv = hv_alba, n=1, method = "weighted bootstrap", weights = weights, to_file = FALSE)
# weighted_hv_alba = weighted_hv_alba[[1]]
# 
# # Perform overlap test to see if the weights changed the distribution
# alba_biased_combined_data = rbind(hv_alba@Data, weighted_hv_alba@Data)
# alba_biased_combined_hv = hypervolume(alba_biased_combined_data)
# alba_biased_combined_samples = hypervolume_resample("combined_biased_resample", hv = alba_biased_combined_hv, method = "bootstrap", n = 128, points_per_resample = nrow(alba_biased_combined_data)/2, verbose = FALSE, cores = 32)
# result = hypervolume_overlap_test(hv_alba, weighted_hv_alba, alba_biased_combined_samples, cores = 32)
# 
# # Show null distribution and observed value
# result$plots$sorensen +
#   theme_bw() +
#   labs(title = NULL) +
#   xlab("Sorensen Distance") +
#   ylab("Density")

## ----out.width="100%", echo = FALSE-------------------------------------------
url = "biased_bootstrap_overlap.png"
knitr::include_graphics(url)

## ----eval=FALSE---------------------------------------------------------------
# hv = hypervolume(na.omit(penguins)[,3:6], verbose = FALSE)
# cols_to_weigh = c("bill_length_mm", "bill_depth_mm")
# 
# # Highest weights are assigned to max bill length and max bill depth
# mu = apply(hv@Data, 2, max)[cols_to_weigh]
# sigma = apply(hv@Data, 2, var)[cols_to_weigh]*2
# biased_path = hypervolume_resample("Bill bias", hv, method = "weighted bootstrap", n = 1, mu = mu, sigma = sigma, cols_to_weigh = cols_to_weigh)
# 
# # Read in hypervolume object from file
# biased_hv = readRDS(file.path(biased_path, "resample 1.rds"))
# 
# combined_dat = data.frame(rbind(hv@Data, biased_hv@Data))
# combined_dat['Type'] = rep(c('original', 'biased'), each = nrow(hv@Data))

## ----eval=FALSE---------------------------------------------------------------
# plot1 = ggplot(combined_dat) +
#   geom_histogram(aes(x = bill_depth_mm, fill = Type), bins = 20) +
#   facet_wrap(~Type) +
#   ggtitle("Distribution of Bill Depth") +
#   xlab("Bill depth (mm)") +
#   ylab("Count") +
#   theme_bw() +
#   theme(legend.position = "none")
# plot2 = ggplot(combined_dat) +
#   geom_histogram(aes(x = bill_length_mm, fill = Type), bins = 20) +
#   facet_wrap(~Type) +
#   ggtitle("Distribution of Bill Length") +
#   xlab("Bill length (mm)") +
#   ylab(NULL) +
#   theme_bw() +
#   theme(legend.position = "none")
# plot3 = ggplot(combined_dat) +
#   geom_histogram(aes(x = flipper_length_mm, fill = Type), bins = 20) +
#   facet_wrap(~Type) +
#   ggtitle("Distribution of Flipper Length") +
#   xlab("Flipper length (mm)") +
#   ylab("Count")+
#   theme_bw() +
#   theme(legend.position = "none")
# plot4 = ggplot(combined_dat) +
#   geom_histogram(aes(x = body_mass_g, fill = Type), bins = 20) +
#   facet_wrap(~Type) +
#   ggtitle("Distribution of Body Mass") +
#   xlab("Body mass (g)") +
#   ylab(NULL) +
#   theme_bw() +
#   theme(legend.position = "none")
# grid.arrange(plot1, plot2, plot3, plot4, nrow = 2)

## ----out.width="100%", echo = FALSE-------------------------------------------
url = "penguin_weights.png"
knitr::include_graphics(url)

## ----eval=FALSE---------------------------------------------------------------
# data("quercus")
# data_alba = subset(quercus, Species=="Quercus alba")[,c("Longitude","Latitude")]
# data_rubra = subset(quercus, Species=="Quercus rubra")[,c("Longitude","Latitude")]
# climatelayers <- getData('worldclim', var='bio', res=10, path=tempdir())
# 
# # z-transform climate layers to make axes comparable
# climatelayers_ss = climatelayers[[c(1,4,12,15)]]
# for (i in 1:nlayers(climatelayers_ss))
# {
#   climatelayers_ss[[i]] <- (climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / cellStats(climatelayers_ss[[i]], 'sd')
# }
# climatelayers_ss_cropped = crop(climatelayers_ss, extent(-150,-50,15,60))
# 
# # extract transformed climate values
# climate_alba = extract(climatelayers_ss_cropped, data_alba)
# climate_rubra = extract(climatelayers_ss_cropped, data_rubra)
# 
# # Generate Hypervolumes
# hv_alba = hypervolume(climate_alba,name='alba',samples.per.point=10)
# hv_rubra = hypervolume(climate_rubra,name='rubra',samples.per.point=10)
# 
# # Method 1: 2hr runtime with 12 threads
# combined_sample = rbind(climate_alba, climate_rubra)
# population_hat = hypervolume(combined_sample)
# 
# # Create bootstrapped hypervolumes of both sample sizes
# method1_path_size_1669 = hypervolume_resample("quercus_1669_boot", population_hat, "bootstrap", n = 100, cores = 12)
# method1_path_size_2110 = hypervolume_resample("quercus_2110_boot", population_hat, "bootstrap", n = 100, cores = 12)
# 
# 
# result1 = hypervolume_overlap_test(hv_alba, hv_rubra, c(method1_path_size_1669, method1_path_size_2110), cores = 12)
# 
# #Method 2: 9hr runtime with 12 threads
# method2_path = hypervolume_permute("rubra_alba_permutation", hv1, hv2, n = 1357, cores = 12)
# 
# result2 = hypervolume_overlap_test(hv1, hv2, method2_path, cores = 2)
# 
# # Graphical Results of null sorensen statistic
# plot1 = result1$plots$sorensen +
#     xlab("Sorensen distance") +
#     ylab("Density") +
#     ggtitle("a)") +
#     xlim(0.7, 1) +
#     theme_bw()
# plot1 = result2$plots$sorensen +
#     xlab("Sorensen distance") +
#     ylab("Density") +
#     ggtitle("b)") +
#     xlim(0.7, 1) +
#     theme_bw()
# grid.arrange(plot1, plot2, ncol=2)

## ----out.width="100%", echo = FALSE-------------------------------------------
url = "overlap_test_demos.png"
knitr::include_graphics(url)

## ----eval=FALSE---------------------------------------------------------------
# library(foreach)
# library(mvtnorm)
# library(doParallel)
# 
# # Choose number of threads to use for parallel computing
# nthreads = detectCores()
# 
# # Load required libraries in the environment of each thread and register cluster to parallel backend
# cl = makeCluster(nthreads)
# clusterEvalQ(cl, {
#   library(hypervolume)
#   library(mvtnorm)
# })
# registerDoParallel(cl)
# 
# # Set shift distance and number of test replications
# N = 40
# M = 20
# offset = # shift distance
# reps = # Number of replications of the test
# 
# # Each replication takes around 4hrs. Total run time equals reps/nthreads * 4
# # Returns a list of p values given by the overlap test
# pvals = foreach(i = 1:reps, .combine = rbind) %dopar% {
#   # Random data N(0, I) and offset mean
#   x1 = rmvnorm(M, rep(0, 5))
#   x2 = rmvnorm(N, rep(offset, 5))
#   hv1 = hypervolume(x1)
#   hv2 = hypervolume(x2)
#   hv_combined = hypervolume(rbind(x1, x2))
# 
#   path_M = hypervolume_resample(paste0("shift_", offset, "_M_", i), hv_combined, method = "bootstrap", n = 20, points_per_resample = M)
#   path_N = hypervolume_resample(paste0("shift_", offset, "_N_", i), hv_combined, method = "bootstrap", n = 20, points_per_resample = N)
#   result = hypervolume_overlap_test(hv1, hv2, c(path_M, path_N))
#   result$p_values
# }
# 
# # Power is calculated as the percent of test replications that result in rejections which depends on the rejection threshold or alpha value.
# # We will use alpha = 0.05
# power = mean(pvals[,"sorensen"] <= 0.05)

## ----eval = FALSE-------------------------------------------------------------
# power0 = mean(pvals0[,"sorensen"] <= 0.05)
# power1 = mean(pvals1[,"sorensen"] <= 0.05)
# power2 = mean(pvals2[,"sorensen"] <= 0.05)
# power3 = mean(pvals3[,"sorensen"] <= 0.05)
# power4 = mean(pvals4[,"sorensen"] <= 0.05)
# power5 = mean(pvals5[,"sorensen"] <= 0.05)
# 
# ggplot(data =NULL, aes(x = c(0, 0.4472136 0.8944272 1.3416408 1.7888544 2.2360680), y = c(power0, power1, power2, power3, power4, power5))) +
#   geom_line() +
#   theme_bw() +
#   ylab("Power") +
#   xlab("Shift distance")

## ----out.width="100%", echo = FALSE-------------------------------------------
url = "shift_powers.png"
knitr::include_graphics(url)

## ----fig.width=7, fig.height=7 , echo=FALSE-----------------------------------
set.seed(123)
library(mvtnorm)
alpha_squared = 1 - 0:5 * 1/6
data1 = data.frame(rmvnorm(1000, mean = c(0, 0)))
data2 = data.frame(rmvnorm(1000, mean = c(0, 0), sigma = diag(c(alpha_squared[2], 1/alpha_squared[2]))))
data3 = data.frame(rmvnorm(1000, mean = c(0, 0), sigma = diag(c(alpha_squared[3], 1/alpha_squared[3]))))
data4 = data.frame(rmvnorm(1000, mean = c(0, 0), sigma = diag(c(alpha_squared[4], 1/alpha_squared[4]))))
data5 = data.frame(rmvnorm(1000, mean = c(0, 0), sigma = diag(c(alpha_squared[5], 1/alpha_squared[5]))))
data6 = data.frame(rmvnorm(1000, mean = c(0, 0), sigma = diag(c(alpha_squared[6], 1/alpha_squared[6]))))
data = rbind(data1, data2, data3, data4, data5, data6)
data["Scale factor"] = rep(c("a = 1", "a = 0.913", "a = 0.816", "a = 0.707", "a = 0.577", "a = 0.408"), each = 1000)

ggplot(data = data, aes(x = X1, y = X2, color = `Scale factor`)) + 
  geom_point(alpha = 0.4) + 
  theme_bw() +
  coord_fixed(ratio = 1) +
  facet_wrap(~`Scale factor`) +
  theme(legend.position = "none")

## ----eval=FALSE---------------------------------------------------------------
# N = 40
# M = 20
# offset = # Scale factor
# reps = # Number of replications of the test
# 
# # Since the data is only 2 dimensions each replication only takes a few minutes
# pvals = foreach(i = 1:reps, .combine = rbind) %dopar% {
#   x1 = rmvnorm(M, rep(0, 2))
#   x2 = rmvnorm(N, rep(0, 2), diag(c(offset^2, 1/(offset^2))))
#   hv1 = hypervolume(x1)
#   hv2 = hypervolume(x2)
#   hv_combined = hypervolume(rbind(x1, x2))
# 
#   path_M = hypervolume_resample(paste0("scale_", offset, "_M_", i), hv_combined, method = "bootstrap", n = 20, points_per_resample = M)
#   path_N = hypervolume_resample(paste0("scale_", offset, "_N_", i), hv_combined, method = "bootstrap", n = 20, points_per_resample = N)
#   result = hypervolume_overlap_test(hv1, hv2, c(path_M, path_N))
#   result$p_values
# }

## ----eval=FALSE---------------------------------------------------------------
# power0 = mean(pvals0[,"sorensen"] <= 0.05)
# power1 = mean(pvals1[,"sorensen"] <= 0.05)
# power2 = mean(pvals2[,"sorensen"] <= 0.05)
# power3 = mean(pvals3[,"sorensen"] <= 0.05)
# power4 = mean(pvals4[,"sorensen"] <= 0.05)
# power5 = mean(pvals5[,"sorensen"] <= 0.05)
# 
# ggplot(data =NULL, aes(x = sqrt(1 - 0:5 * 1/6), y = c(power0, power1, power2, power3, power4, power5))) +
#   geom_line() +
#   theme_bw() +
#   ylab("Power") +
#   xlab("Scale factor") +
#   scale_x_reverse()

## ----out.width="100%", echo = FALSE-------------------------------------------
url = "shape_powers.png"
knitr::include_graphics(url)

