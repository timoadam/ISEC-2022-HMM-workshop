## ----setup, include = FALSE-----------------------------------------------------------------------------------
knitr::opts_chunk$set(
  message = FALSE, error = FALSE, warning = FALSE,
  comment = NA
)


## ----set_directory, include = FALSE---------------------------------------------------------------------------
setwd("/Users/timoadam/Desktop/ISEC 2022 HMM workshop")
library(moveHMM)


## ----load_moveHMM, eval = FALSE-------------------------------------------------------------------------------
## ## Install and load moveHMM
## install.packages("moveHMM")
## library(moveHMM)


## ----data, eval = FALSE---------------------------------------------------------------------------------------
## ## Load data from Movebank
## URL <- paste0("https://www.datarepository.movebank.org/bitstream/handle/",
##               "10255/move.568/At-sea%20distribution%20Antarctic%20Petrel",
##               "%2c%20Antarctica%202012%20%28data%20from%20Descamps%20et%",
##               "20al.%202016%29-gps.csv")
## raw <- read.csv(url(URL))
## 
## ## Keep relevant columns: ID, time, lon, and lat
## data_all <- raw[, c(13, 3, 4, 5)]
## colnames(data_all) <- c("ID", "time", "lon", "lat")
## data_all$time <- as.POSIXct(data_all$time, tz = "MST")
## 
## ## Keep first 10 tracks for this example
## petrel_data <- subset(data_all, ID %in% unique(ID)[1:5])
## 
## ## Split tracks at gaps using the function split_at_gaps
## source(./R/split_at_gap.R)
## petrel_data <- split_at_gap(data = petrel_data, max_gap = 60)
## 
## ## Plot tracks
## ggplot(petrel_data, aes(lon, lat, col = ID)) + geom_path() +
##   geom_point(size = 0.3) + coord_map()
## 
## ## Create distance to centre covariate
## # Define centre for each track as first observation
## i0 <- c(1, which(petrel_data$ID[-1] != petrel_data$ID[-nrow(petrel_data)])
##         + 1)
## centres <- petrel_data[i0, c("ID", "lon", "lat")]
## petrel_data$centre_lon <- rep(centres$lon, rle(petrel_data$ID)$lengths)
## petrel_data$centre_lat <- rep(centres$lat, rle(petrel_data$ID)$lengths)
## 
## # Add distance to centre covariate (based on sp for great circle distance)
## petrel_data$d2c <- sapply(1:nrow(petrel_data), function(i) {
##   spDistsN1(pts = matrix(as.numeric(petrel_data[i, c("lon", "lat")]),
##             ncol = 2), pt = c(petrel_data$centre_lon[i],
##             petrel_data$centre_lat[i]), longlat = TRUE)
## })
## 
## # Divide by 1000 for numerical stability (i.e., unit = 1000 km)
## petrel_data$d2c <- petrel_data$d2c / 1000


## ----load_data------------------------------------------------------------------------------------------------
## Load the data
petrel_data <- read.csv("./data/petrel_data.csv")


## ----prep_data, cache = TRUE----------------------------------------------------------------------------------
## Compute step lengths and turning angles
hmm_data <- prepData(trackData = petrel_data, coordNames = c("lon", "lat"), 
                     type = "LL")


## ----show_data, eval = FALSE----------------------------------------------------------------------------------
## ## Get an overview of the data and to visualise the tracks
## head(hmm_data)
## plot(hmm_data)


## ----plot_histograms, eval = FALSE----------------------------------------------------------------------------
## ## Plot histograms of step length and turning angle
## par(mfrow = c(1, 2))
## hist(hmm_data$step)
## hist(hmm_data$angle)


## ----initial_values-------------------------------------------------------------------------------------------
## Initial values for the parameters of the state-dependent distributions
# For step length
step_mean_par0 <- c(2, 8, 16) # means
step_SD_par0 <- c(4, 5, 6) # SDs
step_zeroprob_par0 <- c(0.01, 0.01, 0.01) # zero probabilities
step_par0 <- c(step_mean_par0, step_SD_par0, step_zeroprob_par0)
# For turning angle
angle_mean_par0 <- c(0, 0, 0) # means
angle_concentration_par0 <- c(0.5, 0.7, 0.9) # concentrations
angle_par0 <- c(angle_mean_par0, angle_concentration_par0)


## ----model_fitting, message = FALSE, results='hide', cache = TRUE---------------------------------------------
## Fit 3-state HMM
mod <- fitHMM(data = hmm_data, nbStates = 3, stepDist = "gamma", 
                 angleDist = "wrpcauchy", stepPar0 = step_par0, 
                 anglePar0 = angle_par0, verbose = 1)


## ----print_model, eval = TRUE---------------------------------------------------------------------------------
## Print estimated model parameters
mod


## ----plot_model, message = FALSE, results='hide', fig.width = 10, fig.height = 5, out.width="100%", fig.align="center"----
## Plot the fitted state-dependent distributions
par(mfrow = c(1, 2))
plot(mod, plotTracks = FALSE, ask = FALSE)


## ----state_decoding, eval = FALSE-----------------------------------------------------------------------------
## ## Global state decoding
## viterbi(mod)
## ## Local state decoding
## stateProbs(mod)
## ## Show both
## plotStates(mod)


## ----add_covariates, message = FALSE, results = 'hide', cache = TRUE------------------------------------------
# Fit 3-state model with linear effect of d2c
mod_d2c <- fitHMM(data = hmm_data, nbStates = 3, stepDist = "gamma", 
                  angleDist = "wrpcauchy", stepPar0 = step_par0, 
                  anglePar0 = angle_par0, formula = ~ d2c, verbose = 1)

# Fit 3-state model with squared effect of d2c
mod_d2c2 <- fitHMM(data = hmm_data, nbStates = 3, stepDist = "gamma", 
                   angleDist = "wrpcauchy", stepPar0 = step_par0, 
                   anglePar0 = angle_par0, formula = ~ d2c + I(d2c ^ 2), 
                   verbose = 1)

# Fit 3-state model cubic effect of d2c
mod_d2c3 <- fitHMM(data = hmm_data, nbStates = 3, stepDist = "gamma", 
                   angleDist = "wrpcauchy", stepPar0 = step_par0, 
                   anglePar0 = angle_par0, formula = ~ d2c + I(d2c ^ 2) 
                   + I(d2c ^ 3), verbose = 1)


## ----stationary_d2c, message = FALSE, results = 'hide', cache = TRUE, fig.width = 6, fig.height = 6, out.width = "50%", fig.align = "center"----
## Plot stationary state probabilities as function of d2c
plotStationary(mod_d2c, plotCI = TRUE)


## ----stationary_d2c2, message = FALSE, results = 'hide', cache = TRUE, fig.width = 6, fig.height = 6, out.width = "50%", fig.align = "center"----
## Plot stationary state probabilities as function of d2c
plotStationary(mod_d2c2, plotCI = TRUE)


## ----stationary_d2c3, message = FALSE, results = 'hide', cache = TRUE, fig.width = 6, fig.height = 6, out.width = "50%", fig.align = "center"----
## Plot stationary state probabilities as function of d2c
plotStationary(mod_d2c3, plotCI = TRUE)


## ----model_selection------------------------------------------------------------------------------------------
## Compute AIC
AIC(mod, mod_d2c, mod_d2c2, mod_d2c3)


## ----model_checking, message = FALSE, results = 'hide', cache = TRUE, fig.width = 10, fig.height = 7.5, out.width = "100%", fig.align = "center"----
## Plot QQ and ACF plots
plotPR(mod_d2c2)


## ----load-packages--------------------------------------------------------------------------------------------
# # Detach moveHMM if loaded
# detach("package:moveHMM", unload = TRUE)
library(momentuHMM)
library(ggplot2)
theme_set(theme_bw())
pal <- c("#E69F00", "#56B4E9", "#009E73")


## ----load-data------------------------------------------------------------------------------------------------
# Load whitetip shark accelerometer data
data <- read.csv("./data/whitetip_data.csv")
data$Time <- as.POSIXct(data$Time)
head(data)
nrow(data)


## ----log-odba, fig.width = 4, fig.height = 3, out.width="50%", fig.align="center", fig.show='hold'------------
# Get logarithm of ODBA for HMM analysis
data$logODBA <- log(data$ODBA)

# Time series plots of ODBA and log-ODBA
ggplot(data, aes(Time, ODBA)) + geom_line()
ggplot(data, aes(Time, logODBA)) + geom_line()


## ----prep-data------------------------------------------------------------------------------------------------
data_hmm <- prepData(data = data, coordNames = NULL, covNames = "TimeOfDay")


## ----2state-init----------------------------------------------------------------------------------------------
# List of observation distributions
dist <- list(logODBA = "norm")

# List of initial parameters
mean0 <- c(-4, -2)
sd0 <- c(0.5, 0.5)
Par0_2s <- list(logODBA = c(mean0, sd0))


## ----2state-fit-----------------------------------------------------------------------------------------------
hmm_2s <- fitHMM(data = data_hmm, nbStates = 2, dist = dist, Par0 = Par0_2s)


## ----2state-plot, fig.width = 6, fig.height = 4, out.width="80%", fig.align="center"--------------------------
plot(hmm_2s, ask = FALSE)


## ----2state-viterbi, fig.width = 5, fig.height = 3, out.width="70%", fig.align="center"-----------------------
data_hmm$viterbi_2s <- factor(viterbi(hmm_2s))

ggplot(tail(data_hmm, 1000), aes(Time, logODBA, col = viterbi_2s, group = NA)) +
  geom_line() +
  scale_color_manual(values = pal, name = "State")


## ----2state-pseudores, fig.width = 7.5, fig.height = 2.5, out.width="100%", fig.align="center", fig.show='hold'----
# Plots of pseudo-residuals
plotPR(hmm_2s)


## ----3state-fit, fig.width = 4, fig.height = 3, out.width="50%", fig.align="center", fig.show='hold'----------
# List of initial parameters
mean0 <- c(-4, -3, -2)
sd0 <- c(0.5, 0.5, 0.5)
Par0_3s <- list(logODBA = c(mean0, sd0))

# Fit HMM
hmm_3s <- fitHMM(data = data_hmm, nbStates = 3, dist = dist, Par0 = Par0_3s)

# Plot state-dependent distributions
plot(hmm_3s, ask = FALSE)

# Get most likely state sequence
data_hmm$viterbi_3s <- factor(viterbi(hmm_3s))

# Plot log-ODBA coloured by states
ggplot(tail(data_hmm, 1000), aes(Time, logODBA, col = viterbi_3s, group = NA)) +
  geom_line() +
  scale_color_manual(values = pal, name = "State")


## ----3state-pseudores, fig.width = 7.5, fig.height = 2.5, out.width="100%", fig.align="center"----------------
plotPR(hmm_3s)


## ----3state-tod-fit, fig.width = 7, fig.height = 5, out.width="90%", fig.align="center"-----------------------
# Fit model with time of day covariate
hmm_3s_tod <- fitHMM(data = data_hmm, nbStates = 3, dist = dist, Par0 = Par0_3s,
                     formula = ~ cosinor(TimeOfDay, period = 24))

# Plot stationary state probabilities as function of time of day
plotStationary(hmm_3s_tod, plotCI = TRUE)

