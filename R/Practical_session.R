## 1. Analysing tracking data with moveHMM

## Install and load moveHMM
install.packages(moveHMM)
library(moveHMM)

## 1.1 Data pre-processing

## Load data
petrel_data <- read.csv("./data/petrel_data.csv")
head(petrel_data)
dim(petrel_data)
unique(petrel_data$ID)

## Compute step lengths and turning angles
hmm_data <- prepData(trackData = petrel_data, coordNames = c("lon", "lat"), type = "LL")
head(hmm_data)
plot(hmm_data)

## 1.2 Model specification/1.3 Model fitting

## Plot histograms of step length and turning angle
par(mfrow = c(1, 2))
hist(hmm_data$step)
hist(hmm_data$angle)
summary(hmm_data$step) # there are step lengths that are equal to zero

## Initial values for the parameters of the state-dependent distributions

## Examples of poor choices:
# For step length
step_mean_par0 <- c(2, 8, 80) # means
step_SD_par0 <- c(1, 4, 1) # SDs
step_zeroprob_par0 <- c(0.01, 0.01, 0.01) # zero probabilities
step_par0 <- c(step_mean_par0, step_SD_par0, step_zeroprob_par0)
# For turning angle
angle_mean_par0 <- c(0, 0, 0) # means
angle_concentration_par0 <- c(0.5, 0.7, 0.9) # concentrations
angle_par0 <- c(angle_mean_par0, angle_concentration_par0)

## Fit HMM
mod <- fitHMM(data = hmm_data, nbStates = 3, stepDist = "gamma", angleDist = "wrpcauchy", 
              stepPar0 = step_par0, anglePar0 = angle_par0, verbose = 2)
mod
mod$mod$minimum # log-likelihood: 13170.3

## Examples of poor choices:
# For step length
step_mean_par0 <- c(2, 8, 16) # means
step_SD_par0 <- c(1, 1, 1) # SDs
step_zeroprob_par0 <- c(0.01, 0.01, 0.01) # zero probabilities
step_par0 <- c(step_mean_par0, step_SD_par0, step_zeroprob_par0)
# For turning angle
angle_mean_par0 <- c(0, 0, 0) # means
angle_concentration_par0 <- c(0.5, 0.8, 0.8) # concentrations
angle_par0 <- c(angle_mean_par0, angle_concentration_par0)

## Fit HMM
mod <- fitHMM(data = hmm_data, nbStates = 3, stepDist = "gamma", angleDist = "wrpcauchy", 
              stepPar0 = step_par0, anglePar0 = angle_par0, verbose = 2)
mod
mod$mod$minimum # log-likelihood: 13170.3

## Good choice:
# For step length
step_mean_par0 <- c(2, 8, 16) # means
step_SD_par0 <- c(1, 4, 8) # SDs
step_zeroprob_par0 <- c(0.01, 0.01, 0.01) # zero probabilities
step_par0 <- c(step_mean_par0, step_SD_par0, step_zeroprob_par0)
# For turning angle
angle_mean_par0 <- c(0, 0, 0) # means
angle_concentration_par0 <- c(0.5, 0.7, 0.9) # concentrations
angle_par0 <- c(angle_mean_par0, angle_concentration_par0)

## Fit HMM
mod <- fitHMM(data = hmm_data, nbStates = 3, stepDist = "gamma", angleDist = "wrpcauchy", 
              stepPar0 = step_par0, anglePar0 = angle_par0, verbose = 2)
mod
mod$mod$minimum # log-likelihood: 10635.4

## 1.4 Results

## Plot the fitted state-dependent distributions
par(mfrow = c(1, 2))
plot(mod, plotTracks = FALSE, ask = FALSE)

## Plot the decoded tracks
plot(mod)

## State decoding
# Global decoding
viterbi(mod)
# Local decoding
stateProbs(mod)
# Show both
plotStates(mod)

## 1.5 Adding covariates on the state transition probabilities

## Fit HMM with linear effect of d2c
mod_d2c <- fitHMM(data = hmm_data, nbStates = 3, stepDist = "gamma", 
                  angleDist = "wrpcauchy", stepPar0 = step_par0, 
                  anglePar0 = angle_par0, formula = ~ d2c, verbose = 1)

## Fit HMM with squared effect of d2c
mod_d2c2 <- fitHMM(data = hmm_data, nbStates = 3, stepDist = "gamma", 
                   angleDist = "wrpcauchy", stepPar0 = step_par0, 
                   anglePar0 = angle_par0, formula = ~ d2c + I(d2c ^ 2), 
                   verbose = 1)

## Fit HMM cubic effect of d2c
mod_d2c3 <- fitHMM(data = hmm_data, nbStates = 3, stepDist = "gamma", 
                   angleDist = "wrpcauchy", stepPar0 = step_par0, 
                   anglePar0 = angle_par0, formula = ~ d2c + I(d2c ^ 2) 
                   + I(d2c ^ 3), verbose = 1)

## Plot stationary state probabilities as function of d2c
plotStationary(mod_d2c, plotCI = TRUE)
plotStationary(mod_d2c2, plotCI = TRUE)
plotStationary(mod_d2c3, plotCI = TRUE)

## 1.6 Model selection

## Compute AIC
AIC(mod, mod_d2c, mod_d2c2, mod_d2c3)

## 1.7 Model checking

## Plot QQ and ACF plots
plotPR(mod_d2c2)

## 2. Analysing accelerometer data with momentuHMM

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

