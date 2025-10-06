#Permutation Null Model

library(coda) #for mcmc tools to use with circglmbayes
library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)

# The link function for a circular GLM is the tan half function:
inv_tan_half <- function(x) {
  return(2 * atan(x))
}

# Converting day to circular
day2circ <- function(day_of_Year) {
  return(2 * pi * day_of_Year / 365 - pi)
}

circ2day <- function(circ) {
  return(365 * (circ + pi) / (2 * pi))
}


df <- read.csv("All_Species.csv")
df$species <- factor(df$species)
index <- df %>% group_by(species) %>%
  summarise(n=n()) %>% filter(n>15)
df <- df %>% filter(species %in% index$species) %>% mutate(JulianDay = as.numeric(JulianDay))

spec <- unique(as.character(df$species))




run_circular_null_model <- function(plant_df, n_iter = 1000) {
  Year <- plant_df$Year
  JulianDay <- as.numeric(plant_df$JulianDay)
  circ <- day2circ(JulianDay)
  
  # Center/scale the year once
  Year_s <- scale(Year)
  Year_center <- attr(Year_s, "scaled:center")
  Year_scale <- attr(Year_s, "scaled:scale")
  sdat <- data.frame(Year, Year_s, circ)
  
  # Fit the model to get observed slope
  fit_obs <- circGLM(circ ~ Year_s, data = sdat, burnin = 200, thin = 30, Q = 2500)
  samples_obs <- data.frame(fit_obs$all_chains)
  names(samples_obs) <- gsub("_chain", "", names(samples_obs))
  observed_slope <- mean(samples_obs$bt)
  
  # Run permutation null model
  null_slopes <- numeric(n_iter)
  for (i in 1:n_iter) {
    shuffled_year <- sample(Year)
    Year_s_perm <- scale(shuffled_year, center = Year_center, scale = Year_scale)
    sdat_perm <- data.frame(Year = shuffled_year, Year_s = Year_s_perm, circ = circ)
    
    fit_null <- circGLM(circ ~ Year_s, data = sdat_perm, burnin = 200, thin = 30, Q = 2500)
    samples_null <- data.frame(fit_null$all_chains)
    names(samples_null) <- gsub("_chain", "", names(samples_null))
    null_slopes[i] <- mean(samples_null$bt)
  }
  
  p_value <- mean(abs(null_slopes) >= abs(observed_slope))
  list(observed_slope = observed_slope, null_slopes = null_slopes, p = p_value)
}
