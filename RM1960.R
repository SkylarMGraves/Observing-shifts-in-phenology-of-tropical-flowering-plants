library(circglmbayes)
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

traceplots_list <- list()
finalplot_list <- list()
coef_list <- data.frame()
preds_list <- list()
plant_list <- list()
sdat_list <- list()
slopes_intercepts <- data.frame()

df <- read.csv("Vangueriella_nigerica_year_julianday.csv")
df <- df %>% filter(Year >= 1960)
df$species <- factor(df$species)
index <- df %>% group_by(species) %>%
  summarise(n=n()) %>% filter(n>15)
df <- df %>% filter(species %in% index$species) %>% mutate(JulianDay = as.numeric(JulianDay))

spec <- unique(as.character(df$species))

for(j in 1:length(spec)){
  plant2 <- na.omit(df[df$species == spec[j],])
  Year <- plant2$Year
  Year_s <- scale(Year) #Year scaled and centered
  Year_center <- attr(Year_s, "scaled:center")
  Year_scale <- attr(Year_s, "scaled:scale")
  plant_list[[j]] <- plant2
  
  circ <- day2circ(plant2$JulianDay)
  sdat <- data.frame(Year, Year_s, circ)
  sdat_list[[j]] <- sdat 
  nchains <- 4
  chains <- list()
  for (k in 1:nchains ) {
    fit <- circGLM(circ ~ Year_s, data=sdat)
    chains[[k]] <- fit$all_chains
  }
  chains <- mcmc.list(chains)
  traceplots_list[[j]] <- chains
  
  fit <- circGLM(circ ~ Year_s, data=sdat, burnin=200, thin=30, Q=2500)
  samples <- fit$all_chains
  samplesdf <- data.frame(samples)
  names(samplesdf) <- names(samplesdf) %>% 
    gsub("_chain", "", .)
  
  Year <- seq(from=min(sdat$Year), to=max(sdat$Year), by=1) 
  Year_s <- scale(Year, center=Year_center, scale=Year_scale)
  n <- length(Year)
  results <- matrix(NA, nrow=n, ncol=5) 
  colnames(results) <- c("mnmu","mulo95","muhi95","ppdlo95","ppdhi95")
  
  for ( i in 1:n ) {
    mu <- samplesdf$b0 + inv_tan_half(samplesdf$bt * Year_s[i])
    results[i,1] <- mean(mu)
    results[i,2:3] <- quantile(mu, prob=c(0.025,0.975)) 
  }
  results <- circ2day(results) 
  preds <- data.frame(Year, Year_s, results)
  preds_list[[spec[j]]] <- preds
  
  coef_spec <- as.data.frame(coef(fit))
  coef_spec$Species <- spec[j]
  coef_spec$Effects <- rownames(coef_spec)
  coef_list <- rbind(coef_list,coef_spec)
  
  # Extract slope, intercept, and endpoint
  slope <- mean(samplesdf$bt)
  intercept <- mean(samplesdf$b0)
  endpoint <- max(preds$results[,1])
  slopes_intercepts <- rbind(slopes_intercepts, data.frame(Species=spec[j], Slope=slope, Intercept=intercept, Endpoint=endpoint))
}





############################

for (j in 1:length(spec)) {
  preds <- preds_list[[spec[j]]]
  
  # Pull date range
  min_year <- min(preds$Year)
  max_year <- max(preds$Year)
  
  print(paste("Species:", spec[j]))
  print(paste("Min Year:", min_year))
  print(paste("Max Year:", max_year))
  
  # Pull the start and end date of the regression line
  pred_min_year <- preds$mnmu[preds$Year == min_year]
  pred_max_year <- preds$mnmu[preds$Year == max_year]
  
  print(paste("Predicted Min Year JulianDay:", pred_min_year))
  print(paste("Predicted Max Year JulianDay:", pred_max_year))
  
  # Pull Slope
  slope_value <- slopes_intercepts$Slope[slopes_intercepts$Species == spec[j]]
  print(paste("Slope:", slope_value))
  
  slope_sd <- sd(samplesdf$bt)
  print(paste("SD:", slope_sd))
  
  
  # Create a figure with Year on x-axis and JulianDay on y-axis
  p <- ggplot(plant_list[[j]], aes(x = Year, y = JulianDay)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = paste("Year vs JulianDay for", spec[j]), x = "Year", y = "JulianDay") +
    ylim(1, 365)
  
  print(p)
}

###########################




# Create a figure with Year on x-axis and JulianDay on y-axis
p <- ggplot(df, aes(x = Year, y = JulianDay)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = " Stylosanthes_erecta Circular Distrobution", x = "Year", y = "JulianDay") +
  ylim(1, 365)
ggsave("Stylosanthes_erecta_Circular.png", plot = p, width = 8, height = 6, dpi = 300)
