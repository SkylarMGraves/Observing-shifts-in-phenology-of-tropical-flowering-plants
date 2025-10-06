### Linear (Non circular) ###


# Load required packages
library(brms)
library(dplyr)
library(ggplot2)

data <- read.csv("ceiba_schottii_year_julianday.csv")




# Function to fit Bayesian GLM for each species
fit_bayesian_glm <- function(data) {
  # Fit a Bayesian GLM (assuming normal distribution for dayOfYear)
  model <- brm(
    formula = JulianDay ~ Year,  # Change in dayOfYear on Year
    data = data,                 # Data for the species
    family = gaussian(),         # Assuming normal distribution for dayOfYear
    prior = c(
      prior(normal(0, 5), class = "b"),    # Prior for the coefficient of 'Year'
      prior(normal(0, 5), class = "Intercept")  # Prior for the intercept
    ),
    chains = 4,                  # Number of Markov chains
    iter = 2000,                 # Number of iterations per chain
    warmup = 1000,               # Number of warmup iterations
    control = list(adapt_delta = 0.95)  # Control parameter for sampling
  )
  
  # Return the model summary as a data frame
  summary_model <- summary(model)
  fixed_effects <- as.data.frame(summary_model$fixed)
  return(fixed_effects)
}

# Loop through each species and fit the model
species_models <- data %>%
  group_by(species) %>%
  group_modify(~ {
    # Fit the model for each species and return the model summary as a data frame
    species_data <- .x
    model_summary <- fit_bayesian_glm(species_data)
    return(model_summary)
  })

# View results for the first species
print(species_models)




# Assuming 'data' contains 'Year' and 'JulianDay'
# Extract model estimates
beta_0 <- species_models$Estimate[1]  # Intercept
beta_1 <- species_models$Estimate[2]  # Slope
ci_lower <- species_models$`l-95% CI`[2]  # Lower CI for slope
ci_upper <- species_models$`u-95% CI`[2]  # Upper CI for slope

# Generate predicted values
data <- data %>%
  mutate(predicted_JulianDay = beta_0 + beta_1 * Year,
         lower_CI = beta_0 + ci_lower * Year,
         upper_CI = beta_0 + ci_upper * Year)

# Pull date range
data <- na.omit(data)

min_year <- min(data$Year)
max_year <- max(data$Year)

print(min_year)
print(max_year)

# Pull the start and end date of the regression line

pred_min_year <- data$predicted_JulianDay[data$Year == min_year]
pred_max_year <- data$predicted_JulianDay[data$Year == max_year]

print(pred_min_year)
print(pred_max_year)

# Pull Slope

print(beta_1)

# Plot observed data and model predictions
p <- ggplot(data, aes(x = Year, y = JulianDay)) +
  geom_point(alpha = 0.5, color = "black") +  # Individual data points
  geom_line(aes(y = predicted_JulianDay), color = "blue", size = 1) +  # Regression line
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2, fill = "blue") +  # Credible interval
  theme_minimal() +
  labs(x = "Year", y = "Julian Day", title = " ceiba_schottii Linear ")

print(p)

ggsave("ceiba_schottii_Linear.png", plot = p, width = 8, height = 6, dpi = 300)







