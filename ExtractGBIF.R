# Load necessary package
library(data.table)
library(ggplot2)

# Read the file with explicit settings
file_path <- "Triclisia_subcordata.csv"
df <- fread(file_path, sep = "\t", header = TRUE, fill = TRUE, quote = "")

# Convert 'eventDate' column to Date class
df$Date <- as.Date(df$eventDate, format = "%Y-%m-%d")

# Create a 'Year' column
df$Year <- format(df$Date, "%Y")  # Extract the year

# Create a 'JulianDay' column
df$JulianDay <- as.numeric(format(df$Date, "%j"))  # Extract day of the year

# View the first few rows of the updated data frame
head(df)

#Create a new table with the selected columns
new_table <- df[, c("species", "Year", "JulianDay")]

# Save the new table as a CSV file
write.csv(new_table, "Triclisia_subcordata_year_julianday.csv", row.names = FALSE)


# Plot species

ggplot(new_table, aes(x = Year, y = JulianDay)) +
  geom_point() +  
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "Julian Day vs Year",
       x = "Year",
       y = "Julian Day") +
  ylim(0, 356) +  # Set y-axis limits
  theme_minimal()


