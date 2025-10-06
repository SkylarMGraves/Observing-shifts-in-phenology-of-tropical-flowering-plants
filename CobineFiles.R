library(dplyr)
library(readr)

# Set the path to your folder with CSV files
file_path <- "C:/Users/skly5321/OneDrive - UCB-O365/Attachments/FlowerOnceAYear/gbif downloads"

# List all CSV files in the folder
csv_files <- list.files(path = file_path, pattern = "\\_year_julianday.csv$", full.names = TRUE)

# Read and bind all CSVs together
combined_df <- csv_files %>%
  lapply(read_csv) %>%      # Or use read.csv if you prefer base R
  bind_rows()

# Optionally, write to a new CSV file
write_csv(combined_df, "All_Species.csv")
