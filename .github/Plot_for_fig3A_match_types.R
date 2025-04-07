# Define adduct information (adduct name and corresponding directory)
adduct_directories <- list(
  c("[M+Na]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[M+K]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[M+NH4]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[M+H-H2O]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[2M+H]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[M+Cl]-", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[M-H-H2O]-", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[2M-H]-", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20")
)

# Directory to save plots and output files
save_directory <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"
output_directory <- file.path(save_directory, "Shifted_or_Shifted_Match")
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

library(dplyr)
library(ggplot2)
library(tidyr)

# Prepare a list to store results for each adduct
all_match_categories <- list()

# Process each adduct
for (adduct_info in adduct_directories) {
  adduct <- adduct_info[1]
  directory <- adduct_info[2]
  
  # Construct filename
  filename <- file.path(directory, paste0(gsub("[\\[\\]+-]", "", adduct), "_pairs_df_noise5.csv"))
  
  # Read the pairs dataframe
  if (!file.exists(filename)) {
    cat("File not found:", filename, "\n")
    next
  }
  
  pairs_df <- read.csv(filename)
  
  # Ensure numeric columns are handled correctly
  pairs_df$total_aligned <- as.numeric(pairs_df$total_aligned)
  pairs_df$Direct_match <- as.numeric(pairs_df$Direct_match)
  pairs_df$Shifted_match <- as.numeric(pairs_df$Shifted_match)
  
  pairs_df$total_aligned[is.na(pairs_df$total_aligned)] <- 0
  pairs_df$Direct_match[is.na(pairs_df$Direct_match)] <- 0
  pairs_df$Shifted_match[is.na(pairs_df$Shifted_match)] <- 0
  
  # Categorize matches with the specified logic
  pairs_df$match_category <- ifelse(
    pairs_df$total_aligned == 0, 
    "No_match",
    ifelse(
      pairs_df$total_aligned != 0 & pairs_df$Direct_match != 0 & pairs_df$Shifted_match == 0, 
      "Only_Direct_match",
      ifelse(
        pairs_df$total_aligned != 0 & pairs_df$Direct_match == 0 & pairs_df$Shifted_match != 0, 
        "Only_Shifted_match",
        ifelse(
          pairs_df$total_aligned != 0 & pairs_df$Direct_match != 0 & pairs_df$Shifted_match != 0, 
          "Both_Direct_and_Shifted_match", 
          "No_match"
        )
      )
    )
  )
  
  # Calculate match category distribution
  match_category_counts <- table(pairs_df$match_category)
  match_category_percentages <- prop.table(match_category_counts) * 100
  
  # Create a DataFrame for this adduct
  adduct_match_df <- data.frame(
    adduct = adduct,
    match_category = names(match_category_counts),
    count = as.vector(match_category_counts),
    percentage = as.vector(match_category_percentages)
  )
  
  # Store results
  all_match_categories[[adduct]] <- adduct_match_df
}

# Combine all results into a single DataFrame
combined_match_df <- do.call(rbind, all_match_categories)

# Reshape the data for stacked percentage plot
plot_data <- combined_match_df %>%
  group_by(adduct) %>%
  mutate(
    percentage = percentage / 100,
    match_category = factor(match_category, 
                            levels = c("Only_Direct_match", 
                                       "Both_Direct_and_Shifted_match", 
                                       "Only_Shifted_match", 
                                       "No_match"))
  ) %>%
  complete(match_category, fill = list(percentage = 0)) %>%
  arrange(adduct, match_category)

# Ensure adduct order matches the input list
plot_data$adduct <- factor(plot_data$adduct, 
                           levels = sapply(adduct_directories, `[[`, 1))

# Create plot
# Create plot
plot <- ggplot(plot_data, aes(x = adduct, y = percentage, fill = match_category)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x axis text
    axis.text.y = element_blank(),  # Remove y axis text
    axis.title.x = element_blank(), # Remove x axis title
    axis.title.y = element_blank(), # Remove y axis title
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),  # Remove plot title
    legend.position = "none"  # Remove the legend
  ) +
  scale_fill_manual(
    values = c(
      "Only_Direct_match" = "#2667FF", 
      "Both_Direct_and_Shifted_match" = "purple", 
      "Only_Shifted_match" = "#E25D75", 
      "No_match" = "gray"
    ),
    labels = c(
      "Only Direct match", 
      "Both Direct and Shifted match", 
      "Only Shifted match", 
      "No match"
    ),
    breaks = c("Only_Direct_match", 
               "Both_Direct_and_Shifted_match", 
               "Only_Shifted_match", 
               "No_match")
  )


# Save the plot
plot_filename <- file.path(save_directory, "Adduct_Match_Categories_Stacked.tiff")
ggsave(plot_filename, plot, width = 4, height = 3, dpi = 200)

# Save the plot data to CSV
write.csv(plot_data, 
          file.path(save_directory, "Adduct_Match_Categories_Data.csv"), 
          row.names = FALSE)

# Print the plot data
print(plot_data)