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
output_directory <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/Shifted_or_Shifted_Match"
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

library(dplyr)
library(ggplot2)

# Prepare an empty data frame to store all processed data
all_data <- data.frame()

# Create a vector to maintain the original order of adducts
adduct_order <- sapply(adduct_directories, `[[`, 1)

# Process each adduct
for (adduct_info in adduct_directories) {
  adduct <- adduct_info[1]
  directory <- adduct_info[2]
  
  cat("\nProcessing adduct:", adduct, "from directory:", directory, "\n")
  
  filename <- file.path(directory, paste0(gsub("[\\[\\]+-]", "", adduct), "_pairs_df_noise5.csv"))
  
  if (!file.exists(filename)) {
    cat("File not found:", filename, "\n")
    next
  }
  
  pairs_df <- read.csv(filename)
  
  if (!all(c("Shifted_match", "Direct_match", "total_aligned") %in% names(pairs_df))) {
    cat("Missing necessary columns in:", filename, "\n")
    next
  }
  
  pairs_df <- pairs_df %>% 
    filter(total_aligned != 0) %>%
    arrange(desc(Shifted_match)) %>%
    mutate(
      ID = row_number() / n(), 
      Adduct = adduct
    )
  
  all_data <- bind_rows(all_data, pairs_df)
}

# Ensure the Adduct column is a factor with levels in the original order
all_data$Adduct <- factor(all_data$Adduct, levels = adduct_order)

# Plot with line plot and facet_wrap
plot <- ggplot(all_data) +
  geom_line(aes(x = ID, y = Shifted_match, color = 'Shifted Match'), linewidth = 1) +
  geom_line(aes(x = ID, y = Direct_match, color = 'Direct Match'), linewidth = 1) +
  facet_wrap(~ Adduct, ncol = 4) +
  labs(
    title = "Distribution of Matches by Adduct",
    x = "Normalized ID", 
    y = "Match Score", 
    color = "Match Type"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 5),  # Adjust facet label size if needed
    legend.position = "none"  # Remove legends
  ) +
  scale_color_manual(values = c('Shifted Match' = "#E25D75", 'Direct Match' = "#2667FF"))

# Save the plot
plot_filename <- file.path(save_directory, "Adduct_Match_Distribution.tiff")
ggsave(plot_filename, plot, width = 5.95, height = 3.3, dpi = 300)
cat("Plot saved to:", plot_filename, "\n")