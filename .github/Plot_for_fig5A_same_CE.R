# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define a polynomial function for curve fitting
polynomial_func <- function(x, a, b, c, d) {
  a * x^3 + b * x^2 + c * x + d
}

# Define adduct information (adduct name and corresponding directory)
adduct_directories <- list(
  c("[M+Na]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[M+K]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[M+NH4]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[M+H-H2O]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[2M+H]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[M+Na]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[M+Cl]-", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[M-H-H2O]-", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"),
  c("[2M-H]-", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20")
)

# Directory to save plots
save_directory <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"

# Ensure the save directory exists
if (!dir.exists(save_directory)) {
  dir.create(save_directory, recursive = TRUE)
}

# Process each adduct
for (adduct_info in adduct_directories) {
  adduct <- adduct_info[1]
  directory <- adduct_info[2]
  
  cat("\nProcessing adduct:", adduct, "from directory:", directory, "\n")
  
  # Construct the filename
  filename <- file.path(directory, paste0(adduct, "_pairs_df_noise5.csv"))
  
  # Check if the file exists
  if (!file.exists(filename)) {
    cat("File not found:", filename, "\n")
    next
  }
  
  # Read the CSV file
  pairs_df <- read.csv(filename)
  
  # Calculate the average CE and round it to the nearest integer
  pairs_df <- pairs_df %>%
    mutate(Average_CE = round((CE1 + CE2) / 2))
  
  # Filter for identical instrument type pairs
  pairs_df <- pairs_df %>%
    filter(INSTRUMENT_TYPE_1 == INSTRUMENT_TYPE_2)
  
  # Calculate the average remddp for each rounded Average_CE
  avg_scores <- pairs_df %>%
    group_by(Average_CE) %>%
    summarise(remddp = mean(remddp, na.rm = TRUE)) %>%
    ungroup()
  
  # Fit a polynomial curve to the data
  x_data <- avg_scores$Average_CE
  y_data <- avg_scores$remddp
  fit <- nls(y_data ~ polynomial_func(x_data, a, b, c, d), start = list(a = 1, b = 1, c = 1, d = 1))
  params <- coef(fit)
  
  # Generate smooth curve data
  x_smooth <- seq(min(x_data), max(x_data), length.out = 100)
  y_smooth <- polynomial_func(x_smooth, params["a"], params["b"], params["c"], params["d"])
  
  # Create the plot
  summary_plot <- ggplot(avg_scores, aes(x = Average_CE, y = remddp)) +
    geom_point(color = "black", size = 0.8) +  # Plot dots
    geom_line(data = data.frame(x = x_smooth, y = y_smooth), aes(x = x, y = y), color = "red", linewidth = 0.8) +  # Plot fitted curve
    theme_minimal() +
    theme(
      legend.position = "none",  # Remove legend
      axis.title = element_blank(),  # Remove axis titles
      axis.text = element_blank(),  # Remove axis labels
      axis.ticks = element_blank(),  # Remove axis ticks
      panel.grid = element_blank(),  # Remove grid
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)  # Add frame with thickness 0.3
    ) +
    xlim(-5, 65) +  # Set x-axis limits
    ylim(-5, 105)  # Set y-axis limits
  
  # Print the plot
  print(summary_plot)
  
  # Clean adduct name for file saving
  clean_adduct_name <- adduct
  
  # Save summary plot as TIFF with transparent background
  tiff_filename <- file.path(save_directory, paste0(clean_adduct_name, "_Average_remddp_vs_Average_CE.tiff"))
  ggsave(
    filename = tiff_filename,
    plot = summary_plot,
    device = "tiff",
    width = 1.5,
    height = 1.5,
    dpi = 150,
    bg = "transparent"  # Set background to transparent
  )
  
  cat("Plot saved as:", tiff_filename, "\n")
}
