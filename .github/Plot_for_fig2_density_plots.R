# Load necessary libraries
library(dplyr)
library(ggplot2)

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

# Define a function to create density plots with normalized histogram
plot_density_with_normalized_histogram <- function(data, title, binwidth = 5) {
  # Filter for identical instrument type pairs
  filtered_data <- data %>%
    filter(INSTRUMENT_TYPE_1 == INSTRUMENT_TYPE_2)
  
  # Function to create density plot for a specific column
  create_single_density_plot <- function(data_column, column_name) {
    # Remove NA values and clip to 0-100 range
    data_column <- data_column[!is.na(data_column)]
    data_column <- pmin(pmax(data_column, 0), 100)
    
    # Breaks from 0 to 100
    breaks <- seq(0, 100, by = binwidth)
    
    mean_value <- mean(data_column)
    median_value <- median(data_column)
    
    # Construct title with mean and median
    updated_title <- paste0(title, " - ", column_name, 
                            "\nMean = ", round(mean_value, 1), 
                            ", Median = ", round(median_value, 1),
                            "\n(n = ", length(data_column), ")")
    
    df <- data.frame(value = data_column)
    
    # Compute histogram data with normalization
    hist_data <- hist(df$value, breaks = breaks, plot = FALSE)
    hist_df <- data.frame(
      mid = hist_data$mids,
      density = hist_data$density / sum(hist_data$density) * 100  # Normalize to ensure total area is 100
    )
    
    # Determine y-axis range for geom_bar
    bar_y_max <- max(hist_df$density, na.rm = TRUE)
    
    # Compute the density plot data and scale it to match the histogram's maximum height
    density_data <- density(data_column)
    scale_factor <- bar_y_max / max(density_data$y)
    
    # Create plot with histogram bars and normalized density line
    ggplot() +
      # Normalized histogram bars
      geom_bar(data = hist_df, 
               aes(x = mid, y = density), 
               stat = "identity", 
               fill = "black", 
               color = "black", 
               alpha = 0.1,
               width = binwidth) +
      # Normalized Density plot
      geom_line(aes(x = density_data$x, y = density_data$y * scale_factor),
                color = "black", size = 1) +
      scale_y_continuous(
        name = "Normalized Density (%)", 
        limits = c(0, bar_y_max * 1.0)
      ) +
      xlim(0, 100) +
      labs(title = updated_title,
           x = "Cosine Similarity") +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 10))
  }
  
  # Create plots for each variable
  plots <- list(
    dp = create_single_density_plot(filtered_data$dp, "dp"),
    remddp = create_single_density_plot(filtered_data$remddp, "remddp"),
    FP = create_single_density_plot(filtered_data$FP, "FP"),
    ts_remddp = create_single_density_plot(filtered_data$ts_remddp, "ts_remddp")
  )
  
  return(plots)
}

# Save plot function (remains the same as in previous script)
save_plot <- function(plot, adduct, variable_name) {
  file_name <- paste0(adduct, "_Identical_Instrument_Density_Plot_", variable_name, ".tiff")
  file_path <- file.path(save_directory, file_name)
  ggsave(file_path, plot, device = "tiff", width = 3, height = 2, dpi = 100)
  cat("  Saved:", file_path, "\n")
}

# Process each adduct (same as previous script)
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
  
  # Generate plots for identical instrument type pairs
  merged_plots <- plot_density_with_normalized_histogram(
    pairs_df, 
    title = paste("Identical Instrument Pairs -", adduct)
  )
  
  # Save plots for each metric
  for (type in c("dp", "remddp", "FP", "ts_remddp")) {
    save_plot(merged_plots[[type]], adduct, type)
  }
}