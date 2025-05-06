# Load necessary libraries
library(dplyr)
library(ggplot2)

# Define adduct info
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

# Save directory
save_directory <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"

# Function to create remddp plot
plot_remddp_density <- function(data, title, binwidth = 5) {
  data_column <- data$remddp
  data_column <- pmin(pmax(na.omit(data_column), 0), 100)
  breaks <- seq(0, 100, by = binwidth)
  hist_data <- hist(data_column, breaks = breaks, plot = FALSE)
  hist_df <- data.frame(mid = hist_data$mids,
                        density = hist_data$density / sum(hist_data$density) * 100)
  bar_y_max <- max(hist_df$density, na.rm = TRUE)
  density_data <- density(data_column)
  scale_factor <- bar_y_max / max(density_data$y)
  
  ggplot() +
    geom_bar(data = hist_df, aes(x = mid, y = density), stat = "identity",
             fill = "black", color = "black", alpha = 0.1, width = binwidth) +
    geom_line(aes(x = density_data$x, y = density_data$y * scale_factor),
              color = "black", size = 1) +
    scale_y_continuous(name = "Normalized Density (%)",
                       limits = c(0, bar_y_max * 1.0)) +
    xlim(0, 100) +
    labs(title = title, x = "Cosine Similarity") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10))
}

# Save plot
save_plot <- function(plot, adduct) {
  file_path <- file.path(save_directory, paste0(adduct, "_Identical_Instrument_Density_Plot_remddp.tiff"))
  ggsave(file_path, plot, device = "tiff", width = 3, height = 2, dpi = 100)
  cat("  Saved:", file_path, "\n")
}

# Main loop
for (adduct_info in adduct_directories) {
  adduct <- adduct_info[1]
  directory <- adduct_info[2]
  filename <- file.path(directory, paste0(adduct, "_pairs_df_noise5.csv"))
  
  cat("\nProcessing:", adduct, "\n")
  if (!file.exists(filename)) {
    cat("  File not found:", filename, "\n")
    next
  }
  
  pairs_df <- read.csv(filename)
  identical_data <- filter(pairs_df, INSTRUMENT_TYPE_1 == INSTRUMENT_TYPE_2)
  remddp_plot <- plot_remddp_density(identical_data,
                                     title = paste("Identical Instrument Pairs -", adduct, "\nremddp"))
  save_plot(remddp_plot, adduct)
}
