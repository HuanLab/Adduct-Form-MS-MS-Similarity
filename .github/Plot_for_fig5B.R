# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

# Set working directory and read files
output_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start"
output_files <- list.files(output_dir, pattern = "_noise5\\.csv$", full.names = TRUE)

combined_data <- output_files |>
  lapply(function(f) {
    cat("Reading:", basename(f), "\n")
    read.csv(f, stringsAsFactors = FALSE)
  }) |>
  do.call(what = rbind)

# Compute weighted intensity for fragments with m/z < 100 and >1% relative intensity
compute_weighted_intensity <- function(mz_str, int_str) {
  mz_vec <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", mz_str), ",")))
  int_vec <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", int_str), ",")))
  if (length(mz_vec) != length(int_vec) || length(mz_vec) == 0) return(NA_real_)
  
  max_int <- max(int_vec, na.rm = TRUE)
  keep_idx <- which(int_vec >= 0.01 * max_int & mz_vec < 100)
  if (length(keep_idx) == 0) return(NA_real_)
  
  mz_final <- mz_vec[keep_idx]
  int_final <- int_vec[keep_idx]
  sum(mz_final * int_final) / sum(int_final)
}

# Apply weighted intensity calculation
combined_data$Weighted_Intensity1 <- mapply(compute_weighted_intensity, combined_data$whole_mz1, combined_data$whole_int1)
combined_data$Weighted_Intensity2 <- mapply(compute_weighted_intensity, combined_data$whole_mz2, combined_data$whole_int2)
combined_data$Weighted_Intensity_Diff <- (combined_data$Weighted_Intensity2 - combined_data$Weighted_Intensity1) / combined_data$Weighted_Intensity1

# Prepare CE bin and filter for complete data
ce_levels <- c("[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+", "[2M+H]+", "[M+Cl]-", "[M-H-H2O]-", "[2M-H]-")

combined_data <- combined_data |>
  filter(complete.cases(.)) |>
  mutate(
    Average_CE = round((CE1 + CE2) / 2),
    CE_bin = case_when(
      Average_CE < 20 ~ "Low",
      Average_CE < 40 ~ "Medium",
      Average_CE <= 60 ~ "High",
      TRUE ~ NA_character_
    ),
    Adduct2 = factor(Adduct2, levels = ce_levels)
  )

# Filter to only Metabolite + Adduct combinations with all 3 CE bins
valid_combinations <- combined_data |>
  filter(!is.na(CE_bin)) |>
  group_by(Metabolite_Name, Adduct2) |>
  summarise(Num_CE_Bins = n_distinct(CE_bin), .groups = "drop") |>
  filter(Num_CE_Bins == 3)

filtered_data <- combined_data |>
  inner_join(valid_combinations, by = c("Metabolite_Name", "Adduct2"))

# Summarize and reshape
intensity_summary_wide <- filtered_data |>
  group_by(Metabolite_Name, Adduct2, CE_bin) |>
  summarise(Avg_Weighted_Intensity_Diff = mean(Weighted_Intensity_Diff, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(
    names_from = CE_bin,
    values_from = Avg_Weighted_Intensity_Diff,
    names_prefix = "Avg_IntensityDiff_"
  )

# Prepare data for plotting
ce_bin_numeric <- c("Low" = 1, "Medium" = 2, "High" = 3)

signed_mean_line_data <- intensity_summary_wide |>
  pivot_longer(cols = starts_with("Avg_IntensityDiff_"), names_to = "CE_bin", values_to = "Mean_Diff") |>
  mutate(
    CE_bin = factor(gsub("Avg_IntensityDiff_", "", CE_bin), levels = names(ce_bin_numeric)),
    CE_numeric = ce_bin_numeric[as.character(CE_bin)]
  ) |>
  group_by(Adduct2, CE_bin, CE_numeric) |>
  summarise(Mean_Diff = round(mean(Mean_Diff, na.rm = TRUE), 2), .groups = "drop") |>
  mutate(
    Adduct2 = factor(Adduct2, levels = ce_levels),
    Point_Color = case_when(
      Mean_Diff > 0.02  ~ "#E25D75",
      Mean_Diff < -0.02 ~ "#2667FF",
      TRUE              ~ "gray70"
    )
  )

# Plot and export
p <- ggplot(signed_mean_line_data, aes(x = CE_numeric, y = Mean_Diff)) +
  geom_point(aes(color = Point_Color), size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2, aes(group = Adduct2)) +
  scale_color_identity() +
  scale_x_continuous(breaks = 1:3, labels = c("Low", "Medium", "High")) +
  facet_wrap(~ Adduct2, nrow = 2, ncol = 4) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(face = "bold"),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.6, "cm"),
    panel.spacing.y = unit(0.8, "cm")
  )

ggsave(
  filename = file.path(output_dir, "Weighted_mz/Mean_Diff_Plot.tiff"),
  plot = p,
  width = 5.78,
  height = 3,
  dpi = 200,
  units = "in",
  bg = "transparent",
  device = "tiff"
)
