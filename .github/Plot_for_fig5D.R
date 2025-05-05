# Set working directories
output_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start"

# List all output CSV files in the output directory
output_files <- list.files(output_dir, pattern = "_noise5\\.csv$", full.names = TRUE)

# Read all processed CSV files into a list of data frames
processed_data_list <- lapply(output_files, function(f) {
  cat("Reading:", basename(f), "\n")
  read.csv(f, stringsAsFactors = FALSE)
})

# Optionally combine into one big data frame (if same structure)
combined_data <- do.call(rbind, processed_data_list)

compute_weighted_intensity <- function(mz_str, int_str, smaller_mz) {
  # Convert string to numeric vectors
  mz_vec <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", mz_str), ",")))
  int_vec <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", int_str), ",")))
  
  if (length(mz_vec) != length(int_vec) || length(mz_vec) == 0) {
    return(NA_real_)
  }
  
  # Apply 1% intensity filter
  max_int <- max(int_vec, na.rm = TRUE)
  threshold <- 0.01 * max_int
  keep_idx <- which(int_vec >= threshold)
  
  if (length(keep_idx) == 0) {
    return(NA_real_)
  }
  
  mz_filtered <- mz_vec[keep_idx]
  int_filtered <- int_vec[keep_idx]
  
  # Further filter: exclude fragments with m/z > smaller_mz - 0.02
  final_idx <- which(mz_filtered <= smaller_mz - 0.02)
  if (length(final_idx) == 0) {
    return(NA_real_)
  }
  
  mz_final <- mz_filtered[final_idx]
  int_final <- int_filtered[final_idx]
  
  # Calculate weighted intensity
  sum(mz_final * int_final) / sum(int_final)
}

# Get the smaller precursor m/z value per row
combined_data <- combined_data %>%
  mutate(
    smaller_mz = pmin(precursor1, precursor2, na.rm = TRUE)
  )

# Recompute weighted intensities
combined_data$Weighted_Intensity1 <- mapply(
  compute_weighted_intensity,
  combined_data$whole_mz1,
  combined_data$whole_int1,
  combined_data$smaller_mz
)

combined_data$Weighted_Intensity2 <- mapply(
  compute_weighted_intensity,
  combined_data$whole_mz2,
  combined_data$whole_int2,
  combined_data$smaller_mz
)

# Recalculate the difference
combined_data$Weighted_Intensity_Diff <- (combined_data$Weighted_Intensity2 - combined_data$Weighted_Intensity1) / combined_data$Weighted_Intensity1


library(dplyr)
library(ggplot2)

# Step 1: Calculate average CE
combined_data <- combined_data %>%
  mutate(
    Average_CE = round((CE1 + CE2) / 2),
    Adduct2 = factor(Adduct2, levels = c(
      "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+",
      "[2M+H]+", "[M+Cl]-", "[M-H-H2O]-", "[2M-H]-"
    ))
  )


library(dplyr)

combined_data <- combined_data[complete.cases(combined_data), ]

# Step 1: Define CE interval membership for each row
combined_data <- combined_data %>%
  mutate(CE_bin = case_when(
    Average_CE >= 0 & Average_CE < 20 ~ "Low",
    Average_CE >= 20 & Average_CE < 40 ~ "Medium",
    Average_CE >= 40 & Average_CE <= 60 ~ "High",
    TRUE ~ NA_character_
  ))

# Step 2: Count unique CE bins per Metabolite_Name + Adduct2
valid_combinations <- combined_data %>%
  filter(!is.na(CE_bin)) %>%
  group_by(Metabolite_Name, Adduct2) %>%
  summarise(Num_CE_Bins = n_distinct(CE_bin), .groups = "drop") %>%
  filter(Num_CE_Bins == 3)  # Only keep combos with all three bins

# Step 3: Filter original data
filtered_data <- combined_data %>%
  inner_join(valid_combinations, by = c("Metabolite_Name", "Adduct2"))

# Calculate average Weighted_Intensity_Diff per CE_bin for each Metabolite_Name and Adduct2
intensity_summary <- filtered_data %>%
  group_by(Metabolite_Name, Adduct2, CE_bin) %>%
  summarise(Avg_Weighted_Intensity_Diff = mean(Weighted_Intensity_Diff, na.rm = TRUE), .groups = "drop")

# Reshape the data to wide format for easier comparison
intensity_summary_wide <- intensity_summary %>%
  tidyr::pivot_wider(
    names_from = CE_bin,
    values_from = Avg_Weighted_Intensity_Diff,
    names_prefix = "Avg_IntensityDiff_"
  )

# View result
print(intensity_summary_wide)

library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)  # for unit()

# Step 1: Prepare signed mean difference data
ce_bin_numeric <- c("Low" = 1, "Medium" = 2, "High" = 3)

signed_mean_line_data <- intensity_summary_wide %>%
  pivot_longer(
    cols = starts_with("Avg_IntensityDiff_"),
    names_to = "CE_bin",
    values_to = "Mean_Diff"
  ) %>%
  mutate(
    CE_bin = gsub("Avg_IntensityDiff_", "", CE_bin),
    CE_bin = factor(CE_bin, levels = c("Low", "Medium", "High")),
    CE_numeric = ce_bin_numeric[as.character(CE_bin)]
  ) %>%
  group_by(Adduct2, CE_bin, CE_numeric) %>%
  summarise(
    Mean_Diff = round(mean(Mean_Diff, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    Adduct2 = factor(Adduct2, levels = c(
      "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+",
      "[2M+H]+", "[M+Cl]-", "[M-H-H2O]-", "[2M-H]-"
    )),
    Point_Color = case_when(
      Mean_Diff > 0.02  ~ "#E25D75",
      Mean_Diff < -0.02 ~ "#2667FF",
      TRUE              ~ "gray70"
    )
  )

# Step 2: Plot with more spacing between panels
p <- ggplot(signed_mean_line_data, aes(x = CE_numeric, y = Mean_Diff)) +
  geom_point(aes(color = Point_Color), size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2, aes(group = Adduct2)) +
  scale_color_identity() +
  scale_x_continuous(breaks = 1:3, labels = c("Low", "Medium", "High")) +
  facet_wrap(~ Adduct2, nrow = 2, ncol = 4) +
  labs(
    title = "Mean Weighted Intensity Diff Across CE Bins",
    x = "Collision Energy Bin",
    y = "Mean Weighted Intensity Diff"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(face = "bold"),
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.6, "cm"),
    panel.spacing.y = unit(0.8, "cm")
  )

# Step 3: Save with transparent background
ggsave(
  filename = "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/Weighted_mz/Mean_Diff_Plot.tiff",
  plot = p,
  width = 5.78,
  height = 3,
  dpi = 200,
  units = "in",
  bg = "transparent",
  device = "tiff"
)
