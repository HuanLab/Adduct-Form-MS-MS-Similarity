library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)

# === Setup ===
output_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start"
files <- list.files(output_dir, pattern = "_noise5\\.csv$", full.names = TRUE)

# Read all CSV files
combined_data <- bind_rows(lapply(files, function(f) {
  cat("Reading:", basename(f), "\n")
  read.csv(f, stringsAsFactors = FALSE)
}))

# === Weighted Intensity Calc (low m/z only, >1% intensity) ===
compute_weighted_intensity <- function(mz_str, int_str) {
  mz <- as.numeric(strsplit(gsub("\\[|\\]", "", mz_str), ",")[[1]])
  int <- as.numeric(strsplit(gsub("\\[|\\]", "", int_str), ",")[[1]])
  if (length(mz) != length(int) || length(mz) == 0) return(NA)
  threshold <- 0.01 * max(int, na.rm = TRUE)
  keep <- which(int >= threshold & mz < 100)
  if (length(keep) == 0) return(NA)
  sum(mz[keep] * int[keep]) / sum(int[keep])
}

combined_data <- combined_data %>%
  mutate(
    Weighted_Intensity1 = mapply(compute_weighted_intensity, whole_mz1, whole_int1),
    Weighted_Intensity2 = mapply(compute_weighted_intensity, whole_mz2, whole_int2),
    Weighted_Intensity_Diff = (Weighted_Intensity2 - Weighted_Intensity1) / Weighted_Intensity1,
    Average_CE = round((CE1 + CE2) / 2),
    Adduct2 = factor(Adduct2, levels = c("[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+", "[2M+H]+", "[M+Cl]-", "[M-H-H2O]-", "[2M-H]-"))
  ) %>%
  filter(complete.cases(.)) %>%
  mutate(
    CE_bin = case_when(
      Average_CE < 20 ~ "Low",
      Average_CE < 40 ~ "Medium",
      Average_CE <= 60 ~ "High",
      TRUE ~ NA_character_
    )
  )

# === Keep Only Metabolite + Adduct with All 3 Bins ===
valid <- combined_data %>%
  filter(!is.na(CE_bin)) %>%
  count(Metabolite_Name, Adduct2, CE_bin) %>%
  count(Metabolite_Name, Adduct2) %>%
  filter(n == 3)

filtered_data <- semi_join(combined_data, valid, by = c("Metabolite_Name", "Adduct2"))

# === Summary (Wide Format) ===
intensity_summary_wide <- filtered_data %>%
  group_by(Metabolite_Name, Adduct2, CE_bin) %>%
  summarise(Avg_Weighted_Intensity_Diff = mean(Weighted_Intensity_Diff, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = CE_bin, values_from = Avg_Weighted_Intensity_Diff, names_prefix = "Avg_IntensityDiff_")

# === 1. Positive % Heatmap ===
positive_percent_plot_data <- intensity_summary_wide %>%
  pivot_longer(starts_with("Avg_IntensityDiff_"), names_to = "CE_bin", values_to = "IntensityDiff") %>%
  mutate(
    CE_bin = factor(gsub("Avg_IntensityDiff_", "", CE_bin), levels = c("Low", "Medium", "High")),
    Is_Positive = IntensityDiff > 0
  ) %>%
  group_by(Adduct2, CE_bin) %>%
  summarise(Positive_Percent = round(mean(Is_Positive, na.rm = TRUE) * 100), .groups = "drop") %>%
  mutate(Adduct2 = factor(Adduct2, levels = rev(levels(combined_data$Adduct2))))

ggplot(positive_percent_plot_data, aes(x = CE_bin, y = Adduct2)) +
  geom_tile(aes(fill = Positive_Percent), color = "white") +
  geom_text(aes(label = paste0(Positive_Percent, "%")), size = 4.5, fontface = "bold") +
  scale_fill_gradientn(colors = c("#2667FF", "white", "#E25D75"),
                       values = scales::rescale(c(0, 50, 100)), limits = c(0, 100),
                       name = "% Positive") +
  labs(title = "Low m/z Fragment Increase Across CE Bins", x = "Collision Energy Bin", y = "Adduct") +
  theme_minimal(base_size = 13) +
  theme(axis.text = element_text(face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold"))

# === 2. Absolute Mean Line Plot ===
ce_bin_numeric <- c("Low" = 1, "Medium" = 2, "High" = 3)
abs_mean_line_data <- intensity_summary_wide %>%
  pivot_longer(starts_with("Avg_IntensityDiff_"), names_to = "CE_bin", values_to = "IntensityDiff") %>%
  mutate(
    CE_bin = factor(gsub("Avg_IntensityDiff_", "", CE_bin), levels = names(ce_bin_numeric)),
    CE_numeric = ce_bin_numeric[as.character(CE_bin)],
    Abs_IntensityDiff = abs(IntensityDiff)
  ) %>%
  group_by(Adduct2, CE_bin, CE_numeric) %>%
  summarise(Abs_Mean_Intensity = round(mean(Abs_IntensityDiff, na.rm = TRUE), 4), .groups = "drop")

ggplot(abs_mean_line_data, aes(CE_numeric, Abs_Mean_Intensity)) +
  geom_point(color = "black", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2) +
  scale_x_continuous(breaks = 1:3, labels = names(ce_bin_numeric)) +
  facet_wrap(~ Adduct2, nrow = 2, ncol = 4) +
  labs(title = "|Mean Weighted Intensity Diff| Across CE", x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "bold"), panel.grid = element_blank(), axis.text = element_text(face = "bold"))

# === 3. Signed Mean Plot with Color ===
signed_mean_line_data <- intensity_summary_wide %>%
  pivot_longer(starts_with("Avg_IntensityDiff_"), names_to = "CE_bin", values_to = "Mean_Diff") %>%
  mutate(
    CE_bin = factor(gsub("Avg_IntensityDiff_", "", CE_bin), levels = names(ce_bin_numeric)),
    CE_numeric = ce_bin_numeric[as.character(CE_bin)],
    Point_Color = case_when(
      Mean_Diff > 0.02 ~ "#E25D75",
      Mean_Diff < -0.02 ~ "#2667FF",
      TRUE ~ "gray70"
    )
  ) %>%
  group_by(Adduct2, CE_bin, CE_numeric) %>%
  summarise(Mean_Diff = round(mean(Mean_Diff, na.rm = TRUE), 2), Point_Color = first(Point_Color), .groups = "drop")

p <- ggplot(signed_mean_line_data, aes(x = CE_numeric, y = Mean_Diff)) +
  geom_point(aes(color = Point_Color), size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2) +
  scale_color_identity() +
  scale_x_continuous(breaks = 1:3, labels = names(ce_bin_numeric)) +
  facet_wrap(~ Adduct2, nrow = 2, ncol = 4) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.6, "cm")
  )

ggsave(
  filename = file.path(output_dir, "Weighted_mz/Mean_Diff_Plot.tiff"),
  plot = p, width = 5.78, height = 3, dpi = 200, units = "in", bg = "transparent", device = "tiff"
)
