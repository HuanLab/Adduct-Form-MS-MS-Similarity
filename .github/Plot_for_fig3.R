library(dplyr)
library(ggplot2)
library(tidyr)

# Define common directory and adducts
base_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"
adducts <- c("[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+", "[2M+H]+", "[M+Cl]-", "[M-H-H2O]-", "[2M-H]-")
adduct_dirs <- setNames(rep(base_dir, length(adducts)), adducts)

# === Function: Histogram + Density plot ===
plot_fp_density <- function(data, adduct, binwidth = 5) {
  df <- data %>% filter(INSTRUMENT_TYPE_1 == INSTRUMENT_TYPE_2, !is.na(FP)) %>%
    mutate(FP = pmin(pmax(FP, 0), 100))
  
  hist_info <- hist(df$FP, breaks = seq(0, 100, binwidth), plot = FALSE)
  hist_df <- data.frame(mid = hist_info$mids, density = hist_info$density / sum(hist_info$density) * 100)
  
  dens <- density(df$FP)
  scale_factor <- max(hist_df$density) / max(dens$y)
  title_text <- sprintf("Identical Instrument Pairs - %s\nMean = %.1f, Median = %.1f\n(n = %d)", 
                        adduct, mean(df$FP), median(df$FP), nrow(df))
  
  ggplot() +
    geom_bar(data = hist_df, aes(mid, density), stat = "identity", fill = "black", alpha = 0.1, width = binwidth) +
    geom_line(aes(dens$x, dens$y * scale_factor), color = "black", size = 1) +
    xlim(0, 100) + ylim(0, max(hist_df$density) * 1.0) +
    labs(title = title_text, x = "Cosine Similarity", y = "Normalized Density (%)") +
    theme_minimal() +
    theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, size = 10))
}

# === Function: Match category classification ===
classify_matches <- function(df) {
  df <- df %>% 
    mutate(across(c(total_aligned, Direct_match, Shifted_match), ~as.numeric(replace_na(.x, 0)))) %>%
    mutate(match_category = case_when(
      total_aligned == 0 ~ "No_match",
      Direct_match > 0 & Shifted_match == 0 ~ "Only_Direct_match",
      Direct_match == 0 & Shifted_match > 0 ~ "Only_Shifted_match",
      Direct_match > 0 & Shifted_match > 0 ~ "Both_Direct_and_Shifted_match",
      TRUE ~ "No_match"
    ))
  return(df)
}

# === 1. Plot FP Histogram & Density ===
for (adduct in unique(adducts)) {
  file <- file.path(adduct_dirs[[adduct]], paste0(adduct, "_pairs_df_noise5.csv"))
  if (!file.exists(file)) next
  df <- read.csv(file)
  plot <- plot_fp_density(df, adduct)
  ggsave(file.path(base_dir, paste0(adduct, "_Identical_Instrument_Density_Plot_MatchedFragmentsRatio.tiff")),
         plot, width = 3, height = 2, dpi = 100)
}

# === 2. Bar plot of match category proportions ===
match_summary <- bind_rows(lapply(adducts, function(adduct) {
  file <- file.path(adduct_dirs[[adduct]], paste0(gsub("[\\[\\]+-]", "", adduct), "_pairs_df_noise5.csv"))
  if (!file.exists(file)) return(NULL)
  df <- read.csv(file)
  df <- classify_matches(df)
  summary <- df %>%
    count(match_category) %>%
    mutate(adduct = adduct, percentage = n / sum(n)) %>%
    select(adduct, match_category, percentage)
}))

match_summary <- match_summary %>%
  mutate(match_category = factor(match_category, levels = c("Only_Direct_match", "Both_Direct_and_Shifted_match", "Only_Shifted_match", "No_match")),
         adduct = factor(adduct, levels = adducts)) %>%
  complete(adduct, match_category, fill = list(percentage = 0))

ggplot(match_summary, aes(adduct, percentage, fill = match_category)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Only_Direct_match" = "#2667FF",
                               "Both_Direct_and_Shifted_match" = "purple",
                               "Only_Shifted_match" = "#E25D75",
                               "No_match" = "gray")) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), legend.position = "none") +
  ggsave(file.path(base_dir, "Adduct_Match_Categories_Stacked.tiff"), width = 4, height = 3, dpi = 200)

write.csv(match_summary, file.path(base_dir, "Adduct_Match_Categories_Data.csv"), row.names = FALSE)

# === 3. Plot match intensity vs. normalized ID ===
all_data <- bind_rows(lapply(adducts, function(adduct) {
  file <- file.path(adduct_dirs[[adduct]], paste0(gsub("[\\[\\]+-]", "", adduct), "_pairs_df_noise5.csv"))
  if (!file.exists(file)) return(NULL)
  df <- read.csv(file)
  if (!all(c("Shifted_match", "Direct_match", "total_aligned") %in% names(df))) return(NULL)
  df %>%
    filter(total_aligned != 0) %>%
    arrange(desc(Shifted_match)) %>%
    mutate(ID = row_number() / n(), Adduct = adduct)
}))

all_data$Adduct <- factor(all_data$Adduct, levels = adducts)

ggplot(all_data) +
  geom_line(aes(ID, Shifted_match, color = "Shifted Match"), linewidth = 1) +
  geom_line(aes(ID, Direct_match, color = "Direct Match"), linewidth = 1) +
  facet_wrap(~Adduct, ncol = 4) +
  labs(title = "Distribution of Matches by Adduct", x = "Normalized ID", y = "Match Score") +
  scale_color_manual(values = c("Shifted Match" = "#E25D75", "Direct Match" = "#2667FF")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 5), legend.position = "none") +
  ggsave(file.path(base_dir, "Adduct_Match_Distribution.tiff"), width = 5.95, height = 3.3, dpi = 300)
