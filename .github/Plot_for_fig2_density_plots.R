library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(parallel)
library(doParallel)

# Common directory paths
main_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start"
nist_dir <- file.path(main_dir, "NIST20")
output_directory <- file.path(main_dir, "Shifted_or_Shifted_Match")
if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)

# Unified adduct list
adducts <- c("[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+", "[2M+H]+", "[M+Cl]-", "[M-H-H2O]-", "[2M-H]-")
adduct_dirs <- setNames(rep(main_dir, length(adducts)), adducts)

# === Part 1: Density plots ===
plot_density_with_normalized_histogram <- function(data, title, binwidth = 5) {
  data <- filter(data, INSTRUMENT_TYPE_1 == INSTRUMENT_TYPE_2)
  make_plot <- function(col, name) {
    values <- pmin(pmax(na.omit(col), 0), 100)
    breaks <- seq(0, 100, by = binwidth)
    hist_vals <- hist(values, breaks = breaks, plot = FALSE)
    hist_df <- data.frame(mid = hist_vals$mids, density = hist_vals$density / sum(hist_vals$density) * 100)
    dens <- density(values)
    scale_factor <- max(hist_df$density) / max(dens$y)
    
    ggplot() +
      geom_bar(data = hist_df, aes(x = mid, y = density), stat = "identity",
               fill = "black", color = "black", alpha = 0.1, width = binwidth) +
      geom_line(aes(x = dens$x, y = dens$y * scale_factor), color = "black", size = 1) +
      labs(
        title = sprintf("%s - %s\nMean = %.1f, Median = %.1f\n(n = %d)",
                        title, name, mean(values), median(values), length(values)),
        x = "Cosine Similarity", y = "Normalized Density (%)"
      ) +
      xlim(0, 100) +
      theme_minimal(base_size = 10) +
      theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, size = 10))
  }
  
  list(
    dp = make_plot(data$dp, "dp"),
    remddp = make_plot(data$remddp, "remddp"),
    FP = make_plot(data$FP, "FP"),
    ts_remddp = make_plot(data$ts_remddp, "ts_remddp")
  )
}

save_plot <- function(plot, adduct, type) {
  ggsave(file.path(nist_dir, sprintf("%s_Identical_Instrument_Density_Plot_%s.tiff", adduct, type)),
         plot, device = "tiff", width = 3, height = 2, dpi = 100)
}

for (adduct in unique(adducts)) {
  file <- file.path(nist_dir, paste0(adduct, "_pairs_df_noise5.csv"))
  if (!file.exists(file)) next
  df <- read.csv(file)
  plots <- plot_density_with_normalized_histogram(df, paste("Identical Instrument Pairs -", adduct))
  lapply(names(plots), function(type) save_plot(plots[[type]], adduct, type))
}

# === Part 2: Match classification and stacked bar ===
NIST20_classified <- read_excel("X:/Users/Botao_Liu/Worklog/Database/NIST20_Pos_classfied.xlsx", col_types = "text")

process_adduct <- function(adduct, directory, class_map) {
  file <- file.path(directory, paste0(gsub("[\\[\\]+-]", "", adduct), "_pairs_df_noise5.csv"))
  df <- read.csv(file, colClasses = "character")
  
  cols <- c("total_aligned", "Direct_match", "Shifted_match")
  df[cols] <- lapply(df[cols], function(x) as.numeric(replace(x, is.na(x), 0)))
  df$superclass <- class_map[df$Metabolite_Name]
  df <- filter(df, !is.na(superclass))
  
  df <- df %>%
    mutate(match_category = case_when(
      total_aligned == 0 ~ "No_match",
      Direct_match > 0 & Shifted_match == 0 ~ "Only_Direct_match",
      Direct_match == 0 & Shifted_match > 0 ~ "Only_Shifted_match",
      Direct_match > 0 & Shifted_match > 0 ~ "Both_Direct_and_Shifted_match",
      TRUE ~ "No_match"
    )) %>%
    count(superclass, match_category, name = "count") %>%
    group_by(superclass) %>%
    mutate(superclass_total = sum(count),
           adduct = adduct,
           percentage = count / superclass_total * 100) %>%
    ungroup()
  
  df
}

class_map <- setNames(NIST20_classified$superclass, NIST20_classified$WHOLE_COMPOUND_NAME)
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
match_results <- foreach(adduct = adducts, .combine = rbind, .packages = c("dplyr")) %dopar% {
  process_adduct(adduct, adduct_dirs[[adduct]], class_map)
}
stopCluster(cl)

write.csv(match_results, file.path(main_dir, "Match_Percentage_by_Superclass.csv"), row.names = FALSE)

plot_data <- match_results %>%
  mutate(adduct = factor(adduct, levels = adducts),
         match_category = factor(match_category, levels = c("Only_Direct_match", "Both_Direct_and_Shifted_match", "Only_Shifted_match", "No_match")),
         superclass_number = as.numeric(factor(superclass)))

ggplot(plot_data, aes(x = factor(superclass_number), y = percentage, fill = match_category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c(
    "Only_Direct_match" = "#2667FF",
    "Both_Direct_and_Shifted_match" = "purple",
    "Only_Shifted_match" = "#E25D75",
    "No_match" = "gray"
  )) +
  facet_wrap(~adduct, scales = "free_y", nrow = 2) +
  ylim(0, 100.1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1),
        strip.text = element_text(size = 10),
        legend.position = "none") +
  ggsave(file.path(main_dir, "Match_Percentage_by_Adducts_Stacked.tiff"), width = 8, height = 4, dpi = 200)

# === Part 3: Line plot for direct vs shifted match ===
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
  geom_line(aes(ID, Shifted_match, color = 'Shifted Match'), linewidth = 1) +
  geom_line(aes(ID, Direct_match, color = 'Direct Match'), linewidth = 1) +
  facet_wrap(~Adduct, ncol = 4) +
  labs(title = "Distribution of Matches by Adduct", x = "Normalized ID", y = "Match Score", color = "Match Type") +
  scale_color_manual(values = c("Shifted Match" = "#E25D75", "Direct Match" = "#2667FF")) +
  theme_minimal() +
  theme(strip.text = element_text(size = 6), panel.grid = element_blank(), legend.position = "none") +
  ggsave(file.path(main_dir, "Adduct_Match_Distribution.tiff"), width = 5.95, height = 3.3, dpi = 300)
