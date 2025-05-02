library(dplyr)
library(tidyr)
library(openxlsx)

# Load data
data <- read.xlsx("C:/Program Files/Botao/Worklog/20240923_RF_POS/WHOLE_dataframe_expanded.xlsx")

# Filter for phosphate-containing compounds
phosphate_data <- data %>%
  filter(Phosphate_SMILE == 1)

# Initialize a data frame to store results for each CE threshold
classification_results <- data.frame(
  CE_Threshold = numeric(),
  Classification = character(),
  Good_Quality_Percentage = numeric(),
  Count = integer()
)

# Loop through different CE thresholds (0, 10, 20, ..., 200)
for (ce_threshold in seq(0, 200, by = 10)) {
  
  # Filter data for rows where WHOLE_COLLISION_ENERGY is greater than the current threshold
  filtered_data <- phosphate_data[phosphate_data$WHOLE_COLLISION_ENERGY > ce_threshold, ]
  
  # Classify adducts into three categories
  filtered_data <- filtered_data %>%
    mutate(Classification = case_when(
      grepl("Na|K", WHOLE_ADDUCT) ~ "Alkali Adducts",
      TRUE ~ "Excluded Alkali Adducts"
    ))
  
  # Add a third classification for "All Adducts"
  filtered_data <- filtered_data %>%
    bind_rows(filtered_data %>% mutate(Classification = "All Adducts"))
  
  # Calculate GQP and counts for each classification
  classification_summary <- filtered_data %>%
    group_by(Classification) %>%
    summarise(
      Good_Quality_Percentage = mean(rowSums(across(c(frag98_presence, frag80_presence, NL80_presence, NL98_presence), ~ . == 1)) > 0) * 100,
      Count = n()
    ) %>%
    mutate(CE_Threshold = ce_threshold)
  
  # Append results to the main data frame
  classification_results <- bind_rows(classification_results, classification_summary)
}

# Reshape the data to have CE threshold in column 1 and GQP for each classification in columns 2-4
reshaped_results <- classification_results %>%
  select(CE_Threshold, Classification, Good_Quality_Percentage) %>%
  pivot_wider(names_from = Classification, values_from = Good_Quality_Percentage)

# View the reshaped results
print(reshaped_results)

write.xlsx(reshaped_results, "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/reshaped_gqp_counts_results.xlsx", rowNames = FALSE)






#### Model training ####










results_lists_1 <- results_lists
results_lists_2 <- results_lists

# Add source column to differentiate data sources
combined_roc_df_1 <- do.call(rbind, lapply(1:length(results_lists_1$prob_predictions), function(i) {
  df <- plot_roc_curve(
    results_lists_1$prob_predictions[[i]], 
    results_lists_1$actual_values[[i]], 
    hyperparameter_tuning$CE_Threshold[i],
    combined = TRUE
  )
  df$Source <- "Data Source 1"
  return(df)
}))

combined_roc_df_2 <- do.call(rbind, lapply(1:length(results_lists_2$prob_predictions), function(i) {
  df <- plot_roc_curve(
    results_lists_2$prob_predictions[[i]], 
    results_lists_2$actual_values[[i]], 
    hyperparameter_tuning$CE_Threshold[i],
    combined = TRUE
  )
  df$Source <- "Data Source 2"
  return(df)
}))


library(openxlsx)

# Save combined_roc_df_1 as an Excel file
output_file_1 <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/combined_roc_df_1.xlsx"
write.xlsx(combined_roc_df_1, output_file_1, rowNames = FALSE)

# Save combined_roc_df_2 as an Excel file
output_file_2 <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/combined_roc_df_2.xlsx"
write.xlsx(combined_roc_df_2, output_file_2, rowNames = FALSE)




# Add the Source column directly when plotting
p_combined <- ggplot() +
  # Plot for Data Source 1
  geom_line(data = combined_roc_df_1, aes(x = FPR, y = TPR, color = "Data Source 1", linetype = factor(CE_Threshold)), size = 2) +
  # Plot for Data Source 2
  geom_line(data = combined_roc_df_2, aes(x = FPR, y = TPR, color = "Data Source 2", linetype = factor(CE_Threshold)), size = 2) +
  # Add diagonal reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", size = 1) +
  # Customize axis and labels
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), name = "1 - Specificity") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), name = "Sensitivity") +
  scale_color_manual(values = c("Data Source 1" = "#3F8EFC", "Data Source 2" = "#ff6060")) +  # Custom colors
  scale_linetype_discrete(name = "CE Threshold") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 2),  # Thicker frame
    aspect.ratio = 1,
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove title
  )

# Save the plot
ggsave(file.path(output_dir, "combined_ROC_curves.tiff"), 
       p_combined, 
       device = "tiff", 
       width = 4, 
       height = 3.2, 
       dpi = 100)

# Print the plot
print(p_combined)
