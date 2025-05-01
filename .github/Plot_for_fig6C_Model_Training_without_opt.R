# Load necessary libraries
library(openxlsx)
library(dplyr)
library(randomForest)
library(tidyr)
library(ggplot2)

# Set working directory
setwd("C:/Program Files/Botao/Worklog/20240923_RF_POS")
# Read input data
WHOLE_dataframe_expanded <- read.xlsx("WHOLE_dataframe_expanded.xlsx")


library(dplyr)

WHOLE_dataframe_expanded <- WHOLE_dataframe_expanded %>%
  mutate(
    NL80_mz_diff = ifelse(NL80_presence == 0, 0, NL80_mz_diff),
    NL98_mz_diff = ifelse(NL98_presence == 0, 0, NL98_mz_diff)
  )


# Set working directory
setwd("X:/Users/Botao_Liu/Worklog/Adduct/20250305_Fig6")
# Read input data
hyperparameter_tuning <- read.xlsx("hyperparameter_tuning_without_opt_Pos.xlsx")

# Set working directory
setwd("X:/Users/Botao_Liu/Worklog/Adduct/20250305_Fig6")
# Set the output directory for figures
output_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20250305_Fig6"


if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ----------------- 🔹 Parameters -----------------
iterations <- 100
ntree <- 200
mtry <- 5

# Initialize lists to store all results
results_lists <- list(
  TPR = list(),
  FPR = list(),
  model_filenames = list(),
  feature_importance = list(),
  conf_matrix = list(),
  prob_predictions = list(),
  actual_values = list()
)

# Selected features for model training
selected_features <- c(
  "WHOLE_PRECURSOR_MZ",
  "frag98_presence", "frag98_mz", "frag98_mz_diff", "frag98_int",
  "frag80_presence", "frag80_mz", "frag80_mz_diff", "frag80_int",
  "NL80_presence", "NL80_Count", "NL80_mz", "NL80_mz_diff", "NL80_parent_int", "NL80_product_int",
  "NL98_presence", "NL98_Count", "NL98_mz", "NL98_mz_diff", "NL98_parent_int", "NL98_product_int",
  "Phosphate_SMILE"
)

# Function to calculate ROC points and AUC
calculate_roc <- function(predictions, labels) {
  ord <- order(predictions, decreasing = TRUE)
  predictions <- predictions[ord]
  labels <- labels[ord]
  tpr <- cumsum(labels) / sum(labels)
  fpr <- cumsum(!labels) / sum(!labels)
  tpr <- c(0, tpr, 1)
  fpr <- c(0, fpr, 1)
  auc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)
  list(tpr = tpr, fpr = fpr, auc = auc)
}

# Main processing loop
set.seed(1234)  # For reproducibility

for (i in 1:nrow(hyperparameter_tuning)) {
  # Extract parameters
  CE_threshold <- hyperparameter_tuning$CE_Threshold[i]
  desired_adducts <- unlist(strsplit(as.character(hyperparameter_tuning$Selected_Adducts[i]), ", "))
  
  # Filter dataset
  filtered_dataframe_expanded <- WHOLE_dataframe_expanded %>%
    filter(WHOLE_ADDUCT %in% desired_adducts, WHOLE_COLLISION_ENERGY >= CE_threshold) %>%
    drop_na()
  
  # Count phosphate cases
  phosphate_cases <- sum(filtered_dataframe_expanded$Phosphate_SMILE == 1)
  
  # Initialize iteration storage
  iteration_results <- list(
    feature_importance = list(),
    conf_matrix = list(),
    predictions = numeric(),
    actuals = numeric()
  )
  
  for (j in 1:iterations) {
    # Balance dataset
    phosphate_data <- filtered_dataframe_expanded %>% filter(Phosphate_SMILE == 1)
    non_phosphate_data <- filtered_dataframe_expanded %>% 
      filter(Phosphate_SMILE == 0) %>%
      sample_n(phosphate_cases)
    
    balanced_data <- bind_rows(phosphate_data, non_phosphate_data)
    balanced_data$Phosphate_SMILE <- as.factor(balanced_data$Phosphate_SMILE)
    
    # Split dataset
    training_indices <- balanced_data %>%
      group_by(Phosphate_SMILE) %>%
      sample_frac(0.8) %>%
      ungroup() %>%
      pull(index)
    
    train_data <- balanced_data[balanced_data$index %in% training_indices, selected_features]
    test_data <- balanced_data[!(balanced_data$index %in% training_indices), selected_features]
    
    # Train model
    model <- randomForest(Phosphate_SMILE ~ ., data = train_data,
                          importance = TRUE, ntree = ntree, mtry = mtry)
    
    # Store feature importance
    feature_importance <- as.data.frame(importance(model))
    feature_importance$Feature <- rownames(feature_importance)
    feature_importance$CE_Threshold <- CE_threshold
    iteration_results$feature_importance[[j]] <- feature_importance
    
    # Make predictions
    predictions <- predict(model, test_data, type = "class")
    prob_predictions <- predict(model, test_data, type = "prob")[,2]
    
    # Store predictions and actual values
    iteration_results$predictions <- c(iteration_results$predictions, prob_predictions)
    iteration_results$actuals <- c(iteration_results$actuals, as.numeric(as.character(test_data$Phosphate_SMILE)))
    
    # Compute confusion matrix
    conf_matrix <- table(Actual = test_data$Phosphate_SMILE, Predicted = predictions)
    conf_matrix_df <- as.data.frame.matrix(conf_matrix)
    
    # Ensure consistent columns
    all_classes <- c("0", "1")
    for (class in all_classes) {
      if (!(class %in% colnames(conf_matrix_df))) {
        conf_matrix_df[[class]] <- 0
      }
    }
    conf_matrix_df <- conf_matrix_df[, all_classes]
    iteration_results$conf_matrix[[j]] <- conf_matrix_df
  }
  
  # Process feature importance
  feature_importance_combined <- do.call(rbind, iteration_results$feature_importance)
  feature_importance_avg <- feature_importance_combined %>%
    group_by(Feature, CE_Threshold) %>%
    summarise(
      Mean_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy),
      SD_MeanDecreaseAccuracy = sd(MeanDecreaseAccuracy),
      Mean_MeanDecreaseGini = mean(MeanDecreaseGini),
      SD_MeanDecreaseGini = sd(MeanDecreaseGini),
      .groups = "drop"
    )
  
  # Process confusion matrix
  conf_matrix_avg <- Reduce("+", iteration_results$conf_matrix) / iterations
  conf_matrix_avg$CE_Threshold <- CE_threshold
  
  # Calculate metrics
  TP <- conf_matrix_avg["1", "1"]
  FN <- conf_matrix_avg["1", "0"]
  TN <- conf_matrix_avg["0", "0"]
  FP <- conf_matrix_avg["0", "1"]
  
  avg_TPR <- TP / (TP + FN)
  avg_FPR <- FP / (FP + TN)
  
  # Save the model with a descriptive filename
  filename <- paste("random_forest_model_CE", CE_threshold, 
                    "_ntree", ntree, "_mtry", mtry, "_iterations", iterations, ".rds", sep = "")
  saveRDS(model, file = filename)
  
  # Store all results
  results_lists$TPR[[i]] <- avg_TPR
  results_lists$FPR[[i]] <- avg_FPR
  results_lists$feature_importance[[i]] <- feature_importance_avg
  results_lists$conf_matrix[[i]] <- conf_matrix_avg
  results_lists$model_filenames[[i]] <- paste0("RF_Model_CE_", CE_threshold, ".rds")
  results_lists$prob_predictions[[i]] <- iteration_results$predictions
  results_lists$actual_values[[i]] <- iteration_results$actuals
}

# Compile results into dataframes
results <- data.frame(
  CE_Threshold = hyperparameter_tuning$CE_Threshold,
  Average_TPR = unlist(results_lists$TPR),
  Average_FPR = unlist(results_lists$FPR),
  Model_Filename = unlist(results_lists$model_filenames)
)

average_feature_importance_df <- do.call(rbind, results_lists$feature_importance)

average_conf_matrix_df <- do.call(rbind, lapply(1:length(results_lists$conf_matrix), function(i) {
  data.frame(
    CE_Threshold = hyperparameter_tuning$CE_Threshold[i],
    TP = results_lists$conf_matrix[[i]][2,2],
    FN = results_lists$conf_matrix[[i]][2,1],
    TN = results_lists$conf_matrix[[i]][1,1],
    FP = results_lists$conf_matrix[[i]][1,2],
    Average_TPR = results_lists$TPR[[i]],
    Average_FPR = results_lists$FPR[[i]]
  )
}))

# Visualization functions
plot_feature_importance <- function(feature_importance_df, CE_value) {
  feature_importance_subset <- feature_importance_df %>%
    filter(CE_Threshold == CE_value) %>%
    arrange(desc(Mean_MeanDecreaseAccuracy))
  
  p <- ggplot(feature_importance_subset, 
              aes(x = reorder(Feature, Mean_MeanDecreaseAccuracy), 
                  y = Mean_MeanDecreaseAccuracy)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_errorbar(aes(ymin = Mean_MeanDecreaseAccuracy - SD_MeanDecreaseAccuracy, 
                      ymax = Mean_MeanDecreaseAccuracy + SD_MeanDecreaseAccuracy),
                  width = 0.5, color = "black") +
    coord_flip() +
    labs(title = paste("Feature Importance (CE Threshold =", CE_value, ")"),
         x = "Features", y = "MeanDecreaseAccuracy ± SD") +
    theme_minimal() +
    theme(panel.grid = element_blank())
  
  # Save the plot
  filename <- file.path(output_dir, sprintf("feature_importance_CE%.1f.tiff", CE_value))
  ggsave(filename, p, device = "tiff", width = 10, height = 8, dpi = 600)
  
  return(p)
}

plot_roc_curve <- function(predictions, actuals, CE_threshold, combined = FALSE) {
  roc_points <- calculate_roc(predictions, actuals)
  
  if (!combined) {
    p <- ggplot(data.frame(FPR = roc_points$fpr, TPR = roc_points$tpr), 
                aes(x = FPR, y = TPR)) +
      geom_line(color = "red", size = 1) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                         name = "1 - Specificity") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                         name = "Sensitivity") +
      ggtitle(paste("ROC Curve (CE Threshold =", CE_threshold, ")")) +
      annotate("text", x = 0.75, y = 0.25, 
               label = sprintf("AUC = %.3f", roc_points$auc),
               size = 4) +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        aspect.ratio = 1
      )
    
    # Save the plot
    filename <- file.path(output_dir, sprintf("ROC_curve_CE%.1f.tiff", CE_threshold))
    ggsave(filename, p, device = "tiff", width = 8, height = 8, dpi = 600)
    
    return(p)
  } else {
    return(data.frame(
      FPR = roc_points$fpr,
      TPR = roc_points$tpr,
      CE_Threshold = CE_threshold
    ))
  }
}

# Generate and save feature importance plots
for (CE_value in unique(average_feature_importance_df$CE_Threshold)) {
  print(plot_feature_importance(average_feature_importance_df, CE_value))
}

# Generate and save individual ROC plots
for (i in 1:length(results_lists$prob_predictions)) {
  print(plot_roc_curve(
    results_lists$prob_predictions[[i]], 
    results_lists$actual_values[[i]], 
    hyperparameter_tuning$CE_Threshold[i]
  ))
}

# Generate and save combined ROC plot
combined_roc_df <- do.call(rbind, lapply(1:length(results_lists$prob_predictions), function(i) {
  plot_roc_curve(
    results_lists$prob_predictions[[i]], 
    results_lists$actual_values[[i]], 
    hyperparameter_tuning$CE_Threshold[i],
    combined = TRUE
  )
}))

p_combined <- ggplot(combined_roc_df, aes(x = FPR, y = TPR, color = factor(CE_Threshold))) +
  geom_line(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     name = "1 - Specificity") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     name = "Sensitivity") +
  scale_color_discrete(name = "CE Threshold") +
  ggtitle("Combined ROC Curves for Different CE Thresholds") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    aspect.ratio = 1,
    legend.position = "right"
  )

# Save combined ROC plot
ggsave(file.path(output_dir, "combined_ROC_curves.tiff"), 
       p_combined, 
       device = "tiff", 
       width = 10, 
       height = 8, 
       dpi = 600)

print(p_combined)





# Set working directory
setwd("X:/Users/Botao_Liu/Worklog/Adduct/20250305_Fig6")

write.xlsx(average_conf_matrix_df, "average_conf_matrix_df.xlsx")
write.xlsx(average_feature_importance_df, "average_feature_importance_df.xlsx")

# Function to convert results_lists to MGF format
convert_to_mgf <- function(results_list, file_path) {
  # Open the file for writing
  mgf_file <- file(file_path, "w")
  
  for (i in 1:length(results_list$prob_predictions)) {
    # Write the BEGIN IONS line
    cat("BEGIN IONS\n", file = mgf_file)
    
    # Assuming results_list has relevant fields like 'm/z' and 'intensity' that we want to save
    # Adjust the following line to match your data structure
    mz_values <- results_list$mz[[i]]  # Modify according to actual structure
    intensity_values <- results_list$intensity[[i]]  # Modify according to actual structure
    
    for (j in 1:length(mz_values)) {
      # Write each m/z and intensity pair
      cat(paste(mz_values[j], intensity_values[j]), "\n", file = mgf_file)
    }
    
    # Write the END IONS line
    cat("END IONS\n", file = mgf_file)
  }
  
  # Close the MGF file
  close(mgf_file)
}

# Save results_lists to an MGF file
convert_to_mgf(results_lists, "path_to_output_directory/results_lists.mgf")
