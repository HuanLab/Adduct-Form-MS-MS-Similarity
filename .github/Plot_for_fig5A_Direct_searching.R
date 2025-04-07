# Load required libraries
library(dplyr)
library(openxlsx)
library(readxl)
library(writexl)
library(ggplot2)
library(parallel)
library(doParallel)
library(tidyr)

# Set working directory
base_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"
adduct_dir <- file.path(base_dir, "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20")
setwd(adduct_dir)


# ---- SAMPLE PREPARATION ----
# Set seed for reproducible sampling
set.seed(1234)

# Reload the classified data
pairs_df <- read_excel("Na_pairs_df_classify.xlsx")

# Filter for identical instrument types
pairs_df <- pairs_df %>% filter(INSTRUMENT_TYPE_1 == INSTRUMENT_TYPE_2)
pairs_df <- pairs_df %>% filter(CE1 > 10 & CE1 < 60)

# Identify superclasses with >100 occurrences
superclass_counts <- table(pairs_df$superclass)
superclass_filtered <- names(superclass_counts[superclass_counts > 100])


# Sample 100 items from each filtered superclass
sampled_data <- pairs_df %>%
  filter(superclass %in% superclass_filtered) %>%
  group_by(superclass) %>%
  filter(!is.na(Metabolite_Name)) %>%
  slice_sample(n = 100, replace = FALSE) %>%
  ungroup() %>%
  mutate(row_index = row_number())

pairs_df <- process_spectral_matching(pairs_df)

# Save the sampled data
write.csv(sampled_data, "sampled_metabolites_by_superclass.csv", row.names = FALSE)

# ---- 4. MATCH H AND NA ADDUCTS ----
# Load H adduct data
adduct_H <- read.xlsx(file.path(base_dir, "adduct_H.xlsx"))

# Create backup of original sampled data
Na_vs_H <- sampled_data

# Find matching H adducts for the sampled Na adducts
H_vs_H <- sampled_data

# Process each row to find H adduct matches
# Process each row to find H adduct matches
for (i in 1:nrow(H_vs_H)) {
  metabolite_name <- H_vs_H$Metabolite_Name[i]
  instrument <- H_vs_H$INSTRUMENT_TYPE_1[i]
  ce2 <- H_vs_H$CE2[i]
  ce1 <- H_vs_H$CE1[i]
  
  # Find matching metabolites with the same instrument, but different CE1
  matched_rows <- adduct_H %>%
    filter(WHOLE_COMPOUND_NAME == metabolite_name, 
           INSTRUMENT_TYPE == instrument,
           WHOLE_COLLISION_ENERGY != ce1) %>%
    filter(abs(WHOLE_COLLISION_ENERGY - ce2) < 3)  # Apply new constraint
  
  if (nrow(matched_rows) > 0) {
    # Find closest collision energy match
    closest_match_index <- which.min(abs(matched_rows$WHOLE_COLLISION_ENERGY - ce2))
    closest_match <- matched_rows[closest_match_index, ]
    
    # Update with H+ adduct data
    H_vs_H$CE2[i] <- closest_match$WHOLE_COLLISION_ENERGY
    H_vs_H$Adduct2[i] <- "[M+H]+"
    H_vs_H$Index_2[i] <- closest_match$WHOLE_Index
    H_vs_H$precursor2[i] <- closest_match$WHOLE_PRECURSOR_mz
    H_vs_H$whole_mz2[i] <- closest_match$WHOLE_mz
    H_vs_H$whole_int2[i] <- closest_match$WHOLE_Int
    H_vs_H$whole_mz_cleaned_2[i] <- closest_match$WHOLE_mz_cleaned
    H_vs_H$whole_int_cleaned_2[i] <- closest_match$WHOLE_Int_cleaned
  } else {
    H_vs_H$Adduct2[i] <- NA
    H_vs_H$CE2[i] <- NA
  }
}

# Remove rows without matches
H_vs_H <- H_vs_H %>% filter(!is.na(Adduct2))

# Identify superclasses with >100 occurrences
superclass_counts <- table(H_vs_H$superclass)
superclass_filtered <- names(superclass_counts[superclass_counts > 50])


# Sample 50 rows from each superclass for final dataset
set.seed(1234)
H_vs_H <- H_vs_H%>%
  filter(superclass %in% superclass_filtered) %>%
  group_by(superclass) %>%
  filter(!is.na(Metabolite_Name)) %>%
  slice_sample(n = 50, replace = FALSE) %>%
  ungroup()

# Filter Na_vs_H to match H_vs_H
Na_vs_H <- Na_vs_H %>% filter(row_index %in% H_vs_H$row_index)

# Save intermediate results
write.xlsx(H_vs_H, "processed_sampled_data.xlsx")

# ---- 5. IMPROVED SPECTRAL SIMILARITY CALCULATION ----
# Convert precursor values to numeric
pairs_df <- H_vs_H
pairs_df$precursor1 <- as.numeric(pairs_df$precursor1)
pairs_df$precursor2 <- as.numeric(pairs_df$precursor2)

# Define improved similarity calculation function
calculate_similarities <- function(ref_fragMz, ref_fragInt, exp_fragMz, exp_fragInt, precursor_mz_diff, mz_tol) {
  # Handle empty data
  if (length(ref_fragMz) == 0 || length(exp_fragMz) == 0) {
    return(list(
      dot_product = 0, 
      modified_cosine_similarity = 0,
      aligned_fragments = 0,
      remddp_aligned_fragments = 0,
      matched_pairs = data.frame()
    ))
  }
  
  # Create a data frame to track all potential matches
  matches <- data.frame(
    ref_idx = integer(),
    exp_idx = integer(),
    ref_mz = numeric(),
    exp_mz = numeric(),
    ref_int = numeric(),
    exp_int = numeric(),
    diff = numeric(),
    is_modified = logical(),
    stringsAsFactors = FALSE
  )
  
  # Find all potential matches within tolerance
  for (i in seq_along(ref_fragMz)) {
    for (j in seq_along(exp_fragMz)) {
      # Direct match
      diff_direct <- abs(ref_fragMz[i] - exp_fragMz[j])
      if (diff_direct <= mz_tol) {
        matches <- rbind(matches, data.frame(
          ref_idx = i,
          exp_idx = j,
          ref_mz = ref_fragMz[i],
          exp_mz = exp_fragMz[j],
          ref_int = ref_fragInt[i],
          exp_int = exp_fragInt[j],
          diff = diff_direct,
          is_modified = FALSE,
          stringsAsFactors = FALSE
        ))
      }
      
      # Modified match (accounting for precursor mass shift)
      diff_modified <- abs((ref_fragMz[i] + precursor_mz_diff) - exp_fragMz[j])
      if (diff_modified <= mz_tol) {
        matches <- rbind(matches, data.frame(
          ref_idx = i,
          exp_idx = j,
          ref_mz = ref_fragMz[i],
          exp_mz = exp_fragMz[j],
          ref_int = ref_fragInt[i],
          exp_int = exp_fragInt[j],
          diff = diff_modified,
          is_modified = TRUE,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # If no matches found, return zeros
  if (nrow(matches) == 0) {
    return(list(
      dot_product = 0, 
      modified_cosine_similarity = 0,
      aligned_fragments = 0,
      remddp_aligned_fragments = 0,
      matched_pairs = data.frame()
    ))
  }
  
  # Sort matches by difference (smallest first)
  matches <- matches[order(matches$diff), ]
  
  # Select best matches (lowest mass difference first, then highest intensity product)
  selected_matches <- data.frame(
    ref_idx = integer(),
    exp_idx = integer(),
    ref_mz = numeric(),
    exp_mz = numeric(),
    ref_int = numeric(),
    exp_int = numeric(),
    diff = numeric(),
    is_modified = logical(),
    stringsAsFactors = FALSE
  )
  
  used_ref <- logical(length(ref_fragMz))
  used_exp <- logical(length(exp_fragMz))
  
  # First handle direct matches
  direct_matches <- matches[!matches$is_modified, ]
  if (nrow(direct_matches) > 0) {
    for (i in 1:nrow(direct_matches)) {
      ref_idx <- direct_matches$ref_idx[i]
      exp_idx <- direct_matches$exp_idx[i]
      
      if (!used_ref[ref_idx] && !used_exp[exp_idx]) {
        used_ref[ref_idx] <- TRUE
        used_exp[exp_idx] <- TRUE
        
        selected_matches <- rbind(selected_matches, direct_matches[i, ])
      }
    }
  }
  
  # Then handle modified matches for remaining unmatched fragments
  modified_matches <- matches[matches$is_modified, ]
  if (nrow(modified_matches) > 0) {
    for (i in 1:nrow(modified_matches)) {
      ref_idx <- modified_matches$ref_idx[i]
      exp_idx <- modified_matches$exp_idx[i]
      
      if (!used_ref[ref_idx] && !used_exp[exp_idx]) {
        used_ref[ref_idx] <- TRUE
        used_exp[exp_idx] <- TRUE
        
        selected_matches <- rbind(selected_matches, modified_matches[i, ])
      }
    }
  }
  
  # Create vectors for dot product calculation
  dot_exp_int_original <- numeric(length(ref_fragMz))
  dot_exp_int_modified <- numeric(length(ref_fragMz))
  
  # Fill in matched intensities for direct matches
  direct_selected <- selected_matches[!selected_matches$is_modified, ]
  if (nrow(direct_selected) > 0) {
    for (i in 1:nrow(direct_selected)) {
      ref_idx <- direct_selected$ref_idx[i]
      dot_exp_int_original[ref_idx] <- direct_selected$exp_int[i]
    }
  }
  
  # Fill in matched intensities for modified matches
  modified_selected <- selected_matches[selected_matches$is_modified, ]
  if (nrow(modified_selected) > 0) {
    for (i in 1:nrow(modified_selected)) {
      ref_idx <- modified_selected$ref_idx[i]
      dot_exp_int_modified[ref_idx] <- modified_selected$exp_int[i]
    }
  }
  
  # Combined matches for modified cosine
  combined_exp_int <- pmax(dot_exp_int_original, dot_exp_int_modified)
  
  # Calculate normalization factors
  norm_ref <- sqrt(sum(ref_fragInt^2))
  norm_exp <- sqrt(sum(exp_fragInt^2))
  
  # Calculate dot product
  dot_product <- if (norm_ref == 0 || norm_exp == 0) 0 else 
    100 * sum(ref_fragInt * dot_exp_int_original) / (norm_ref * norm_exp)
  
  # Calculate modified cosine similarity
  modified_cosine <- if (norm_ref == 0 || norm_exp == 0) 0 else 
    100 * sum(ref_fragInt * combined_exp_int) / (norm_ref * norm_exp)
  
  # Count matched fragments
  aligned_fragments <- sum(dot_exp_int_original > 0)
  remddp_aligned_fragments <- sum(dot_exp_int_modified > 0)
  
  # Return all results
  return(list(
    dot_product = dot_product,
    modified_cosine_similarity = modified_cosine,
    aligned_fragments = aligned_fragments,
    remddp_aligned_fragments = remddp_aligned_fragments,
    matched_pairs = selected_matches
  ))
}

# Process row function for parallel computation
process_row <- function(i, pairs_data, noise_threshold, mz_tol) {
  if (i %% 100 == 0) cat(paste("Processing row", i, "\n"))
  
  # Extract precursor m/z values
  precursor1 <- pairs_data$precursor1[i]
  precursor2 <- pairs_data$precursor2[i]
  
  # Parse spectral data - using whole_mz and whole_int directly
  # instead of cleaned versions to avoid filtering out important peaks
  whole_mz1 <- as.numeric(unlist(strsplit(as.character(pairs_data$whole_mz1[i]), ",")))
  whole_int1 <- as.numeric(unlist(strsplit(as.character(pairs_data$whole_int1[i]), ",")))
  whole_mz2 <- as.numeric(unlist(strsplit(as.character(pairs_data$whole_mz2[i]), ",")))
  whole_int2 <- as.numeric(unlist(strsplit(as.character(pairs_data$whole_int2[i]), ",")))
  
  # Apply noise threshold
  valid_indices_1 <- whole_int1 >= noise_threshold
  valid_indices_2 <- whole_int2 >= noise_threshold
  
  whole_mz1 <- whole_mz1[valid_indices_1]
  whole_int1 <- whole_int1[valid_indices_1]
  whole_mz2 <- whole_mz2[valid_indices_2]
  whole_int2 <- whole_int2[valid_indices_2]
  
  # Filter fragments below precursor (apply less stringent threshold)
  larger_precursor_mz <- min(precursor1, precursor2)
  mz_threshold_dot <- larger_precursor_mz + 0.02  # Include precursors
  
  ref_filter_dot <- whole_mz2 <= mz_threshold_dot
  exp_filter_dot <- whole_mz1 <= mz_threshold_dot
  
  ref_fragMz_dot <- whole_mz2[ref_filter_dot]
  ref_fragInt_dot <- whole_int2[ref_filter_dot]
  exp_fragMz_dot <- whole_mz1[exp_filter_dot]
  exp_fragInt_dot <- whole_int1[exp_filter_dot]
  
  # Calculate precursor difference
  precursor_mz_diff <- precursor2 - precursor1
  
  # Calculate similarities
  result <- calculate_similarities(
    ref_fragMz_dot, ref_fragInt_dot,
    exp_fragMz_dot, exp_fragInt_dot,
    precursor_mz_diff, mz_tol
  )
  
  # Prepare full results
  return(list(
    dp = result$dot_product, 
    remddp = result$modified_cosine_similarity,
    dp_aligned_fragments = result$aligned_fragments,
    remddp_aligned_fragments = result$remddp_aligned_fragments,
    matched_pairs = result$matched_pairs
  ))
}

# Main processing function
process_spectral_matching <- function(pairs_df, noise_threshold = 5, mz_tol = 0.02) {
  # Setup parallel processing
  num_cores <- max(1, detectCores() - 1)
  cat(paste("Setting up parallel processing with", num_cores, "cores\n"))
  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Export necessary libraries and functions
  clusterEvalQ(cl, { library(dplyr) })
  clusterExport(cl, c("process_row", "calculate_similarities"))
  
  # Process in parallel
  start_time <- Sys.time()
  cat("Starting parallel processing...\n")
  
  results <- parLapply(cl, 1:nrow(pairs_df), function(i) {
    process_row(i, pairs_df, noise_threshold, mz_tol)
  })
  
  stopCluster(cl)
  
  end_time <- Sys.time()
  cat(paste("Processing completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n"))
  
  # Extract matched pairs for each row (for debugging)
  matched_pairs_list <- lapply(results, function(res) res$matched_pairs)
  
  # Add results to dataframe
  pairs_df$dp <- sapply(results, `[[`, "dp")
  pairs_df$remddp <- sapply(results, `[[`, "remddp")
  pairs_df$dp_aligned_fragments <- sapply(results, `[[`, "dp_aligned_fragments")
  pairs_df$remddp_aligned_fragments <- sapply(results, `[[`, "remddp_aligned_fragments")
  
  # Save detailed match info to a separate file (optional)
  saveRDS(matched_pairs_list, "detailed_fragment_matches.rds")
  
  return(pairs_df)
}

# Function to test a specific row (for debugging)
test_single_row <- function(pairs_df, row_id, noise_threshold = 5, mz_tol = 0.02) {
  result <- process_row(row_id, pairs_df, noise_threshold, mz_tol)
  
  # Print detailed results
  cat("Row ID:", row_id, "\n")
  cat("Metabolite:", pairs_df$Metabolite_Name[row_id], "\n")
  cat("Dot Product:", result$dp, "\n")
  cat("Modified Cosine Similarity:", result$remddp, "\n")
  cat("Aligned Fragments:", result$dp_aligned_fragments, "\n")
  cat("Modified Aligned Fragments:", result$remddp_aligned_fragments, "\n")
  
  # Print match details
  if (nrow(result$matched_pairs) > 0) {
    cat("\nMatched Fragment Pairs:\n")
    print(result$matched_pairs)
  } else {
    cat("\nNo matched fragments found.\n")
  }
  
  return(result)
}

# Run spectral matching
H_vs_H <- process_spectral_matching(pairs_df)

# Save results
write_xlsx(H_vs_H, "H_vs_H.xlsx")
write_xlsx(Na_vs_H, "Na_vs_H.xlsx")

# ---- 6. VISUALIZATION AND ANALYSIS ----
# Reload saved results
H_vs_H <- read.xlsx("H_vs_H.xlsx")
Na_vs_H <- read.xlsx("Na_vs_H.xlsx")

# Remove duplicate columns
H_vs_H <- H_vs_H[, !duplicated(colnames(H_vs_H))]

# Assign sequential IDs by superclass
H_vs_H <- H_vs_H %>%
  arrange(superclass, desc(remddp))

H_vs_H$ID <- NA
counter <- 1

for (sc in unique(H_vs_H$superclass)) {
  group_indices <- which(H_vs_H$superclass == sc)
  H_vs_H$ID[group_indices] <- counter:(counter + length(group_indices) - 1)
  counter <- max(H_vs_H$ID, na.rm = TRUE) + 5  # Add space between groups
}

# Add matching IDs to Na_vs_H based on Index_1
Na_vs_H$ID <- H_vs_H$ID[match(Na_vs_H$Index_1, H_vs_H$Index_1)]

# Add type identifiers and combine datasets
H_vs_H$Type <- "H_vs_H"
Na_vs_H$Type <- "Na_vs_H"
combined_data <- bind_rows(H_vs_H, Na_vs_H)

# Calculate label positions for plot
label_positions <- combined_data %>%
  group_by(superclass) %>%
  summarize(pos = mean(ID, na.rm = TRUE), 
            y_pos = min(remddp, na.rm = TRUE) - 0.1, 
            .groups = "drop")

# Calculate statistics
dot_counts <- combined_data %>%
  group_by(superclass, Type) %>%
  summarize(total_cases = n(),
            below_70 = sum(remddp < 70, na.rm = TRUE),
            percentage_below_70 = round((sum(remddp < 70, na.rm = TRUE) / n()) * 100, 2),
            .groups = "drop")

# Export statistics
write_xlsx(dot_counts, file.path(adduct_dir, "dot_counts.xlsx"))

# Create dot plot
plot <- ggplot(combined_data, aes(x = ID, y = remddp, color = Type)) +
  geom_point(size = 1.5) +
  geom_vline(data = combined_data %>% 
               group_by(superclass) %>% 
               summarize(pos = min(ID), .groups = "drop"), 
             aes(xintercept = pos), linetype = "dashed", color = "black") +
  geom_text(data = label_positions, 
            aes(x = pos, y = y_pos, label = superclass), 
            angle = 0, vjust = 6, hjust = 1, size = 1, inherit.aes = FALSE) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  labs(x = "ID", y = "Average Dot Product", title = "Average Dot Product by Superclass") +
  scale_shape_manual(values = c("H_vs_H" = 16, "Na_vs_H" = 16)) +
  scale_color_manual(values = c("H_vs_H" = "#2667FF", "Na_vs_H" = "#E25D75"))

# Save the plot
ggsave(file.path(adduct_dir, "average_dot_product_plot_no_legend.tiff"), 
       plot = plot, width = 5, height = 3.8, dpi = 150)

# Calculate differences between H_vs_H and Na_vs_H
dp_diff <- H_vs_H %>%
  left_join(Na_vs_H, by = "ID", suffix = c("_H", "_Na")) %>%
  mutate(dp_difference = remddp_H - remddp_Na,
         CE = (CE1_H + CE2_H + CE1_Na + CE2_Na) / 4)

# Define color gradient based on CE values
dp_diff$CE_color <- cut(dp_diff$CE, 
                        breaks = c(-Inf, 10, 30, Inf), 
                        labels = c("pink", "lightcoral", "red"))

# Calculate statistics by CE level
mean_dp_diff <- dp_diff %>%
  group_by(CE_color) %>%
  summarize(mean_dp = mean(dp_difference, na.rm = TRUE),
            num_cases = n())

# Find superclass boundaries
superclass_positions <- dp_diff %>%
  group_by(superclass_H) %>%
  summarize(pos = min(ID), .groups = "drop")

# Create difference plot
diff_plot <- ggplot(dp_diff, aes(x = ID, y = dp_difference, color = CE_color)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(data = superclass_positions, 
             aes(xintercept = pos), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("pink" = "pink", "lightcoral" = "lightcoral", "red" = "red")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  labs(x = "ID", y = "Difference in Average Dot Product", 
       title = "Difference in Average Dot Product by ID")

# Save the difference plot
ggsave(file.path(adduct_dir, "difference_average_dot_product_with_superclass_lines.tiff"), 
       plot = diff_plot, width = 5, height = 3.8, dpi = 150)