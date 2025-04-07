setwd("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20")

library(openxlsx)
adduct_H <- read.xlsx("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/adduct_[M+H]+.xlsx")
adduct_Na <- read.xlsx("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/adduct_[M+Na]+.xlsx")


library(dplyr)

# Filter rows where WHOLE_SMILES is present in both datasets
common_smiles <- intersect(adduct_H$WHOLE_SMILES, adduct_Na$WHOLE_SMILES)

# Subset both data frames
adduct_H_filtered <- adduct_H %>% filter(WHOLE_SMILES %in% common_smiles)
adduct_Na_filtered <- adduct_Na %>% filter(WHOLE_SMILES %in% common_smiles)


library(readxl)

# Load the Excel file
NIST20_Pos_classfied <- read_excel(
  "X:/Users/Botao_Liu/Worklog/Database/NIST20_Pos_classfied.xlsx",
  col_types = "text"
)

# Add new columns with initial NA values
adduct_H_filtered[, c("superclass", "class", "subclass", "canonical")] <- NA

# For each row in adduct_H_filtered
for(i in 1:nrow(adduct_H_filtered)) {
  # Find matching rows in NIST20_Pos_classfied
  matches <- which(NIST20_Pos_classfied$WHOLE_COMPOUND_NAME == adduct_H_filtered$WHOLE_COMPOUND_NAME[i])
  
  # If matches found, use the first match
  if(length(matches) > 0) {
    adduct_H_filtered$WHOLE_SMILES[i] <- NIST20_Pos_classfied$WHOLE_SMILES[matches[1]]
    adduct_H_filtered$superclass[i] <- NIST20_Pos_classfied$superclass[matches[1]]
    adduct_H_filtered$class[i] <- NIST20_Pos_classfied$class[matches[1]]
    adduct_H_filtered$subclass[i] <- NIST20_Pos_classfied$subclass[matches[1]]
    adduct_H_filtered$canonical[i] <- NIST20_Pos_classfied$canonical[matches[1]]
  }
}


# Add new columns with initial NA values
adduct_Na_filtered[, c("superclass", "class", "subclass", "canonical")] <- NA

# For each row in adduct_Na_filtered
for(i in 1:nrow(adduct_Na_filtered)) {
  # Find matching rows in NIST20_Pos_classfied
  matches <- which(NIST20_Pos_classfied$WHOLE_COMPOUND_NAME == adduct_Na_filtered$WHOLE_COMPOUND_NAME[i])
  
  # If matches found, use the first match
  if(length(matches) > 0) {
    adduct_Na_filtered$WHOLE_SMILES[i] <- NIST20_Pos_classfied$WHOLE_SMILES[matches[1]]
    adduct_Na_filtered$superclass[i] <- NIST20_Pos_classfied$superclass[matches[1]]
    adduct_Na_filtered$class[i] <- NIST20_Pos_classfied$class[matches[1]]
    adduct_Na_filtered$subclass[i] <- NIST20_Pos_classfied$subclass[matches[1]]
    adduct_Na_filtered$canonical[i] <- NIST20_Pos_classfied$canonical[matches[1]]
  }
}



# Load necessary library
library(writexl)

# Set the working directory
setwd("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20")

# Write the data frame to an Excel file
write_xlsx(adduct_H_filtered, "adduct_H_filtered.xlsx")

# Write the data frame to an Excel file
write_xlsx(adduct_Na_filtered, "adduct_Na_filtered.xlsx")


library(openxlsx)

adduct_H <- read.xlsx("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/adduct_H_filtered.xlsx")

adduct_Na <- read.xlsx("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/adduct_Na_filtered.xlsx")


library(dplyr)

# Filter superclasses with count above 500
superclass_counts <- adduct_H %>%
  group_by(superclass) %>%
  tally() %>%
  filter(n > 500)

# Get the filtered data with only selected superclasses
adduct_H_filtered <- adduct_H %>%
  filter(superclass %in% superclass_counts$superclass)

# Randomly select 300 rows from each superclass to compose adduct_H_1
set.seed(123)  # Set seed for reproducibility
adduct_H_1 <- adduct_H_filtered %>%
  group_by(superclass) %>%
  sample_n(size = 300, replace = FALSE) %>%
  ungroup()

# Randomly select another 200 rows from each superclass to compose adduct_H_2
set.seed(456)  # Use a different seed
adduct_H_2 <- adduct_H_filtered %>%
  group_by(superclass) %>%
  sample_n(size = 300, replace = FALSE) %>%
  ungroup()

# Save the results
write.xlsx(adduct_H_1, "adduct_H_1.xlsx")
write.xlsx(adduct_H_2, "adduct_H_2.xlsx")

# Load the readxl package for reading the Excel files
library(readxl)

# Set the working directory
setwd("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20")

# Read the Excel files back into R
adduct_H_1 <- read_excel("adduct_H_1.xlsx")
adduct_H_2 <- read_excel("adduct_H_2.xlsx")


library(dplyr)     # For data manipulation
library(data.table) # For fast data manipulation

# Define adduct name
adduct <- "[M+H]+"

# Load and filter data for [M+H]+
adduct_H <- adduct_H_1

# Load and filter data for the current adduct
adduct_other <- adduct_H_2

# Preallocate a list to store pairs
pairs_list <- vector("list", nrow(adduct_H) * nrow(adduct_other))
count <- 0

# Loop through each metabolite from adduct_H
for (i in 1:nrow(adduct_H)) {
  
  for (j in seq_len(nrow(adduct_other))) {
    
    # Check the CE difference
    ce_diff <- abs(adduct_H$WHOLE_COLLISION_ENERGY[i] - adduct_other$WHOLE_COLLISION_ENERGY[j])
    
    # Check if the CE difference is smaller than 3 
    if (ce_diff < 3 ) {
      count <- count + 1
      pairs_list[[count]] <- data.frame(
        Metabolite_Name_1 = adduct_H$WHOLE_COMPOUND_NAME[i],
        Metabolite_Name_2 = adduct_other$WHOLE_COMPOUND_NAME[j],
        SMILE_1 = adduct_H$WHOLE_SMILES[i],
        SMILE_2 = adduct_other$WHOLE_SMILES[j],
        CE1 = adduct_H$WHOLE_COLLISION_ENERGY[i],
        CE2 = adduct_other$WHOLE_COLLISION_ENERGY[j],
        Adduct1 = adduct,
        Adduct2 = adduct,
        Index_1 = adduct_H$Index[i],
        Index_2 = adduct_other$Index[j],
        precursor1 = adduct_H$WHOLE_PRECURSOR_mz[i],
        precursor2 = adduct_other$WHOLE_PRECURSOR_mz[j],
        whole_mz1 = adduct_H$WHOLE_mz[i],
        whole_mz2 = adduct_other$WHOLE_mz[j],
        whole_int1 = adduct_H$WHOLE_Int[i],
        whole_int2 = adduct_other$WHOLE_Int[j],
        whole_mz_cleaned_1 = adduct_H$WHOLE_mz_cleaned[i],
        whole_mz_cleaned_2 = adduct_other$WHOLE_mz_cleaned[j],
        whole_int_cleaned_1 = adduct_H$WHOLE_Int_cleaned[i],
        whole_int_cleaned_2 = adduct_other$WHOLE_Int_cleaned[j],
        INSTRUMENT_TYPE_1 = adduct_H$INSTRUMENT_TYPE[i],
        INSTRUMENT_TYPE_2 = adduct_other$INSTRUMENT_TYPE[j]
      )
    }
  }
}

# Combine all pairs into a data frame
if (count > 0) {
  pairs_df <- do.call(rbind, pairs_list[1:count]) # Only take the filled parts
  print(paste("Pairs found for adduct:", adduct, ":", nrow(pairs_df)))
  
  # Write the dataframe to a CSV file
  write.csv(pairs_df, paste0(adduct, "_pairs_df.csv"), row.names = FALSE)
} else {
  print(paste("No pairs found for adduct:", adduct))
}

pairs_df$precursor1 <- as.numeric(pairs_df$precursor1)
pairs_df$precursor2 <- as.numeric(pairs_df$precursor2)

# Load required libraries for parallel processing
library(parallel)
library(doParallel)

process_spectral_matching <- function(pairs_df, noise_threshold = 5, mz_tol = 0.02) {
  # Detect the number of cores and set up cluster
  library(parallel)
  library(doParallel)
  library(dplyr)
  
  num_cores <- detectCores() - 1  # Leave one core free for system processes
  if (num_cores < 1) num_cores <- 1
  
  cat(paste("Setting up parallel processing with", num_cores, "cores\n"))
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Function to calculate dot product and modified cosine similarity
  calculate_dot_product_and_modified_cosine <- function(ref_fragMz, ref_fragInt, exp_fragMz, exp_fragInt, precursor_mz_diff, mz_tol = 0.02) {
    # Pre-allocate vectors and variables
    dot_exp_int_original <- numeric(length(ref_fragMz))
    dot_exp_int_modified <- numeric(length(ref_fragMz))
    used_indices <- logical(length(exp_fragMz))
    
    
    original_aligned_count <- 0
    modified_aligned_count <- 0
    
    # Pre-calculate distance matrices
    dist_matrix_original <- outer(exp_fragMz, ref_fragMz, function(x, y) abs(x - y))
    dist_matrix_modified <- outer(exp_fragMz, ref_fragMz + precursor_mz_diff, function(x, y) abs(x - y))
    
    # Perform alignment and record results
    for (j in seq_along(ref_fragMz)) {
      # Check for original matches
      fragMz_indices <- which(dist_matrix_original[, j] <= mz_tol & !used_indices)
      
      if (length(fragMz_indices) > 0) {
        selected_index <- fragMz_indices[which.max(exp_fragInt[fragMz_indices])]
        dot_exp_int_original[j] <- exp_fragInt[selected_index]
        used_indices[selected_index] <- TRUE
        original_aligned_count <- original_aligned_count + 1
      }
      
      # Check for modified matches
      modified_indices <- which(dist_matrix_modified[, j] <= mz_tol & !used_indices)
      
      if (length(modified_indices) > 0) {
        selected_index <- modified_indices[which.max(exp_fragInt[modified_indices])]
        dot_exp_int_modified[j] <- exp_fragInt[selected_index]
        used_indices[selected_index] <- TRUE
        modified_aligned_count <- modified_aligned_count + 1
      }
    }
    
    # Combine both original and modified matched intensities for modified cosine similarity
    combined_exp_int <- pmax(dot_exp_int_original, dot_exp_int_modified)
    
    # Calculate norms
    norm_ref <- sqrt(sum(ref_fragInt^2))
    norm_exp <- sqrt(sum(exp_fragInt^2))
    
    # Calculate dot products
    dot_product <- if (norm_ref == 0) 0 else 100 * sum(ref_fragInt * dot_exp_int_original) / (norm_ref * norm_exp)
    modified_cosine_similarity <- if (norm_ref == 0) 0 else 100 * sum(ref_fragInt * combined_exp_int) / (norm_ref * norm_exp)
    
    total_aligned_fragments <- original_aligned_count + modified_aligned_count
    
    return(list(
      dot_product = dot_product,
      modified_cosine_similarity = modified_cosine_similarity,
      original_aligned_fragments = original_aligned_count,
      modified_aligned_fragments = modified_aligned_count,
      total_aligned_fragments = total_aligned_fragments
    ))
  }
  
  # Completely define the row processing function within the main function
  process_row <- function(i, pairs_data, noise_threshold, mz_tol) {
    if (i %% 100 == 0) cat(paste("Processing row", i, "\n"))
    
    # Extract precursor m/z
    precursor1 <- pairs_data$precursor1[i]
    precursor2 <- pairs_data$precursor2[i]
    
    # Extract and clean fragment m/z and intensities
    whole_mz_cleaned_1 <- as.numeric(unlist(strsplit(as.character(pairs_data$whole_mz1[i]), ",")))
    whole_int_cleaned_1 <- as.numeric(unlist(strsplit(as.character(pairs_data$whole_int1[i]), ",")))
    whole_mz_cleaned_2 <- as.numeric(unlist(strsplit(as.character(pairs_data$whole_mz2[i]), ",")))
    whole_int_cleaned_2 <- as.numeric(unlist(strsplit(as.character(pairs_data$whole_int2[i]), ",")))
    
    # Filter intensities above noise threshold
    valid_indices_1 <- whole_int_cleaned_1 >= noise_threshold
    valid_indices_2 <- whole_int_cleaned_2 >= noise_threshold
    
    # Clean the m/z and intensity values based on noise threshold
    whole_mz_cleaned_1 <- whole_mz_cleaned_1[valid_indices_1]
    whole_int_cleaned_1 <- whole_int_cleaned_1[valid_indices_1]
    whole_mz_cleaned_2 <- whole_mz_cleaned_2[valid_indices_2]
    whole_int_cleaned_2 <- whole_int_cleaned_2[valid_indices_2]
    
    # Calculate lower precursor m/z for thresholding
    lower_precursor_mz <- min(precursor1, precursor2)
    
    # For dot product calculation
    mz_threshold_dot <- lower_precursor_mz - 0.02
    ref_filter_dot <- whole_mz_cleaned_2 <= mz_threshold_dot
    exp_filter_dot <- whole_mz_cleaned_1 <= mz_threshold_dot
    
    ref_fragMz_dot <- whole_mz_cleaned_2[ref_filter_dot]
    ref_fragInt_dot <- whole_int_cleaned_2[ref_filter_dot]
    exp_fragMz_dot <- whole_mz_cleaned_1[exp_filter_dot]
    exp_fragInt_dot <- whole_int_cleaned_1[exp_filter_dot]
    
    # For fragment alignment calculation
    mz_threshold_align <- lower_precursor_mz + 0.02
    ref_filter_align <- whole_mz_cleaned_2 <= mz_threshold_align
    exp_filter_align <- whole_mz_cleaned_1 <= mz_threshold_align
    
    Frag_Number_MH <- sum(exp_filter_align)
    Frag_Number_Na <- sum(ref_filter_align)
    
    # Calculate precursor m/z difference
    precursor_mz_diff <- precursor2 - precursor1
    
    # Calculate using filtered fragments
    result <- calculate_dot_product_and_modified_cosine(
      ref_fragMz = ref_fragMz_dot, 
      ref_fragInt = ref_fragInt_dot,
      exp_fragMz = exp_fragMz_dot, 
      exp_fragInt = exp_fragInt_dot,
      precursor_mz_diff = precursor_mz_diff,
      mz_tol = mz_tol
    )
    
    # Calculate derived metrics
    total_aligned <- result$original_aligned_fragments + result$modified_aligned_fragments
    total_fragments <- Frag_Number_MH + Frag_Number_Na
    
    Direct_match <- ifelse(total_aligned == 0, 0, 100 * result$original_aligned_fragments / total_aligned)
    Shifted_match <- ifelse(total_aligned == 0, 0, 100 * result$modified_aligned_fragments / total_aligned)
    
    FP <- ifelse(total_fragments == 0, 0, 100 * 2 * total_aligned / total_fragments)
    
    ts_dp <- (result$dot_product + FP) / 2
    ts_remddp <- (result$modified_cosine_similarity + FP) / 2
    
    # Return results as a list
    return(list(
      dp = round(result$dot_product),
      remddp = round(result$modified_cosine_similarity),
      dp_aligned_fragments = round(result$original_aligned_fragments),
      remddp_aligned_fragments = round(result$modified_aligned_fragments),
      Frag_Number_MH = Frag_Number_MH,
      Frag_Number_Na = Frag_Number_Na,
      total_aligned = round(total_aligned),
      Direct_match = round(Direct_match),
      Shifted_match = round(Shifted_match),
      total_fragments = round(total_fragments),
      FP = round(FP),
      ts_dp = round(ts_dp),
      ts_remddp = round(ts_remddp)
    ))
  }
  
  # Export all required functions and libraries to the cluster
  clusterEvalQ(cl, {
    library(dplyr)
  })
  
  # Process in parallel
  start_time <- Sys.time()
  cat("Starting parallel processing...\n")
  
  # Process all rows
  results <- parLapply(cl, 1:nrow(pairs_df), function(i) {
    process_row(i, pairs_df, noise_threshold, mz_tol)
  })
  
  # Stop the cluster
  stopCluster(cl)
  
  end_time <- Sys.time()
  cat(paste("Processing completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n"))
  
  # Add results to the dataframe
  pairs_df <- pairs_df %>%
    mutate(
      dp = sapply(results, `[[`, "dp"),
      remddp = sapply(results, `[[`, "remddp"),
      dp_aligned_fragments = sapply(results, `[[`, "dp_aligned_fragments"),
      remddp_aligned_fragments = sapply(results, `[[`, "remddp_aligned_fragments"),
      Frag_Number_MH = sapply(results, `[[`, "Frag_Number_MH"),
      Frag_Number_Na = sapply(results, `[[`, "Frag_Number_Na"),
      total_aligned = sapply(results, `[[`, "total_aligned"),
      Direct_match = sapply(results, `[[`, "Direct_match"),
      Shifted_match = sapply(results, `[[`, "Shifted_match"),
      total_fragments = sapply(results, `[[`, "total_fragments"),
      FP = sapply(results, `[[`, "FP"),
      ts_dp = sapply(results, `[[`, "ts_dp"),
      ts_remddp = sapply(results, `[[`, "ts_remddp")
    )
  
  return(pairs_df)
}

# Prepare and process function 
prepare_and_process <- function(pairs_df, noise_threshold = 5, mz_tol = 0.02) {
  # Initialize result columns if they don't exist
  result_columns <- c("dp", "remddp", "dp_aligned_fragments", "remddp_aligned_fragments",
                      "Frag_Number_MH", "Frag_Number_Na", "total_aligned", "Direct_match", 
                      "Shifted_match", "total_fragments", "FP", "ts_dp", "ts_remddp")
  
  for (col in result_columns) {
    if (!(col %in% names(pairs_df))) {
      pairs_df[[col]] <- numeric(nrow(pairs_df))
    }
  }
  
  # Run the processing
  return(process_spectral_matching(pairs_df, noise_threshold, mz_tol))
}

# Run the parallelized processing
pairs_df <- prepare_and_process(pairs_df)

# Display the results
cat("Processing complete. Here's a sample of the results:\n")
print(head(pairs_df))
library(tidyr)
# Replace NA values with 0 in all numeric columns
pairs_df <- pairs_df %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))







# Set the working directory
setwd("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20")

write.xlsx(pairs_df, "pairs_df_H_vs_H_dp.xlsx")


##############################
# fingerprint-based structural similarity
##############################
adduct <- "pairs_df_H_vs_H_dp"

library(openxlsx)
library(rJava)
library(rcdk)
library(fingerprint)

# Set the working directory
setwd("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20")

# Read the data
pairs_df <- read.xlsx("pairs_df_H_vs_H_dp.xlsx")

pairs_df <- pairs_df[ , c("SMILE_1", "SMILE_2")]

# Split the data into 15 parts
split_size <- ceiling(nrow(pairs_df) / 15)
split_list <- split(pairs_df, ceiling(seq_len(nrow(pairs_df)) / split_size))


library(openxlsx)
library(rcdk)
library(fingerprint)
library(parallel)
library(doParallel)

# Determine number of cores to use (max - 1)
num_cores <- detectCores() - 1
cat(paste("Using", num_cores, "cores for parallel processing\n"))

# Structural Similarity Function
structureSmi <- function(smiles_A, smiles_B) {
  if (is.na(smiles_A) || is.na(smiles_B) || smiles_A == "" || smiles_B == "") {
    return(NA)  # Skip if either SMILES is missing or empty
  }
  
  tryCatch({
    mols <- parse.smiles(c(smiles_A, smiles_B))
    
    if (any(sapply(mols, is.null))) {
      return(NA)  # Skip if parsing fails
    }
    
    fps <- lapply(mols, get.fingerprint, type = "pubchem")
    fp.sim.matrix(fps, method = 'tanimoto')[2, 1]
  }, error = function(e) {
    return(NA)  # Skip if any unexpected error occurs
  })
}

# Main Processing Function
process_files <- function(num_parts = 15) {
  # Initialize an empty list to store cleaned data frames
  cleaned_list <- list()
  
  # Process Each File One by One
  for (i in 1:num_parts) {
    cat(paste("Processing Part", i, "\n"))
    
    # Read Excel file
    part_df <- read.xlsx(paste0("pairs_df_part_", i, ".xlsx"))
    
    # Set up parallel processing
    cl <- makeCluster(num_cores)
    
    # Export necessary functions and libraries to all worker nodes
    clusterExport(cl, "structureSmi")
    clusterEvalQ(cl, {
      library(rcdk)
      library(fingerprint)
    })
    
    # Process in parallel using clusterApply
    # Split data into chunks
    chunk_size <- ceiling(nrow(part_df) / num_cores)
    row_indices <- split(1:nrow(part_df), ceiling(seq_along(1:nrow(part_df)) / chunk_size))
    
    # Function to process a chunk
    process_chunk <- function(indices) {
      chunk <- part_df[indices, ]
      chunk$Structural_Similarity <- mapply(structureSmi, 
                                            chunk$SMILE_1, 
                                            chunk$SMILE_2)
      # Return only non-NA rows
      return(chunk[!is.na(chunk$Structural_Similarity), ])
    }
    
    # Process chunks in parallel
    processed_chunks <- clusterApply(cl, row_indices, process_chunk)
    
    # Stop cluster
    stopCluster(cl)
    
    # Combine processed chunks
    part_df_cleaned <- do.call(rbind, processed_chunks)
    
    # Save cleaned dataset to list
    cleaned_list[[i]] <- part_df_cleaned
    
    # Save the cleaned part separately 
    write.xlsx(part_df_cleaned, 
               paste0("pairs_df_part_", i, "_with_similarity_cleaned.xlsx"), 
               row.names = FALSE)
    
    # Garbage collection to manage memory
    gc()
  }
  
  # Merge all parts into one dataset
  merged_df <- do.call(rbind, cleaned_list)
  
  # Save the merged dataset to an Excel file
  write.xlsx(merged_df, 
             "pairs_df_merged_with_similarity.xlsx", 
             row.names = FALSE)
  
  cat("All parts have been processed, cleaned, and merged successfully.\n")
  
  return(merged_df)
}

# Execute the processing
merged_data <- process_files()

library(openxlsx)

merged_data <- read.xlsx("pairs_df_merged_with_similarity.xlsx")


library(readxl)

# Load the Excel file
NIST20_Pos_classfied <- read_excel(
  "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/NIST20_Pos_classfied.xlsx",
  col_types = "text"
)
NIST20_Pos_classfied <- NIST20_Pos_classfied[, c("WHOLE_SMILES", "superclass", "class", "subclass")]

# Remove duplicated WHOLE_SMILES in NIST20_Pos_classfied, keeping only the first match
NIST20_Pos_classfied_unique <- NIST20_Pos_classfied[!duplicated(NIST20_Pos_classfied$WHOLE_SMILES), ]

# Merge the two data frames based on SMILES
merged_data <- merge(merged_data, NIST20_Pos_classfied_unique[, c("WHOLE_SMILES", "superclass", "class", "subclass")], 
                     by.x = "SMILE_1", by.y = "WHOLE_SMILES", all.x = TRUE)

# Reorder rows to the original order if needed
merged_data <- merged_data[order(as.numeric(row.names(merged_data))), ]


# Load necessary library
library(writexl)

# Set the working directory
setwd("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20")

# Write the data frame to an Excel file
write_xlsx(merged_data, "pairs_df_H_vs_H_dp_with_similarity_classed.xlsx")









####Mannually merge together

merged_data <- read.xlsx("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/pairs_df_H_vs_H_all_together.xlsx")

merged_data <- merged_data %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))


library(ggplot2)

# Plot distribution colored by superclass
ggplot(merged_data, aes(x = remddp, y = Structural_Similarity, color = superclass)) +
  geom_point(alpha = 0.6) + # Transparent scatter points
  labs(title = "Distribution of Structural Similarity vs. remddp",
       x = "remddp",
       y = "Structural Similarity",
       color = "Superclass") + # Legend title
  theme_minimal() + # Clean theme
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) 
 scale_color_manual(values = c("#E25D75", "#2667FF", "green", "orange", "purple", "black", "grey", "pink")) # Custom colors


# Filter data where both remddp 
filtered_pairs_df <- merged_data[merged_data$remddp > 70 & merged_data$Structural_Similarity > 0.7, ]

# Filter for identical instrument type pairs
filtered_pairs_df <- filtered_pairs_df %>%
  filter(INSTRUMENT_TYPE_1 == INSTRUMENT_TYPE_2)

# Add a 'highlight' column to label the points as 'TP', 'FP', or 'Other'
merged_data$highlight <- ifelse(merged_data$remddp > 70 & merged_data$Structural_Similarity > 0.7, "TP", 
                                ifelse(merged_data$remddp > 70 & merged_data$Structural_Similarity <= 0.7, "FP", "Other"))

# Plot the distribution with TP in blue, FP in red, and other points in grey, removing grid, legend, and title
ggplot(merged_data, aes(x = remddp, y = Structural_Similarity, color = highlight)) +
  geom_point(alpha = 0.6) + # Transparent scatter points
  scale_color_manual(values = c("TP" = "#2667FF", "FP" = "grey", "Other" = "grey")) + # Custom colors for TP, FP, and others
  theme_minimal() + # Clean theme
  theme(
    plot.title = element_blank(), # Remove the title
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none", # Remove legend
    panel.grid = element_blank() # Remove grid
  ) +
  coord_flip() + # Flip the axes
  ylim(0, 1) + # Setting limits for better visualization
  xlim(0, 100) # Setting limits for better visualization

# Save the plot as a TIFF file
ggsave("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/H_vs_H70.tiff", width = 8, height = 8, dpi = 300)

# Add a 'highlight' column to label the points as 'TP', 'FP', or 'Other'
merged_data$highlight <- ifelse(merged_data$remddp > 70 & merged_data$Structural_Similarity > 0.7, "TP", 
                                ifelse(merged_data$remddp > 70 & merged_data$Structural_Similarity <= 0.7, "FP", "Other"))

# Plot the distribution with TP in blue, FP in red, and other points in grey, removing grid, legend, and title
ggplot(merged_data, aes(x = remddp, y = Structural_Similarity, color = highlight)) +
  geom_point(alpha = 0.6) + # Transparent scatter points
  scale_color_manual(values = c("TP" = "#2667FF", "FP" = "grey", "Other" = "grey")) + # Custom colors for TP, FP, and others
  theme_minimal() + # Clean theme
  theme(
    plot.title = element_blank(), # Remove the title
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none", # Remove legend
    panel.grid = element_blank() # Remove grid
  ) +
  coord_flip() + # Flip the axes
  ylim(0, 1) + # Setting limits for better visualization
  xlim(0, 100) # Setting limits for better visualization



setwd("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20")

library(openxlsx)
adduct_Na <- read.xlsx("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/adduct_Na_filtered.xlsx")

H_vs_H <- filtered_pairs_df

# Add Index column to both filtered_pairs_df and H_vs_H
filtered_pairs_df$Index <- 1:nrow(filtered_pairs_df)
H_vs_H$Index <- 1:nrow(H_vs_H)

# Check if row counts match
if (nrow(filtered_pairs_df) != nrow(H_vs_H)) {
  stop("The number of rows in filtered_pairs_df and H_vs_H do not match!")
}

# Iterate over each row in filtered_pairs_df
for (i in 1:nrow(filtered_pairs_df)) {
  # Extract Metabolite_Name and CE2 for the current row in filtered_pairs_df
  metabolite_name <- filtered_pairs_df$Metabolite_Name_2[i]
  ce2 <- filtered_pairs_df$CE2[i]
  
  # Find the matching rows in adduct_Na based on Metabolite_Name
  matched_rows <- adduct_Na[adduct_Na$WHOLE_COMPOUND_NAME == metabolite_name, ]
  
  # If there are matched rows, proceed
  if (nrow(matched_rows) > 0) {
    
    # Find the closest WHOLE_COLLISION_ENERGY to CE2 in matched rows
    closest_match <- matched_rows[which.min(abs(matched_rows$WHOLE_COLLISION_ENERGY - ce2)), ]
    
    # Update filtered_pairs_df with the matched values from adduct_Na
    filtered_pairs_df$Metabolite_Name_2[i] <- closest_match$WHOLE_COMPOUND_NAME
    filtered_pairs_df$SMILE_2[i] <- closest_match$WHOLE_SMILES
    filtered_pairs_df$CE2[i] <- closest_match$WHOLE_COLLISION_ENERGY
    filtered_pairs_df$Adduct2[i] <- "[M+Na]+"
    filtered_pairs_df$Index_2[i] <- closest_match$WHOLE_Index
    filtered_pairs_df$precursor2[i] <- closest_match$WHOLE_PRECURSOR_mz
    filtered_pairs_df$whole_mz2[i] <- closest_match$WHOLE_mz
    filtered_pairs_df$whole_int2[i] <- closest_match$WHOLE_Int
    filtered_pairs_df$whole_mz_cleaned_2[i] <- closest_match$WHOLE_mz_cleaned
    filtered_pairs_df$whole_int_cleaned_2[i] <- closest_match$WHOLE_Int_cleaned
    filtered_pairs_df$INSTRUMENT_TYPE_2[i] <- closest_match$INSTRUMENT        
  }
}

# View the updated filtered_pairs_df
head(filtered_pairs_df)

# Find the rows where Adduct2 is "[M+H]+"
rows_to_remove <- which(filtered_pairs_df$Adduct2 == "[M+H]+")

# Check if rows_to_remove are valid (i.e., no empty indices)
if (length(rows_to_remove) == 0) {
  stop("No rows found with Adduct2 == '[M+H]+'")
}

# Remove the corresponding rows from filtered_pairs_df and H_vs_H based on Index
filtered_pairs_df <- filtered_pairs_df[!filtered_pairs_df$Index %in% filtered_pairs_df$Index[rows_to_remove], ]
H_vs_H <- H_vs_H[!H_vs_H$Index %in% filtered_pairs_df$Index[rows_to_remove], ]

# View the updated filtered_pairs_df and H_vs_H
head(filtered_pairs_df)
head(H_vs_H)




pairs_df <- filtered_pairs_df


pairs_df$precursor1 <- as.numeric(pairs_df$precursor1)
pairs_df$precursor2 <- as.numeric(pairs_df$precursor2)


process_spectral_matching <- function(pairs_df, noise_threshold = 5, mz_tol = 0.02) {
  # Detect the number of cores and set up cluster
  library(parallel)
  library(doParallel)
  library(dplyr)
  
  num_cores <- detectCores() - 1  # Leave one core free for system processes
  if (num_cores < 1) num_cores <- 1
  
  cat(paste("Setting up parallel processing with", num_cores, "cores\n"))
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Function to calculate dot product and modified cosine similarity
  calculate_dot_product_and_modified_cosine <- function(ref_fragMz, ref_fragInt, exp_fragMz, exp_fragInt, precursor_mz_diff, mz_tol = 0.02) {
    # Pre-allocate vectors and variables
    dot_exp_int_original <- numeric(length(ref_fragMz))
    dot_exp_int_modified <- numeric(length(ref_fragMz))
    used_indices <- logical(length(exp_fragMz))
    
    
    original_aligned_count <- 0
    modified_aligned_count <- 0
    
    # Pre-calculate distance matrices
    dist_matrix_original <- outer(exp_fragMz, ref_fragMz, function(x, y) abs(x - y))
    dist_matrix_modified <- outer(exp_fragMz, ref_fragMz + precursor_mz_diff, function(x, y) abs(x - y))
    
    # Perform alignment and record results
    for (j in seq_along(ref_fragMz)) {
      # Check for original matches
      fragMz_indices <- which(dist_matrix_original[, j] <= mz_tol & !used_indices)
      
      if (length(fragMz_indices) > 0) {
        selected_index <- fragMz_indices[which.max(exp_fragInt[fragMz_indices])]
        dot_exp_int_original[j] <- exp_fragInt[selected_index]
        used_indices[selected_index] <- TRUE
        original_aligned_count <- original_aligned_count + 1
      }
      
      # Check for modified matches
      modified_indices <- which(dist_matrix_modified[, j] <= mz_tol & !used_indices)
      
      if (length(modified_indices) > 0) {
        selected_index <- modified_indices[which.max(exp_fragInt[modified_indices])]
        dot_exp_int_modified[j] <- exp_fragInt[selected_index]
        used_indices[selected_index] <- TRUE
        modified_aligned_count <- modified_aligned_count + 1
      }
    }
    
    # Combine both original and modified matched intensities for modified cosine similarity
    combined_exp_int <- pmax(dot_exp_int_original, dot_exp_int_modified)
    
    # Calculate norms
    norm_ref <- sqrt(sum(ref_fragInt^2))
    norm_exp <- sqrt(sum(exp_fragInt^2))
    
    # Calculate dot products
    dot_product <- if (norm_ref == 0) 0 else 100 * sum(ref_fragInt * dot_exp_int_original) / (norm_ref * norm_exp)
    modified_cosine_similarity <- if (norm_ref == 0) 0 else 100 * sum(ref_fragInt * combined_exp_int) / (norm_ref * norm_exp)
    
    total_aligned_fragments <- original_aligned_count + modified_aligned_count
    
    return(list(
      dot_product = dot_product,
      modified_cosine_similarity = modified_cosine_similarity,
      original_aligned_fragments = original_aligned_count,
      modified_aligned_fragments = modified_aligned_count,
      total_aligned_fragments = total_aligned_fragments
    ))
  }
  
  # Completely define the row processing function within the main function
  process_row <- function(i, pairs_data, noise_threshold, mz_tol) {
    if (i %% 100 == 0) cat(paste("Processing row", i, "\n"))
    
    # Extract precursor m/z
    precursor1 <- pairs_data$precursor1[i]
    precursor2 <- pairs_data$precursor2[i]
    
    # Extract and clean fragment m/z and intensities
    whole_mz_cleaned_1 <- as.numeric(unlist(strsplit(as.character(pairs_data$whole_mz1[i]), ",")))
    whole_int_cleaned_1 <- as.numeric(unlist(strsplit(as.character(pairs_data$whole_int1[i]), ",")))
    whole_mz_cleaned_2 <- as.numeric(unlist(strsplit(as.character(pairs_data$whole_mz2[i]), ",")))
    whole_int_cleaned_2 <- as.numeric(unlist(strsplit(as.character(pairs_data$whole_int2[i]), ",")))
    
    # Filter intensities above noise threshold
    valid_indices_1 <- whole_int_cleaned_1 >= noise_threshold
    valid_indices_2 <- whole_int_cleaned_2 >= noise_threshold
    
    # Clean the m/z and intensity values based on noise threshold
    whole_mz_cleaned_1 <- whole_mz_cleaned_1[valid_indices_1]
    whole_int_cleaned_1 <- whole_int_cleaned_1[valid_indices_1]
    whole_mz_cleaned_2 <- whole_mz_cleaned_2[valid_indices_2]
    whole_int_cleaned_2 <- whole_int_cleaned_2[valid_indices_2]
    
    # Calculate lower precursor m/z for thresholding
    lower_precursor_mz <- min(precursor1, precursor2)
    
    # For dot product calculation
    mz_threshold_dot <- lower_precursor_mz - 0.02
    ref_filter_dot <- whole_mz_cleaned_2 <= mz_threshold_dot
    exp_filter_dot <- whole_mz_cleaned_1 <= mz_threshold_dot
    
    ref_fragMz_dot <- whole_mz_cleaned_2[ref_filter_dot]
    ref_fragInt_dot <- whole_int_cleaned_2[ref_filter_dot]
    exp_fragMz_dot <- whole_mz_cleaned_1[exp_filter_dot]
    exp_fragInt_dot <- whole_int_cleaned_1[exp_filter_dot]
    
    # For fragment alignment calculation
    mz_threshold_align <- lower_precursor_mz + 0.02
    ref_filter_align <- whole_mz_cleaned_2 <= mz_threshold_align
    exp_filter_align <- whole_mz_cleaned_1 <= mz_threshold_align
    
    Frag_Number_MH <- sum(exp_filter_align)
    Frag_Number_Na <- sum(ref_filter_align)
    
    # Calculate precursor m/z difference
    precursor_mz_diff <- precursor2 - precursor1
    
    # Calculate using filtered fragments
    result <- calculate_dot_product_and_modified_cosine(
      ref_fragMz = ref_fragMz_dot, 
      ref_fragInt = ref_fragInt_dot,
      exp_fragMz = exp_fragMz_dot, 
      exp_fragInt = exp_fragInt_dot,
      precursor_mz_diff = precursor_mz_diff,
      mz_tol = mz_tol
    )
    
    # Calculate derived metrics
    total_aligned <- result$original_aligned_fragments + result$modified_aligned_fragments
    total_fragments <- Frag_Number_MH + Frag_Number_Na
    
    Direct_match <- ifelse(total_aligned == 0, 0, 100 * result$original_aligned_fragments / total_aligned)
    Shifted_match <- ifelse(total_aligned == 0, 0, 100 * result$modified_aligned_fragments / total_aligned)
    
    FP <- ifelse(total_fragments == 0, 0, 100 * 2 * total_aligned / total_fragments)
    
    ts_dp <- (result$dot_product + FP) / 2
    ts_remddp <- (result$modified_cosine_similarity + FP) / 2
    
    # Return results as a list
    return(list(
      dp = round(result$dot_product),
      remddp = round(result$modified_cosine_similarity),
      dp_aligned_fragments = round(result$original_aligned_fragments),
      remddp_aligned_fragments = round(result$modified_aligned_fragments),
      Frag_Number_MH = Frag_Number_MH,
      Frag_Number_Na = Frag_Number_Na,
      total_aligned = round(total_aligned),
      Direct_match = round(Direct_match),
      Shifted_match = round(Shifted_match),
      total_fragments = round(total_fragments),
      FP = round(FP),
      ts_dp = round(ts_dp),
      ts_remddp = round(ts_remddp)
    ))
  }
  
  # Export all required functions and libraries to the cluster
  clusterEvalQ(cl, {
    library(dplyr)
  })
  
  # Process in parallel
  start_time <- Sys.time()
  cat("Starting parallel processing...\n")
  
  # Process all rows
  results <- parLapply(cl, 1:nrow(pairs_df), function(i) {
    process_row(i, pairs_df, noise_threshold, mz_tol)
  })
  
  # Stop the cluster
  stopCluster(cl)
  
  end_time <- Sys.time()
  cat(paste("Processing completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n"))
  
  # Add results to the dataframe
  pairs_df <- pairs_df %>%
    mutate(
      dp = sapply(results, `[[`, "dp"),
      remddp = sapply(results, `[[`, "remddp"),
      dp_aligned_fragments = sapply(results, `[[`, "dp_aligned_fragments"),
      remddp_aligned_fragments = sapply(results, `[[`, "remddp_aligned_fragments"),
      Frag_Number_MH = sapply(results, `[[`, "Frag_Number_MH"),
      Frag_Number_Na = sapply(results, `[[`, "Frag_Number_Na"),
      total_aligned = sapply(results, `[[`, "total_aligned"),
      Direct_match = sapply(results, `[[`, "Direct_match"),
      Shifted_match = sapply(results, `[[`, "Shifted_match"),
      total_fragments = sapply(results, `[[`, "total_fragments"),
      FP = sapply(results, `[[`, "FP"),
      ts_dp = sapply(results, `[[`, "ts_dp"),
      ts_remddp = sapply(results, `[[`, "ts_remddp")
    )
  
  return(pairs_df)
}

# Prepare and process function 
prepare_and_process <- function(pairs_df, noise_threshold = 5, mz_tol = 0.02) {
  # Initialize result columns if they don't exist
  result_columns <- c("dp", "remddp", "dp_aligned_fragments", "remddp_aligned_fragments",
                      "Frag_Number_MH", "Frag_Number_Na", "total_aligned", "Direct_match", 
                      "Shifted_match", "total_fragments", "FP", "ts_dp", "ts_remddp")
  
  for (col in result_columns) {
    if (!(col %in% names(pairs_df))) {
      pairs_df[[col]] <- numeric(nrow(pairs_df))
    }
  }
  
  # Run the processing
  return(process_spectral_matching(pairs_df, noise_threshold, mz_tol))
}

# Run the parallelized processing
pairs_df <- prepare_and_process(pairs_df)

# Display the results
cat("Processing complete. Here's a sample of the results:\n")
print(head(pairs_df))

# Replace NA values with 0 in all numeric columns
pairs_df <- pairs_df %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))

H_vs_Na <- pairs_df
# Add a new column to each data frame to indicate the source of the data
H_vs_H$Source <- "H_vs_H"
H_vs_Na$Source <- "H_vs_Na"

# Combine the two data frames into one
combined_data <- rbind(H_vs_H, H_vs_Na)

# Plot distribution with x and y axes flipped, removing grid, adding lines at x=0.7 and y=0.7, and removing the legend and title
plot <- ggplot(combined_data, aes(x = Structural_Similarity, y = remddp, color = Source)) +
  geom_point(alpha = 0.6) + # Transparent scatter points
  scale_color_manual(values = c("H_vs_H" = "#2667FF", "H_vs_Na" = "#E25D75")) + # Custom colors
  theme_minimal() + # Clean theme
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid = element_blank(),  # Remove grid
        legend.position = "none") +   # Remove legend
  xlim(0.69, 1.0) + ylim(0, 100) 

plot

# Save the plot as a TIFF file
ggsave("Na_demage.tiff", plot = plot, device = "tiff", dpi = 100, width = 4, height = 4)




# Ensure necessary libraries are loaded
library(dplyr)
library(ggplot2)

# Add a new column to each data frame to indicate the source of the data
H_vs_H$Source <- "H_vs_H"
H_vs_Na$Source <- "H_vs_Na"

# Sort H_vs_H by remddp
H_vs_H_sorted <- H_vs_H %>%
  filter(Metabolite_Name_1 != Metabolite_Name_2) %>%
  arrange(desc(remddp))

# Find common indices between H_vs_H and H_vs_Na
common_indices <- intersect(H_vs_H_sorted$Index, H_vs_Na$Index)

# Filter both data frames to retain only rows with common indices
H_vs_H_filtered <- H_vs_H_sorted %>% filter(Index %in% common_indices)
H_vs_Na_filtered <- H_vs_Na %>% filter(Index %in% common_indices)

# Combine the two data frames into one
combined_data <- rbind(H_vs_H_filtered, H_vs_Na_filtered)

# Plot distribution with x and y axes flipped, removing grid, adding lines at x=0.7 and y=0.7, and removing the legend and title
plot <- ggplot(combined_data, aes(x = Structural_Similarity, y = remddp, color = Source)) +
  geom_point(size = 2)  +
  scale_color_manual(values = c("H_vs_H" = "#2667FF", "H_vs_Na" = "#E25D75")) +  # Custom colors
  theme_minimal() +  # Clean theme
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid = element_blank(),  # Remove grid
        legend.position = "none") +   # Remove legend
  xlim(0.69, 1.0) + ylim(0, 100) 

# Save the plot as a TIFF file
ggsave("Na_demage.tiff", plot = plot, device = "tiff", dpi = 300, width = 8, height = 8)




# Ensure necessary libraries are loaded
library(dplyr)
library(ggplot2)

# Add a new column to each data frame to indicate the source of the data
H_vs_H$Source <- "H_vs_H"
H_vs_Na$Source <- "H_vs_Na"

# Sort from high to low
H_vs_H_sorted <- H_vs_H %>%
  filter(Metabolite_Name_1 != Metabolite_Name_2) %>%
  arrange(desc(remddp))


# Find common indices between H_vs_H and H_vs_Na
common_indices <- intersect(H_vs_H_sorted$Index, H_vs_Na$Index)

# Filter both data frames to retain only rows with common indices
H_vs_H_filtered <- H_vs_H_sorted %>% filter(Index %in% common_indices)
H_vs_Na_filtered <- H_vs_Na %>% filter(Index %in% common_indices)

# Reorder H_vs_Na by the order of H_vs_H_filtered based on Index
H_vs_Na_filtered <- H_vs_Na_filtered %>%
  arrange(match(Index, H_vs_H_filtered$Index))

# Reassign the row index to both data frames (ensuring they are now aligned)
H_vs_H_filtered$New_Index <- seq_along(H_vs_H_filtered$Index)
H_vs_Na_filtered$New_Index <- seq_along(H_vs_Na_filtered$Index)

# Combine the two data frames into one, keeping the reassigned index
combined_data <- rbind(H_vs_H_filtered, H_vs_Na_filtered)

# Plot with reassigned Index on the x-axis and remddp on the y-axis
plot <- ggplot(combined_data, aes(x = New_Index, y = remddp, color = Source)) +
  geom_point(size = 1) +  # Transparent scatter points
  scale_shape_manual(values = c("H_vs_H" = 15, "Na_vs_H" = 15)) +  # Square for H_vs_H, Circle for Na_vs_H
  scale_color_manual(values = c("H_vs_H" = "#2667FF", "H_vs_Na" = "#E25D75")) +  # Custom colors
  theme_minimal() +  # Clean theme
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid = element_blank(),  # Remove grid
        legend.position = "none") +   # Remove legend
  xlab("Reassigned Index") +  # Label for x-axis
  ylab("Average dp") +  # Label for y-axis
  ylim(0, 100)  # Adjust y-axis limits for remddp

# Save the plot as a TIFF file
ggsave("Na_demage_sorted_ID.tiff", plot = plot, device = "tiff", width = 3.8, height = 3.8, dpi = 300)

# Calculate the count of red dots with remddp < 70
red_dots_under_70 <- combined_data %>%
  filter(Source == "H_vs_Na" & remddp < 70) %>%
  nrow()

# Total number of red dots
total_red_dots <- combined_data %>%
  filter(Source == "H_vs_Na") %>%
  nrow()

# Calculate the percentage
percentage_red_dots_under_70 <- (red_dots_under_70 / total_red_dots) * 100

# Print results
cat("Count of red dots with remddp < 70:", red_dots_under_70, "\n")
cat("Total red dots:", total_red_dots, "\n")
cat("Percentage of red dots with remddp < 70:", round(percentage_red_dots_under_70, 2), "%\n")







# Ensure necessary libraries are loaded
library(dplyr)
library(ggplot2)

# Add a new column to each data frame to indicate the source of the data
H_vs_H$Source <- "H_vs_H"
H_vs_Na$Source <- "H_vs_Na"

# Sort from high to low
H_vs_H_sorted <- H_vs_H %>%
  filter(Metabolite_Name_1 != Metabolite_Name_2) %>%
  arrange(desc(remddp))


# Find common indices between H_vs_H and H_vs_Na
common_indices <- intersect(H_vs_H_sorted$Index, H_vs_Na$Index)

# Filter both data frames to retain only rows with common indices
H_vs_H_filtered <- H_vs_H_sorted %>% filter(Index %in% common_indices)
H_vs_Na_filtered <- H_vs_Na %>% filter(Index %in% common_indices)

# Reorder H_vs_Na by the order of H_vs_H_filtered based on Index
H_vs_Na_filtered <- H_vs_Na_filtered %>%
  arrange(match(Index, H_vs_H_filtered$Index))

# Reassign the row index to both data frames (ensuring they are now aligned)
H_vs_H_filtered$New_Index <- seq_along(H_vs_H_filtered$Index)
H_vs_Na_filtered$New_Index <- seq_along(H_vs_Na_filtered$Index)

# Combine the two data frames into one, keeping the reassigned index
combined_data <- rbind(H_vs_H_filtered, H_vs_Na_filtered)

# Plot with reassigned Index on the x-axis and remddp on the y-axis
plot <- ggplot(combined_data, aes(x = New_Index, y = remddp, color = Source)) +
  geom_point(size = 1.4) +  # Transparent scatter points
  scale_color_manual(values = c("H_vs_H" = "#2667FF", "H_vs_Na" = "#E25D75")) +  # Custom colors
  theme_minimal() +  # Clean theme
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid = element_blank(),  # Remove grid
        legend.position = "none") +   # Remove legend
  xlab("Reassigned Index") +  # Label for x-axis
  ylab("Average dp") +  # Label for y-axis
  ylim(0, 100)  # Adjust y-axis limits for remddp

ggsave(
  filename = "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/Na_demage_sorted_ID.tiff",
  plot = plot,
  device = "tiff",
  width = 3.8,
  height = 3.8,
  dpi = 150
)


# Calculate the count of red dots with remddp < 70
red_dots_under_70 <- combined_data %>%
  filter(Source == "H_vs_Na" & remddp < 70) %>%
  nrow()

# Total number of red dots
total_red_dots <- combined_data %>%
  filter(Source == "H_vs_Na") %>%
  nrow()

# Calculate the percentage
percentage_red_dots_under_70 <- (red_dots_under_70 / total_red_dots) * 100

# Print results
cat("Count of red dots with remddp < 70:", red_dots_under_70, "\n")
cat("Total red dots:", total_red_dots, "\n")
cat("Percentage of red dots with remddp < 70:", round(percentage_red_dots_under_70, 2), "%\n")