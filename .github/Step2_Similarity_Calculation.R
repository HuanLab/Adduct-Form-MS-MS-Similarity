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


# Wrapper function for reversed modified dot product calculation
calculate_reversed_modified_dot_product <- function(ref_fragMz, ref_fragInt, exp_fragMz, exp_fragInt, precursor_mz_diff) {
  result <- calculate_dot_product_and_modified_cosine(
    ref_fragMz, ref_fragInt, exp_fragMz, exp_fragInt, precursor_mz_diff
  )
  
  return(list(
    dp = result$dot_product,
    remddp = result$modified_cosine_similarity,
    original_aligned_fragments = result$original_aligned_fragments,
    remddp_aligned_fragments = result$modified_aligned_fragments
  ))
}

#------------------------------------------------------------
# Single Row Processing
#------------------------------------------------------------

# Function to process a single row
process_single_row <- function(i, pairs_df, file_name, noise_threshold) {
  if (i %% 100 == 0) print(paste("Processing row", i, "of", nrow(pairs_df), "in", file_name))
  
  # Extract data for the row
  precursor1 <- pairs_df$precursor1[i]
  precursor2 <- pairs_df$precursor2[i]
  
  # Split strings once and convert to numeric
  whole_mz_cleaned_1 <- as.numeric(unlist(strsplit(as.character(pairs_df$whole_mz1[i]), ",")))
  whole_int_cleaned_1 <- as.numeric(unlist(strsplit(as.character(pairs_df$whole_int1[i]), ",")))
  whole_mz_cleaned_2 <- as.numeric(unlist(strsplit(as.character(pairs_df$whole_mz2[i]), ",")))
  whole_int_cleaned_2 <- as.numeric(unlist(strsplit(as.character(pairs_df$whole_int2[i]), ",")))
  
  # Filter by noise threshold
  keep_indices_1 <- whole_int_cleaned_1 >= noise_threshold
  keep_indices_2 <- whole_int_cleaned_2 >= noise_threshold
  
  whole_mz_cleaned_1 <- whole_mz_cleaned_1[keep_indices_1]
  whole_int_cleaned_1 <- whole_int_cleaned_1[keep_indices_1]
  whole_mz_cleaned_2 <- whole_mz_cleaned_2[keep_indices_2]
  whole_int_cleaned_2 <- whole_int_cleaned_2[keep_indices_2]
  
  # Check if we have valid data after filtering
  if (length(whole_mz_cleaned_1) == 0 || length(whole_mz_cleaned_2) == 0) {
    return(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
  }
  
  lower_precursor_mz <- min(precursor1, precursor2)
  
  # Filter for dot product calculation
  mz_threshold_dot <- lower_precursor_mz - 0.02
  ref_filter_dot <- whole_mz_cleaned_2 <= mz_threshold_dot
  exp_filter_dot <- whole_mz_cleaned_1 <= mz_threshold_dot
  
  ref_fragMz_dot <- whole_mz_cleaned_2[ref_filter_dot]
  ref_fragInt_dot <- whole_int_cleaned_2[ref_filter_dot]  # Fixed: was using whole_mz_cleaned_2 incorrectly
  exp_fragMz_dot <- whole_mz_cleaned_1[exp_filter_dot]
  exp_fragInt_dot <- whole_int_cleaned_1[exp_filter_dot]  # Fixed: was using exp_filter_dot incorrectly
  
  # Filter for alignment calculation
  mz_threshold_align <- lower_precursor_mz + 0.02
  ref_filter_align <- whole_mz_cleaned_2 <= mz_threshold_align
  exp_filter_align <- whole_mz_cleaned_1 <= mz_threshold_align
  
  Frag_Number_MH <- sum(exp_filter_align)
  Frag_Number_Na <- sum(ref_filter_align)
  
  # Check if we have valid fragments for dot product
  if (length(ref_fragMz_dot) == 0 || length(exp_fragMz_dot) == 0) {
    return(c(0, 0, 0, 0, Frag_Number_MH, Frag_Number_Na, 0, 0, 0,
             Frag_Number_MH + Frag_Number_Na, 0, 0))
  }
  
  # Calculate dot products
  precursor_mz_diff <- precursor2 - precursor1
  
  reverse_dot_product_dp_result <- calculate_reversed_modified_dot_product(
    ref_fragMz_dot, ref_fragInt_dot, exp_fragMz_dot, exp_fragInt_dot, precursor_mz_diff
  )
  
  dp <- reverse_dot_product_dp_result$dp
  remddp <- reverse_dot_product_dp_result$remddp
  dp_aligned_fragments <- reverse_dot_product_dp_result$original_aligned_fragments
  remddp_aligned_fragments <- reverse_dot_product_dp_result$remddp_aligned_fragments
  
  # Calculate derived metrics
  total_aligned <- dp_aligned_fragments + remddp_aligned_fragments
  Direct_match <- ifelse(total_aligned == 0, 0, 100 * dp_aligned_fragments / total_aligned)
  Shifted_match <- ifelse(total_aligned == 0, 0, 100 * remddp_aligned_fragments / total_aligned)
  total_fragments <- Frag_Number_MH + Frag_Number_Na
  FP <- ifelse(total_fragments == 0, 0, 100 * 2 * total_aligned / total_fragments)
  
  ts_dp <- (dp + FP) / 2
  ts_remddp <- (remddp + FP) / 2
  
  return(c(
    round(dp), round(remddp), round(dp_aligned_fragments), round(remddp_aligned_fragments),
    Frag_Number_MH, Frag_Number_Na, round(total_aligned), round(Direct_match), round(Shifted_match), 
    round(total_fragments), round(FP), round(ts_dp), round(ts_remddp)
  ))
}

#------------------------------------------------------------
# Main Processing Function
#------------------------------------------------------------

# Main file processing function (Windows-compatible)
process_adduct_file <- function(file_path, output_dir, noise_threshold = 5) {
  file_name <- basename(file_path)
  print(paste("Reading file:", file_name))
  pairs_df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  
  
  # Choose appropriate parallel processing method based on platform
  if (requireNamespace("parallel", quietly = TRUE)) {
    num_cores <- max(1, parallel::detectCores() - 1)
    print(paste("Using", num_cores, "cores for processing"))
    
    if (.Platform$OS.type == "windows") {
      # For Windows, use parLapply with socket cluster
      cl <- parallel::makeCluster(num_cores)
      
      # Export required functions and data to the cluster
      parallel::clusterExport(cl, c("calculate_dot_product_and_modified_cosine", 
                                    "calculate_reversed_modified_dot_product",
                                    "process_single_row"))
      
      # Also need to export the full dataframe for Windows clusters
      parallel::clusterExport(cl, "pairs_df", envir = environment())
      
      results <- parallel::parLapply(cl, 1:nrow(pairs_df), function(i) {
        process_single_row(i, pairs_df, file_name, noise_threshold)
      })
      
      parallel::stopCluster(cl)
    } else {
      # For Unix/Linux, use mclapply
      results <- parallel::mclapply(1:nrow(pairs_df), function(i) {
        process_single_row(i, pairs_df, file_name, noise_threshold)
      }, mc.cores = num_cores)
    }
  } else {
    # Regular sequential processing
    print("Parallel package not available, using sequential processing")
    results <- lapply(1:nrow(pairs_df), function(i) {
      process_single_row(i, pairs_df, file_name, noise_threshold)
    })
  }
  
  # Combine results
  print("Combining results")
  result_matrix <- do.call(rbind, results)
  colnames(result_matrix) <- c("dp", "remddp", "dp_aligned_fragments", "remddp_aligned_fragments",
                               "Frag_Number_MH", "Frag_Number_Na", "total_aligned", "Direct_match", "Shifted_match",
                               "total_fragments", "FP", "ts_dp", "ts_remddp")
  
  # Update dataframe with results
  pairs_df[, colnames(result_matrix)] <- result_matrix
  
  # Replace NA values with 0
  pairs_df[is.na(pairs_df)] <- 0
  
  # Write output file
  output_file <- file.path(output_dir, gsub("\\.csv$", "_noise5.csv", file_name))
  print(paste("Writing results to:", output_file))
  write.csv(pairs_df, output_file, row.names = FALSE)
  return(output_file)
}

#------------------------------------------------------------
# Main Script 
#------------------------------------------------------------

# Set working directories
base_dir_pos <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"
base_dir_neg <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"
output_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Process all positive adduct files
pos_adduct_files <- c(
  file.path(base_dir_pos, "[M+Na]+_pairs_df.csv"),
  file.path(base_dir_pos, "[M+K]+_pairs_df.csv"),
  file.path(base_dir_pos, "[2M+H]+_pairs_df.csv"),
  file.path(base_dir_pos, "[M+NH4]+_pairs_df.csv"),
  file.path(base_dir_pos, "[M+H-H2O]+_pairs_df.csv")
)

# Process all negative adduct files
neg_adduct_files <- c(
  file.path(base_dir_neg, "[M+Cl]-_pairs_df.csv"),
  file.path(base_dir_neg, "[M-H-H2O]-_pairs_df.csv"),
  file.path(base_dir_neg, "[2M-H]-_pairs_df.csv")
)

# Process all files
all_files <- c(pos_adduct_files, neg_adduct_files)
results <- lapply(all_files, function(file) {
  tryCatch({
    print(paste("===== Processing file:", basename(file), "====="))
    process_adduct_file(file, output_dir, noise_threshold = 5)
  }, error = function(e) {
    message("Error processing file ", file, ": ", e$message)
    return(NULL)
  })
})

# Print summary of processed files
successful_files <- results[!sapply(results, is.null)]
cat("\n===== PROCESSING SUMMARY =====\n")
cat("Successfully processed", length(successful_files), "out of", length(all_files), "files.\n")
cat("New files with noise threshold of 5 have been saved to:\n")
cat(paste("  ", output_dir), "\n\n")
cat("Files processed:\n")
for (i in seq_along(successful_files)) {
  cat(paste("  ", i, ". ", basename(successful_files[[i]])), "\n")
}
