### Section 1
# Load necessary libraries
library(dplyr)
library(data.table)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(foreach)
library(doParallel)

# Force single-core mode for compatibility on Windows
options(mc.cores = 1)
setDTthreads(1)

# List of adducts to process
adducts <- c("[M+NH4]+", "[M+Na]+", "[M+K]+", "[M+H-H2O]+", "[2M+H]+")

# Base directory path
base_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"

# Create output directories for each adduct if they don't exist
for (adduct in adducts) {
  adduct_folder <- gsub("[\\[\\]\\+\\-]", "", adduct)
  output_dir <- file.path(base_dir, adduct_folder)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
}

# Define CE difference intervals (-30 to 30, step = 3)
ce_intervals <- list()
for (i in seq(-30, 27, by = 3)) {
  ce_intervals[[length(ce_intervals) + 1]] <- c(i, i + 3)
}

# Outer loop: Process one adduct at a time
for (adduct in adducts) {
  print(paste("Processing adduct:", adduct))
  
  # Set working directory to the location of files
  setwd(base_dir)
  
  # Determine adduct folder name (removing special characters)
  adduct_folder <- gsub("[\\[\\]\\+\\-]", "", adduct)
  output_dir <- file.path(base_dir, adduct_folder)
  
  # Load reference adduct [M+H]+
  tryCatch({
    adduct_H <- read.xlsx("adduct_[M+H]+.xlsx") %>% as.data.table()
    
    # Load data for the current adduct
    adduct_filename <- paste0("adduct_", gsub("[\\[\\]\\+\\-]", "", adduct), ".xlsx")
    adduct_other <- read.xlsx(adduct_filename) %>% as.data.table()
    
    # Retain only rows where compound names exist in both data frames
    adduct_H <- adduct_H[WHOLE_COMPOUND_NAME %in% adduct_other$WHOLE_COMPOUND_NAME]
    adduct_other <- adduct_other[WHOLE_COMPOUND_NAME %in% adduct_H$WHOLE_COMPOUND_NAME]
    
    # Initialize cluster for parallel processing
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    
    # Inner loop: Process all CE intervals for the current adduct
    for (interval in ce_intervals) {
      print(paste("Processing CE interval:", interval[1], "to", interval[2]))
      
      # Initialize an empty list to store matched pairs
      pairs_list <- list()
      
      # Perform parallel processing using foreach
      pairs_list <- foreach(i = 1:nrow(adduct_H), .combine = rbind, .errorhandling = 'remove') %dopar% {
        temp_list <- list()
        for (j in 1:nrow(adduct_other)) {
          if (adduct_H$WHOLE_COMPOUND_NAME[i] == adduct_other$WHOLE_COMPOUND_NAME[j]) {
            ce_diff <- (adduct_H$WHOLE_COLLISION_ENERGY[i] - adduct_other$WHOLE_COLLISION_ENERGY[j])
            if (ce_diff >= interval[1] & ce_diff < interval[2]) {
              temp_list <- append(temp_list, list(data.frame(
                Metabolite_Name = adduct_H$WHOLE_COMPOUND_NAME[i],
                CE1 = adduct_H$WHOLE_COLLISION_ENERGY[i],
                CE2 = adduct_other$WHOLE_COLLISION_ENERGY[j],
                Adduct1 = "[M+H]+",
                Adduct2 = adduct,
                Index_1 = adduct_H$Index[i],
                Index_2 = adduct_other$Index[j],
                precursor1 = adduct_H$WHOLE_PRECURSOR_mz[i],
                precursor2 = adduct_other$WHOLE_PRECURSOR_mz[j],
                whole_mz1 = adduct_H$WHOLE_mz[i],
                whole_mz2 = adduct_other$WHOLE_mz[j],
                whole_int1 = adduct_H$WHOLE_Int[i],
                whole_int2 = adduct_other$WHOLE_Int[j]
              )))
            }
          }
        }
        do.call(rbind, temp_list)
      }
      
      # Combine and save results
      if (nrow(pairs_list) > 0) {
        filename <- paste0(adduct, "_pairs_", interval[1], "_to_", interval[2], ".csv")
        write.csv(pairs_list, file.path(output_dir, filename), row.names = FALSE)
        print(paste("Saved", nrow(pairs_list), "pairs to", file.path(output_dir, filename)))
      } else {
        print(paste("No pairs found for adduct", adduct, "in interval", interval[1], "to", interval[2]))
      }
    }
    
    # Stop cluster
    stopCluster(cl)
    
  }, error = function(e) {
    print(paste("Error processing adduct", adduct, ":", e$message))
  })
}




### Section 1
# Load necessary libraries
library(dplyr)
library(data.table)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(foreach)
library(doParallel)

# Force single-core mode for compatibility on Windows
options(mc.cores = 1)
setDTthreads(1)

# List of adducts to process
adducts <- c("[M+Cl]-", "[2M-H]-", "[M-H-H2O]-")

# Base directory path
base_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"

# Create output directories for each adduct if they don't exist
for (adduct in adducts) {
  adduct_folder <- gsub("[\\[\\]\\+\\-]", "", adduct)
  output_dir <- file.path(base_dir, adduct_folder)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
}

# Define CE difference intervals (-30 to 30, step = 3)
ce_intervals <- list()
for (i in seq(-30, 27, by = 3)) {
  ce_intervals[[length(ce_intervals) + 1]] <- c(i, i + 3)
}

# Outer loop: Process one adduct at a time
for (adduct in adducts) {
  print(paste("Processing adduct:", adduct))
  
  # Set working directory to the location of files
  setwd(base_dir)
  
  # Determine adduct folder name (removing special characters)
  adduct_folder <- gsub("[\\[\\]\\+\\-]", "", adduct)
  output_dir <- file.path(base_dir, adduct_folder)
  
  # Load reference adduct [M+H]+
  tryCatch({
    adduct_H <- read.xlsx("adduct_[M-H]-.xlsx") %>% as.data.table()
    
    # Load data for the current adduct
    adduct_filename <- paste0("adduct_", gsub("[\\[\\]\\+\\-]", "", adduct), ".xlsx")
    adduct_other <- read.xlsx(adduct_filename) %>% as.data.table()
    
    # Retain only rows where compound names exist in both data frames
    adduct_H <- adduct_H[WHOLE_COMPOUND_NAME %in% adduct_other$WHOLE_COMPOUND_NAME]
    adduct_other <- adduct_other[WHOLE_COMPOUND_NAME %in% adduct_H$WHOLE_COMPOUND_NAME]
    
    # Initialize cluster for parallel processing
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    
    # Inner loop: Process all CE intervals for the current adduct
    for (interval in ce_intervals) {
      print(paste("Processing CE interval:", interval[1], "to", interval[2]))
      
      # Initialize an empty list to store matched pairs
      pairs_list <- list()
      
      # Perform parallel processing using foreach
      pairs_list <- foreach(i = 1:nrow(adduct_H), .combine = rbind, .errorhandling = 'remove') %dopar% {
        temp_list <- list()
        for (j in 1:nrow(adduct_other)) {
          if (adduct_H$WHOLE_COMPOUND_NAME[i] == adduct_other$WHOLE_COMPOUND_NAME[j]) {
            ce_diff <- (adduct_H$WHOLE_COLLISION_ENERGY[i] - adduct_other$WHOLE_COLLISION_ENERGY[j])
            if (ce_diff >= interval[1] & ce_diff < interval[2]) {
              temp_list <- append(temp_list, list(data.frame(
                Metabolite_Name = adduct_H$WHOLE_COMPOUND_NAME[i],
                CE1 = adduct_H$WHOLE_COLLISION_ENERGY[i],
                CE2 = adduct_other$WHOLE_COLLISION_ENERGY[j],
                Adduct1 = "[M-H]-",
                Adduct2 = adduct,
                Index_1 = adduct_H$Index[i],
                Index_2 = adduct_other$Index[j],
                INSTRUMENT_1 = adduct_H$INSTRUMENT[i],
                INSTRUMENT_2 = adduct_other$Index[j],
                precursor1 = adduct_H$WHOLE_PRECURSOR_mz[i],
                precursor2 = adduct_other$WHOLE_PRECURSOR_mz[j],
                whole_mz1 = adduct_H$WHOLE_mz[i],
                whole_mz2 = adduct_other$WHOLE_mz[j],
                whole_int1 = adduct_H$WHOLE_Int[i],
                whole_int2 = adduct_other$WHOLE_Int[j]
              )))
            }
          }
        }
        do.call(rbind, temp_list)
      }
      
      # Combine and save results
      if (nrow(pairs_list) > 0) {
        filename <- paste0(adduct, "_pairs_", interval[1], "_to_", interval[2], ".csv")
        write.csv(pairs_list, file.path(output_dir, filename), row.names = FALSE)
        print(paste("Saved", nrow(pairs_list), "pairs to", file.path(output_dir, filename)))
      } else {
        print(paste("No pairs found for adduct", adduct, "in interval", interval[1], "to", interval[2]))
      }
    }
    
    # Stop cluster
    stopCluster(cl)
    
  }, error = function(e) {
    print(paste("Error processing adduct", adduct, ":", e$message))
  })
}















# Load required libraries
library(parallel)
library(doParallel)
library(dplyr)

# Core spectral matching function
process_spectral_matching <- function(pairs_df, noise_threshold = 5, mz_tol = 0.02) {
  # Detect the number of cores and set up cluster
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
  
  # Export parameters and functions to the cluster
  clusterExport(cl, c('noise_threshold', 'mz_tol'), envir = environment())
  
  # Completely define the row processing function within the main function
  process_row <- function(i) {
    if (i %% 100 == 0) cat(paste("Processing row", i, "\n"))
    
    # Extract precursor m/z
    precursor1 <- pairs_df$precursor1[i]
    precursor2 <- pairs_df$precursor2[i]
    
    # Extract and clean fragment m/z and intensities
    whole_mz_cleaned_1 <- as.numeric(unlist(strsplit(as.character(pairs_df$whole_mz1[i]), ",")))
    whole_int_cleaned_1 <- as.numeric(unlist(strsplit(as.character(pairs_df$whole_int1[i]), ",")))
    whole_mz_cleaned_2 <- as.numeric(unlist(strsplit(as.character(pairs_df$whole_mz2[i]), ",")))
    whole_int_cleaned_2 <- as.numeric(unlist(strsplit(as.character(pairs_df$whole_int2[i]), ",")))
    
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
    remddp <- (result$modified_cosine_similarity + FP) / 2
    
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
      remddp = round(remddp)
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
  results <- parLapply(cl, 1:nrow(pairs_df), process_row)
  
  # Stop the cluster
  stopCluster(cl)
  
  end_time <- Sys.time()
  cat(paste("Processing completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n"))
  
  # Add results to the dataframe, replacing NA with 0
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
      remddp = sapply(results, `[[`, "remddp")
    ) %>%
    # Replace NA with 0 for all columns
    mutate(across(c(dp, remddp, dp_aligned_fragments, remddp_aligned_fragments, 
                    Frag_Number_MH, Frag_Number_Na, total_aligned, 
                    Direct_match, Shifted_match, total_fragments, 
                    FP, ts_dp, remddp), ~replace_na(., 0)))
  
  return(pairs_df)
}

# Function to process files in a directory with NA to 0 replacement
prepare_and_process <- function(file_path, noise_threshold = 5, mz_tol = 0.02) {
  # Read the CSV file
  pairs_df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Process the spectral matching
  processed_df <- process_spectral_matching(pairs_df, noise_threshold, mz_tol)
  
  # Replace any remaining NA values with 0 before saving
  processed_df <- processed_df %>%
    mutate(across(everything(), ~replace_na(., 0)))
  
  # Optional: Save the processed file with NA replaced by 0
  output_file <- sub("\\.csv$", "_processed.csv", file_path)
  write.csv(processed_df, output_file, row.names = FALSE)
  
  return(processed_df)
}

# Remaining functions stay the same as in the original script
process_adduct_files <- function(base_dir, adduct, noise_threshold = 5, mz_tol = 0.02) {
  # Determine adduct folder name (removing special characters)
  adduct_folder <- gsub("[\\[\\]\\+\\-]", "", adduct)
  input_dir <- file.path(base_dir, adduct_folder)
  
  # Check if input directory exists
  if (!dir.exists(input_dir)) {
    cat(paste("Directory not found:", input_dir, "\n"))
    return(NULL)
  }
  
  # Get list of CSV files in the directory
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    cat(paste("No CSV files found in", input_dir, "\n"))
    return(NULL)
  }
  
  # Process each file
  processed_files <- lapply(csv_files, function(file_path) {
    cat(paste("Processing file:", file_path, "\n"))
    tryCatch({
      processed_df <- prepare_and_process(file_path, noise_threshold, mz_tol)
      cat(paste("Processed file:", file_path, "successfully\n"))
      return(processed_df)
    }, error = function(e) {
      cat(paste("Error processing file", file_path, ":", e$message, "\n"))
      return(NULL)
    })
  })
  
  # Filter out any NULL results
  processed_files <- processed_files[!sapply(processed_files, is.null)]
  
  return(list(
    adduct = adduct,
    folder = adduct_folder,
    total_files = length(csv_files),
    processed_files = length(processed_files)
  ))
}

# Main processing function remains the same
process_all_adducts <- function(base_dir, adducts, noise_threshold = 5, mz_tol = 0.02) {
  # Process each adduct
  results <- lapply(adducts, function(adduct) {
    cat(paste("\n==== Processing adduct:", adduct, "====\n"))
    process_adduct_files(base_dir, adduct, noise_threshold, mz_tol)
  })
  
  # Generate summary
  summary <- do.call(rbind, results)
  print(summary)
  
  return(results)
}

# Set up processing parameters
adducts_negative <- c("[M+Cl]-", "[2M-H]-", "[M-H-H2O]-")
adducts_positive <- c("[M+NH4]+", "[M+Na]+", "[M+K]+", "[M+H-H2O]+", "[2M+H]+")
base_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"

# Processing parameters
noise_threshold <- 5
mz_tol <- 0.02

# Run the analysis
tryCatch({
  start_time <- Sys.time()
  
  # Process negative mode adducts
  results_negative <- process_all_adducts(
    base_dir, 
    adducts_negative, 
    noise_threshold = noise_threshold, 
    mz_tol = mz_tol
  )
  
  # Process positive mode adducts
  results_positive <- process_all_adducts(
    base_dir, 
    adducts_positive, 
    noise_threshold = noise_threshold, 
    mz_tol = mz_tol
  )
  
  end_time <- Sys.time()
  
  # Summary
  cat("Processing completed successfully.\n")
  cat(sprintf("Total processing time: %s minutes\n", 
              round(difftime(end_time, start_time, units = "mins"), 2)))
  
}, error = function(e) {
  cat("An error occurred during processing:\n")
  cat(e$message, "\n")
})














# CE difference intervals to explore
# Define CE difference intervals (-30 to 30, step = 3)
ce_intervals <- list()
for (i in seq(-30, 27, by = 3)) {
  ce_intervals[[length(ce_intervals) + 1]] <- c(i, i + 3)
}

# Output directory for all plots
output_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Load required libraries
library(ggplot2)
library(reshape2)
library(tidyr)

# Define adduct information (adduct name and corresponding directory)
adduct_directories <- list(
  c("[M+NH4]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/[M+NH4]+"),
  c("[M+H-H2O]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/[M+H-H2O]+"),
  c("[2M+H]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/[2M+H]+"),
  c("[M+Na]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/[M+Na]+"),
  c("[M+K]+", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/[M+K]+"),
  c("[M+Cl]-", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/[M+Cl]-"),
  c("[M-H-H2O]-", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/[M-H-H2O]-"),
  c("[2M-H]-", "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20/[2M-H]-")
)

# Initialize a dataframe to store summary statistics for all adducts
all_adducts_summary <- data.frame(
  Adduct = character(),
  CE_Interval = character(),
  Interval_Midpoint = numeric(),
  Avg_dp = numeric(),
  Avg_FP = numeric(),
  Avg_remddp = numeric(),
  Avg_ts_remddp = numeric(),
  Count = numeric(),
  stringsAsFactors = FALSE
)

# Create a function to generate density data
create_density_df <- function(data_column) {
  if (length(data_column) == 0 || all(is.na(data_column))) {
    return(NULL)
  }
  density_data <- density(data_column, na.rm = TRUE)
  density_df <- data.frame(x = density_data$x, y = density_data$y)
  return(density_df)
}

# Process each adduct
for (adduct_info in adduct_directories) {
  adducts <- adduct_info[1]
  directory <- adduct_info[2]
  
  cat("\nProcessing adduct:", adducts, "from directory:", directory, "\n")
  
  # Set working directory for saving
  setwd(directory)
  
  # Generate input files list dynamically
  input_files <- paste0(adducts, "_pairs_", 
                        seq(-30, 27, by = 3), "_to_", 
                        seq(-27, 30, by = 3), "_processed.csv")
  
  # Explicitly define CE intervals based on input files
  defined_intervals <- lapply(seq_along(input_files), function(i) {
    c(seq(-30, 27, by = 3)[i], seq(-27, 30, by = 3)[i])
  })
  
  # Initialize an empty data frame to store all combined density data
  combined_density_df <- data.frame()
  
  # Initialize an empty data frame to store summary statistics
  summary_df <- data.frame(
    CE_Interval = character(),
    Avg_dp = numeric(),
    Avg_FP = numeric(),
    Avg_remddp = numeric(),
    Avg_ts_remddp = numeric(),
    Count = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each input file
  for (i in seq_along(input_files)) {
    filename <- input_files[i]
    interval <- defined_intervals[[i]]
    
    # Check if the file exists before reading
    if (!file.exists(filename)) {
      cat("File not found:", filename, "\n")
      next
    }
    
    # Read the filtered_cosine_similarity_df
    filtered_cosine_similarity_df <- read.csv(filename, stringsAsFactors = FALSE)
    
    # Ensure numeric conversion with error handling
    filtered_cosine_similarity_df$dp <- as.numeric(filtered_cosine_similarity_df$dp)
    filtered_cosine_similarity_df$FP <- as.numeric(filtered_cosine_similarity_df$FP)
    filtered_cosine_similarity_df$remddp <- as.numeric(filtered_cosine_similarity_df$remddp)
    filtered_cosine_similarity_df$ts_dp <- as.numeric(filtered_cosine_similarity_df$ts_dp)
    
    # Filter for identical instrument type pairs
    filtered_cosine_similarity_df <- filtered_cosine_similarity_df %>%
      filter(INSTRUMENT_TYPE_1 == INSTRUMENT_TYPE_2)

    # Create density data for each variable
    density_dp <- create_density_df(filtered_cosine_similarity_df$dp)
    if (!is.null(density_dp)) density_dp$Variable <- "dp"
    
    density_FP <- create_density_df(filtered_cosine_similarity_df$FP)
    if (!is.null(density_FP)) density_FP$Variable <- "FP"
    
    density_remddp <- create_density_df(filtered_cosine_similarity_df$remddp)
    if (!is.null(density_remddp)) density_remddp$Variable <- "remddp"
    
    density_ts_remddp <- create_density_df(filtered_cosine_similarity_df$ts_dp)
    if (!is.null(density_ts_remddp)) density_ts_remddp$Variable <- "ts_dp"
    
    # Combine density data for this interval
    interval_density_df <- rbind(
      if (!is.null(density_dp)) density_dp,
      if (!is.null(density_FP)) density_FP,
      if (!is.null(density_remddp)) density_remddp,
      if (!is.null(density_ts_remddp)) density_ts_remddp
    )
    
    # Add interval information for grouping
    if (!is.null(interval_density_df)) {
      interval_density_df$CE_Interval <- paste(interval[1], "to", interval[2])
      
      # Append to the combined density data frame
      combined_density_df <- rbind(combined_density_df, interval_density_df)
    }
    
    # Calculate averages for this interval
    avg_dp <- mean(filtered_cosine_similarity_df$dp, na.rm = TRUE)
    avg_FP <- mean(filtered_cosine_similarity_df$FP, na.rm = TRUE)
    avg_remddp <- mean(filtered_cosine_similarity_df$remddp, na.rm = TRUE)
    avg_ts_remddp <- mean(filtered_cosine_similarity_df$ts_remddp, na.rm = TRUE)
    
    # Count the number of cases in this interval
    count <- nrow(filtered_cosine_similarity_df)
    
    # Append the calculated averages and count to the summary data frame
    summary_df <- rbind(summary_df, data.frame(
      CE_Interval = paste(interval[1], "to", interval[2]),
      Avg_dp = avg_dp,
      Avg_FP = avg_FP,
      Avg_remddp = avg_remddp,
      Avg_ts_remddp = avg_ts_remddp,
      Count = count,
      stringsAsFactors = FALSE
    ))
  }
  
  # Skip if no data was found for this adduct
  if (nrow(combined_density_df) == 0) {
    cat("No data found for adduct:", adducts, "\n")
    next
  }
  
  # Add a column for the midpoint of the CE intervals
  summary_df$Interval_Midpoint <- sapply(strsplit(summary_df$CE_Interval, " to "), function(x) {
    mean(as.numeric(x))
  })
  
  summary_df$Interval_Midpoint <- 0 - summary_df$Interval_Midpoint
  
  # Add adduct information to the summary dataframe
  summary_df$Adduct <- adducts
  
  # Append this adduct's summary to the all_adducts_summary dataframe
  all_adducts_summary <- rbind(all_adducts_summary, summary_df)
  
  # Reshape the summary_df to a long format for easier plotting
  long_summary_df <- melt(
    summary_df,
    id.vars = "Interval_Midpoint",
    measure.vars = c("Avg_remddp"),
    variable.name = "Metric",
    value.name = "Value"
  )
  
  # Find the coordinates of the maximum point
  max_point <- long_summary_df[which.max(long_summary_df$Value), ]
  
  # Modify the plotting section to handle color consistently
  summary_plot <- ggplot(long_summary_df, aes(x = Interval_Midpoint, y = Value)) +
    geom_line(color = "black", size = 0.8) +   # Black line
    geom_point(
      data = long_summary_df,
      color = "black", 
      size = 0.8
    ) +  # Black points
    geom_point(
      data = max_point,
      color = "red", 
      size = 1.2
    ) +  # Red point for max value
    xlim(-35, 35) +  # Set x-axis limits
    ylim(0, 100) +  # Set y-axis limits
    labs(
      title = "",  # Remove title
      x = "",      # Remove x-axis label
      y = ""       # Remove y-axis label
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",       # Remove legend
      plot.title = element_blank(),   # Remove title
      axis.title = element_blank(),   # Remove axis titles
      axis.text = element_blank(),    # Remove axis labels
      axis.ticks = element_blank(),   # Remove axis ticks
      panel.grid = element_blank(),   # Remove grid
      panel.background = element_rect(fill = "transparent", colour = NA),  # Transparent background
      plot.background = element_rect(fill = "transparent", colour = NA),   # Transparent background
      panel.border = element_rect(color = "black", fill = NA, size = 0.3)  # Add full frame (right, left, top, bottom)
    )
  
  print(summary_plot)
  
  # Save summary plot as TIFF with transparent background
  ggsave(
    filename = file.path(output_dir, paste0(adducts, "_Summary_Statistics.tiff")),
    plot = summary_plot,
    device = "tiff",
    width = 1.5,
    height = 1.5,
    dpi = 150,
    bg = "transparent"  # Set background to transparent
  )
  
  cat("Completed processing for adduct:", adducts, "\n")
}

# Save the all_adducts_summary dataframe to a CSV file
write.csv(all_adducts_summary, file.path(output_dir, "all_adducts_summary.csv"), row.names = FALSE)

cat("All adducts processed successfully\n")
