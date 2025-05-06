# Specify the directory path
directory_path <- "X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20"

# List all CSV files in the directory
csv_files <- list.files(path = directory_path, pattern = "\\.csv$", full.names = TRUE)

# Filter files with '3eV_md_dp' in their names
filtered_files <- csv_files[grepl("3eV_md_dp", csv_files)]

# Read the filtered CSV files into a list
csv_data <- lapply(filtered_files, read.csv)

# Optionally, name the list with the file names (excluding the directory path)
names(csv_data) <- basename(filtered_files)

# Print the list of data frames (if needed)
print(names(csv_data))



# Extract the Metabolite_Name column from each file into a list
metabolite_names <- lapply(csv_data, function(df) df[["Metabolite_Name"]])

# Find the unique Metabolite_Names in each file
unique_metabolite_names <- lapply(metabolite_names, unique)

# Combine all unique Metabolite_Names into a single vector
all_metabolites <- unlist(unique_metabolite_names)

# Count the occurrences of each Metabolite_Name across files
metabolite_counts <- table(all_metabolites)

# Find the metabolite(s) that exist in most of the files
most_common_metabolites <- names(metabolite_counts[metabolite_counts == max(metabolite_counts)])

# Print the result
print(most_common_metabolites)

# Target metabolite name
target_metabolite <- "3-[(4-Fluorophenyl)sulfonyl]propanoic acid"

# Get the full set of column names from all files
all_columns <- unique(unlist(lapply(csv_data, names)))

# Initialize an empty list to store filtered data
filtered_data <- list()

# Loop through each file and filter rows with the target Metabolite_Name
for (file_name in names(csv_data)) {
  df <- csv_data[[file_name]]
  
  # Add missing columns to the current DataFrame (if any)
  for (col in setdiff(all_columns, names(df))) {
    df[[col]] <- NA
  }
  
  # Reorder columns to match the full set of columns
  df <- df[all_columns]
  
  # Filter rows for the target metabolite
  filtered_rows <- subset(df, Metabolite_Name == target_metabolite)
  
  # If there are matching rows, add them to the list
  if (nrow(filtered_rows) > 0) {
    filtered_rows$file_name <- file_name # Add file name as a column
    filtered_data[[file_name]] <- filtered_rows
  }
}

# Combine all filtered data into a single DataFrame
combined_data <- do.call(rbind, filtered_data)

# Print the combined DataFrame
print(head(combined_data))

combined_data$Average_CE <- (combined_data$CE1 + combined_data$CE2) / 2

print(head(combined_data))

# Group by Adduct2 and filter rows with the closest Average_CE to 15
filtered_data <- combined_data %>%
  group_by(Adduct2) %>%
  filter(abs(Average_CE - 15) == min(abs(Average_CE - 15))) %>%
  ungroup()

# View the filtered data
print(filtered_data)

filtered_data$Metabolites_CE <- paste0(filtered_data$Metabolite_Name, "_", filtered_data$Average_CE)

# Set metabolite name
metabolite_name <- "3-[(4-Fluorophenyl)sulfonyl]propanoic acid_15.18"

# Manually set the row index 'i'
i <- 2  # You can manually change this value to select the row
ms2_tol <- 0.02

# Extract and process whole_mz_cleaned_1 from the i-th matching row
whole_mz_cleaned_1_str <- filtered_data %>%
  filter(Metabolites_CE == metabolite_name) %>%
  slice(i) %>%  # Manually set 'i'
  pull(whole_mz_cleaned_1)
whole_mz_cleaned_1 <- as.numeric(strsplit(whole_mz_cleaned_1_str, ",")[[1]])

# Extract and process whole_int_cleaned_1 from the i-th matching row
whole_int_cleaned_1_str <- filtered_data %>%
  filter(Metabolites_CE == metabolite_name) %>%
  slice(i) %>%
  pull(whole_int_cleaned_1)
whole_int_cleaned_1 <- as.numeric(strsplit(whole_int_cleaned_1_str, ",")[[1]])

# Extract and process whole_mz_cleaned_2 from the i-th matching row
whole_mz_cleaned_2_str <- filtered_data %>%
  filter(Metabolites_CE == metabolite_name) %>%
  slice(i) %>%
  pull(whole_mz_cleaned_2)
whole_mz_cleaned_2 <- as.numeric(strsplit(whole_mz_cleaned_2_str, ",")[[1]])

# Extract and process whole_int_cleaned_2 from the i-th matching row
whole_int_cleaned_2_str <- filtered_data %>%
  filter(Metabolites_CE == metabolite_name) %>%
  slice(i) %>%
  pull(whole_int_cleaned_2)
whole_int_cleaned_2 <- as.numeric(strsplit(whole_int_cleaned_2_str, ",")[[1]])

# Normalize intensities
exp_fragInt <- whole_int_cleaned_2
ref_fragMzInt <- whole_int_cleaned_1

# Round normalized intensities to 2 decimals
exp_fragInt <- round(exp_fragInt, 2)
ref_fragMzInt <- round(ref_fragMzInt, 2)

# Create data frames for experimental and reference spectra
exp_data <- data.frame(mz = round(whole_mz_cleaned_2, 2), intensity = exp_fragInt, type = "Experimental")
ref_data <- data.frame(mz = round(whole_mz_cleaned_1, 2), intensity = ref_fragMzInt, type = "Reference")

# Extract precursor m/z values
precursor1 <- filtered_data %>%
  filter(Metabolites_CE == metabolite_name) %>%
  slice(i) %>%
  pull(precursor1)

precursor2 <- filtered_data %>%
  filter(Metabolites_CE == metabolite_name) %>%
  slice(i) %>%
  pull(precursor2)

# Create all_data data frame for plotting
all_data <- rbind(
  data.frame(mz = exp_data$mz, intensity = exp_data$intensity, type = "Experimental"),
  data.frame(mz = ref_data$mz, intensity = -ref_data$intensity, type = "Reference")
)


# Plot original data with custom colors
plot(all_data$mz, all_data$intensity, type = 'h', xlim = c(50, 220), 
     col = ifelse(all_data$type == "Experimental", "blue", "red"), 
     xlab = "m/z", ylab = "Relative Intensity", lwd = 10)

# Add black line at y = 0
abline(h = 0, col = "black", lwd = 8)

# Add text annotations for [M+NH4]+ and [M+H]+
text(x = min(all_data$mz) + 215, y = max(all_data$intensity) - 100, labels = "MS/MS of [M+Na]+", pos = 4, col = "blue", cex = 1)
text(x = min(all_data$mz) + 215, y = min(all_data$intensity) + 100, labels = "MS/MS of [M+H]+", pos = 4, col = "red", cex = 1)

# Function to find the closest match within tolerance
find_closest <- function(mz_vector, mz_value, tolerance) {
  abs_diff <- abs(mz_vector - mz_value)
  min_diff <- min(abs_diff)
  if (min_diff <= tolerance) {
    return(mz_vector[which.min(abs_diff)])
  } else {
    return(NA)
  }
}

# Find the closest m/z values for the precursor ions in each spectrum
closest_mz_MH <- find_closest(all_data$mz, precursor1, ms2_tol)
closest_mz_Na <- find_closest(all_data$mz, precursor2, ms2_tol)

# Add a red upward-pointing triangle for the precursor in the [M+H]+ spectrum
if (!is.na(closest_mz_MH)) {
  index_MH <- which(all_data$mz == closest_mz_MH & all_data$type == "Reference")  # Reference is [M+H]+
  max_intensity_MH <- max(abs(all_data$intensity[index_MH]))  # Use abs to get the height
  points(closest_mz_MH, -max_intensity_MH, pch = 17, col = "red", cex = 0.8)  # Upward red triangle for [M+H]+
}

# Add a blue downward-pointing triangle for the precursor in the [2M+H]+ spectrum
if (!is.na(closest_mz_Na)) {
  index_Na <- which(all_data$mz == closest_mz_Na & all_data$type == "Experimental")  # Experimental is [2M+H]+
  max_intensity_Na <- max(all_data$intensity[index_Na])  # Get the height for [2M+H]+
  points(closest_mz_Na, max_intensity_Na, pch = 17, col = "red", cex = 0.8)  # Downward red triangle for [2M+H]+
}

# Make the outline of the plot thicker
box(lwd = 10)  # Adjust the line width for the plot border




# Experimental spectrum plot
plot(exp_data$mz, exp_data$intensity, type = 'h', col = "#9400D3", 
     xlim = c(50, 220), ylim = c(0, max(exp_data$intensity) * 1.2),
     xlab = "m/z", ylab = "Relative Intensity", lwd = 17,
     main = "Experimental Spectrum [2M+H]+")

# Make the plot frame thicker
box(lwd = 25)  # Adjust the value of lwd for desired thickness

# Add annotations for the precursor
if (!is.na(closest_mz_Na)) {
  points(closest_mz_Na, max(exp_data$intensity), pch = 17, col = "red", cex = 1.2)  # Red triangle for precursor
  text(closest_mz_Na, max(exp_data$intensity) * 1.1, labels = "Precursor [2M+H]+", col = "red", pos = 3, cex = 0.8)
}




# Experimental spectrum plot
plot(exp_data$mz, exp_data$intensity, type = 'h', col = "#FF7F00", 
     xlim = c(50, 220), ylim = c(0, max(exp_data$intensity) * 1.2),
     xlab = "m/z", ylab = "Relative Intensity", lwd = 17,
     main = "Experimental Spectrum [M+H-H2O]+")

# Make the plot frame thicker
box(lwd = 25)  # Adjust the value of lwd for desired thickness

# Add annotations for the precursor
if (!is.na(closest_mz_Na)) {
  points(closest_mz_Na, max(exp_data$intensity), pch = 17, col = "red", cex = 1.2)  # Red triangle for precursor
  text(closest_mz_Na, max(exp_data$intensity) * 1.1, labels = "Precursor [M+H-H2O]+", col = "red", pos = 3, cex = 0.8)
}



# Experimental spectrum plot
plot(exp_data$mz, exp_data$intensity, type = 'h', col = "#FF0000", 
     xlim = c(50, 220), ylim = c(0, max(exp_data$intensity) * 1.2),
     xlab = "m/z", ylab = "Relative Intensity", lwd = 17,
     main = "Experimental Spectrum [M+Na]+")

# Make the plot frame thicker
box(lwd = 25)  # Adjust the value of lwd for desired thickness

# Add annotations for the precursor
if (!is.na(closest_mz_Na)) {
  points(closest_mz_Na, max(exp_data$intensity), pch = 17, col = "red", cex = 1.2)  # Red triangle for precursor
  text(closest_mz_Na, max(exp_data$intensity) * 1.1, labels = "Precursor [M+Na]+", col = "red", pos = 3, cex = 0.8)
}


# Experimental spectrum plot
plot(exp_data$mz, exp_data$intensity, type = 'h', col = "#00A859", 
     xlim = c(50, 220), ylim = c(0, max(exp_data$intensity) * 1.2),
     xlab = "m/z", ylab = "Relative Intensity", lwd = 17,
     main = "Experimental Spectrum [M+NH4]+")

# Make the plot frame thicker
box(lwd = 25)  # Adjust the value of lwd for desired thickness

# Add annotations for the precursor
if (!is.na(closest_mz_Na)) {
  points(closest_mz_Na, max(exp_data$intensity), pch = 17, col = "red", cex = 1.2)  # Red triangle for precursor
  text(closest_mz_Na, max(exp_data$intensity) * 1.1, labels = "Precursor [M+NH4]+", col = "red", pos = 3, cex = 0.8)
}


# Reference spectrum plot
plot(ref_data$mz, ref_data$intensity, type = 'h', col = "black",
     xlim = c(50, 220), ylim = c(0, max(ref_data$intensity) * 1.2),
     xlab = "m/z", ylab = "Relative Intensity", lwd = 17,
     main = "Reference Spectrum [M+H]+")

# Make the plot frame thicker
box(lwd = 25)  # Adjust the value of lwd for desired thickness

# Add annotations for the precursor
if (!is.na(closest_mz_MH)) {
  points(closest_mz_MH, max(ref_data$intensity), pch = 17, col = "blue", cex = 1.2)  # Upward blue triangle for precursor
  text(closest_mz_MH, max(ref_data$intensity) * 1.1, labels = "Precursor [M+H]+", col = "blue", pos = 3, cex = 0.8)
}














library(dplyr)

# Group by Adduct2 and filter rows with Adduct2 == "[M+Na]+"
filtered_data <- combined_data %>%
  group_by(Adduct2) %>%
  filter(Adduct2 == "[M+Na]+") %>%
  ungroup()

# View the filtered data
print(filtered_data)

filtered_data$Metabolites_CE <- paste0(filtered_data$Metabolite_Name, "_", filtered_data$Average_CE)

# Set metabolite name
metabolite_name <- "3-[(4-Fluorophenyl)sulfonyl]propanoic acid_21.87"

# Manually set the row index 'i'
i <- 1  # You can manually change this value to select the row
ms2_tol <- 0.02

# Extract and process whole_mz_cleaned_1 from the i-th matching row
whole_mz_cleaned_1_str <- filtered_data %>%
  filter(Metabolites_CE == metabolite_name) %>%
  slice(i) %>%  # Manually set 'i'
  pull(whole_mz_cleaned_1)
whole_mz_cleaned_1 <- as.numeric(strsplit(whole_mz_cleaned_1_str, ",")[[1]])

# Extract and process whole_int_cleaned_1 from the i-th matching row
whole_int_cleaned_1_str <- filtered_data %>%
  filter(Metabolites_CE == metabolite_name) %>%
  slice(i) %>%
  pull(whole_int_cleaned_1)
whole_int_cleaned_1 <- as.numeric(strsplit(whole_int_cleaned_1_str, ",")[[1]])

# Extract and process whole_mz_cleaned_2 from the i-th matching row
whole_mz_cleaned_2_str <- filtered_data %>%
  filter(Metabolites_CE == metabolite_name) %>%
  slice(i) %>%
  pull(whole_mz_cleaned_2)
whole_mz_cleaned_2 <- as.numeric(strsplit(whole_mz_cleaned_2_str, ",")[[1]])

# Extract and process whole_int_cleaned_2 from the i-th matching row
whole_int_cleaned_2_str <- filtered_data %>%
  filter(Metabolites_CE == metabolite_name) %>%
  slice(i) %>%
  pull(whole_int_cleaned_2)
whole_int_cleaned_2 <- as.numeric(strsplit(whole_int_cleaned_2_str, ",")[[1]])

# Normalize intensities
exp_fragInt <- whole_int_cleaned_2
ref_fragMzInt <- whole_int_cleaned_1

# Round normalized intensities to 2 decimals
exp_fragInt <- round(exp_fragInt, 2)
ref_fragMzInt <- round(ref_fragMzInt, 2)

# Create data frames for experimental and reference spectra
exp_data <- data.frame(mz = round(whole_mz_cleaned_2, 2), intensity = exp_fragInt, type = "Experimental")
ref_data <- data.frame(mz = round(whole_mz_cleaned_1, 2), intensity = ref_fragMzInt, type = "Reference")

# Extract precursor m/z values
precursor1 <- filtered_data %>%
  filter(Metabolites_CE == metabolite_name) %>%
  slice(i) %>%
  pull(precursor1)

precursor2 <- filtered_data %>%
  filter(Metabolites_CE == metabolite_name) %>%
  slice(i) %>%
  pull(precursor2)

# Create all_data data frame for plotting
all_data <- rbind(
  data.frame(mz = exp_data$mz, intensity = exp_data$intensity, type = "Experimental"),
  data.frame(mz = ref_data$mz, intensity = -ref_data$intensity, type = "Reference")
)

# Plot original data with custom colors
plot(all_data$mz, all_data$intensity, type = 'h', xlim = c(50, 160), 
     col = ifelse(all_data$type == "Experimental", "#EA4C46", "black"), 
     xlab = "m/z", ylab = "Relative Intensity", lwd = 17)


#### 15 eV #F6BDC0, 22 eV  #F1959B, 42 eV #EA4C46, 55 eV #FF0000


# Add black line at y = 0
abline(h = 0, col = "black", lwd = 20)

# Add text annotations for [M+NH4]+ and [M+H]+
text(x = min(all_data$mz) + 115, y = max(all_data$intensity) - 100, labels = "MS/MS of [M+Na]+", pos = 4, col = "blue", cex = 1)
text(x = min(all_data$mz) + 115, y = min(all_data$intensity) + 100, labels = "MS/MS of [M+H]+", pos = 4, col = "red", cex = 1)

# Function to find the closest match within tolerance
find_closest <- function(mz_vector, mz_value, tolerance) {
  abs_diff <- abs(mz_vector - mz_value)
  min_diff <- min(abs_diff)
  if (min_diff <= tolerance) {
    return(mz_vector[which.min(abs_diff)])
  } else {
    return(NA)
  }
}

# Find the closest m/z values for the precursor ions in each spectrum
closest_mz_MH <- find_closest(all_data$mz, precursor1, ms2_tol)
closest_mz_Na <- find_closest(all_data$mz, precursor2, ms2_tol)

# Add a red upward-pointing triangle for the precursor in the [M+H]+ spectrum
if (!is.na(closest_mz_MH)) {
  index_MH <- which(all_data$mz == closest_mz_MH & all_data$type == "Reference")  # Reference is [M+H]+
  max_intensity_MH <- max(abs(all_data$intensity[index_MH]))  # Use abs to get the height
  points(closest_mz_MH, -max_intensity_MH, pch = 17, col = "red", cex = 0.8)  # Upward red triangle for [M+H]+
}

# Add a blue downward-pointing triangle for the precursor in the [2M+H]+ spectrum
if (!is.na(closest_mz_Na)) {
  index_Na <- which(all_data$mz == closest_mz_Na & all_data$type == "Experimental")  # Experimental is [2M+H]+
  max_intensity_Na <- max(all_data$intensity[index_Na])  # Get the height for [2M+H]+
  points(closest_mz_Na, max_intensity_Na, pch = 17, col = "red", cex = 0.8)  # Downward red triangle for [2M+H]+
}

# Make the outline of the plot thicker
box(lwd = 20)  # Adjust the line width for the plot border
