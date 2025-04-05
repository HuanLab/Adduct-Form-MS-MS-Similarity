# Load libraries
library(openxlsx)
library(dplyr)

# Set working directory
setwd("X:/Users/Botao_Liu/Worklog/Adduct/20250320_New_start/NIST20")

# Define function to process and save adduct files
process_adducts <- function(input_file, adducts) {
  df <- read.xlsx(input_file)
  
  for (adduct in adducts) {
    filtered <- df %>% 
      filter(WHOLE_ADDUCT == adduct) %>%
      mutate(across(where(is.list), ~sapply(., toString)))
    
    output_filename <- paste0("adduct_", adduct, ".xlsx")
    write.xlsx(filtered, output_filename, rowNames = FALSE)
  }
}

# Negative adducts
neg_adducts <- c("[M-H]-", "[M+Cl]-", "[M-H-H2O]-", "[2M-H]-")
process_adducts("NIST20_NEG.xlsx", neg_adducts)

# Positive adducts
pos_adducts <- c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+", "[2M+H]+")
process_adducts("NIST20_POS.xlsx", pos_adducts)

cat("All adduct Excel files saved successfully.\n")