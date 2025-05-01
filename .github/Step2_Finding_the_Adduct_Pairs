library(dplyr)
library(openxlsx)
library(data.table)
library(tidyr)
library(jsonlite)
library(writexl)
library(rcdk)
library(classyfireR)

# === PATHS ===
base_dir <- "X:/Users/Botao_Liu/Worklog/Adduct/20241029"
setwd(base_dir)

# === 1. Load and prepare data ===
pos_raw <- read.xlsx(file.path("C:/Program Files/Botao/Worklog/20240923_RF_POS", "NIST20_POS.xlsx"))
colnames(pos_raw)[c(5, 9)] <- c("WHOLE_Int", "WHOLE_PRECURSOR_mz")

# === Unique SMILES for classification ===
smiles_df <- distinct(pos_raw, WHOLE_SMILES, .keep_all = TRUE) %>% select(WHOLE_SMILES)

# Classify with classyfire
classify_smiles_compounds <- function(df) {
  df %>% rowwise() %>%
    mutate(
      classification = tryCatch({
        mol <- parse.smiles(WHOLE_SMILES)[[1]]
        inchikey <- get.inchi.key(mol)
        get_classification(inchikey)$classification
      }, error = function(e) NULL),
      superclass = classification$superclass$name,
      class = classification$class$name,
      subclass = classification$subclass$name
    ) %>% ungroup() %>% select(-classification)
}

classified_smiles <- classify_smiles_compounds(smiles_df)
write.csv(classified_smiles, "positive_dataframe_SMILES.csv", row.names = FALSE)

# === Merge classifications with full data ===
json_data <- read_json("ALL_NIST20_output.json")
json_df <- bind_rows(lapply(names(json_data), function(k) {
  val <- json_data[[k]]
  if (is.list(val) && length(val) == 3) {
    df <- as.data.frame(t(unlist(val)), stringsAsFactors = FALSE)
    df$Key <- k
    df
  }
})) %>% select(Key, superclass, class, subclass)

positive_df <- merge(pos_raw, json_df, by.x = "WHOLE_SMILES", by.y = "Key", all.x = TRUE)
write.csv(positive_df, "positive_dataframe.csv", row.names = FALSE)

# === 2. Clean the data ===
positive_df <- fread("positive_dataframe.csv") %>%
  mutate(
    WHOLE_COLLISION_ENERGY = as.numeric(WHOLE_COLLISION_ENERGY),
    WHOLE_Int = as.character(WHOLE_Int),
    WHOLE_mz_cleaned = Map(function(mz, int) {
      mz_vals <- as.numeric(unlist(strsplit(mz, ",")))
      int_vals <- as.numeric(unlist(strsplit(int, ",")))
      keep <- int_vals >= 30
      if (any(keep)) paste(mz_vals[keep], collapse = ",") else NA
    }, WHOLE_mz, WHOLE_Int),
    WHOLE_Int_cleaned = Map(function(mz, int) {
      int_vals <- as.numeric(unlist(strsplit(int, ",")))
      keep <- int_vals >= 30
      if (any(keep)) paste(int_vals[keep], collapse = ",") else NA
    }, WHOLE_mz, WHOLE_Int)
  ) %>%
  filter(
    !is.na(WHOLE_COLLISION_ENERGY),
    WHOLE_COLLISION_ENERGY >= 0 & WHOLE_COLLISION_ENERGY <= 60,
    WHOLE_ADDUCT %in% c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+", "[2M+H]+")
  ) %>%
  mutate(Index = row_number())

# === 3. Find Adduct Pairs ===
adducts <- c("[M+K]+", "[M+NH4]+", "[M+Na]+", "[M+H-H2O]+", "[2M+H]+")
adduct_H <- positive_df %>% filter(WHOLE_ADDUCT == "[M+H]+") %>% as.data.table()

find_adduct_pairs <- function(adduct_name) {
  adduct_O <- fread(paste0("adduct_", gsub("[\\[\\]+-]", "", adduct_name), ".xlsx")) %>% as.data.table()
  common_names <- intersect(adduct_H$WHOLE_COMPOUND_NAME, adduct_O$WHOLE_COMPOUND_NAME)
  df_H <- adduct_H[WHOLE_COMPOUND_NAME %in% common_names]
  df_O <- adduct_O[WHOLE_COMPOUND_NAME %in% common_names]
  
  pairs <- merge(df_H, df_O, by = "WHOLE_COMPOUND_NAME", suffixes = c("_H", "_O")) %>%
    filter(abs(WHOLE_COLLISION_ENERGY_H - WHOLE_COLLISION_ENERGY_O) < 3) %>%
    mutate(
      Adduct1 = "[M+H]+", Adduct2 = adduct_name,
      CE1 = WHOLE_COLLISION_ENERGY_H, CE2 = WHOLE_COLLISION_ENERGY_O
    ) %>%
    select(
      Metabolite_Name = WHOLE_COMPOUND_NAME,
      CE1, CE2, Adduct1, Adduct2,
      Index_1 = Index_H, Index_2 = Index_O,
      precursor1 = WHOLE_PRECURSOR_mz_H, precursor2 = WHOLE_PRECURSOR_mz_O,
      whole_mz1 = WHOLE_mz_H, whole_mz2 = WHOLE_mz_O,
      whole_int1 = WHOLE_Int_H, whole_int2 = WHOLE_Int_O,
      whole_mz_cleaned_1, whole_mz_cleaned_2,
      whole_int_cleaned_1, whole_int_cleaned_2,
      INSTRUMENT_TYPE_1 = INSTRUMENT_TYPE_H, INSTRUMENT_TYPE_2 = INSTRUMENT_TYPE_O
    )
  
  write.csv(pairs, paste0(adduct_name, "_pairs_df.csv"), row.names = FALSE)
  cat(adduct_name, ": Pairs found =", nrow(pairs), "\n")
  return(pairs)
}

lapply(adducts, find_adduct_pairs)

# === 4. Filter by Instrument Type ===
lapply(adducts, function(adduct) {
  df <- fread(paste0(adduct, "_pairs_df.csv"))
  df_filtered <- df[df$INSTRUMENT_TYPE_1 == df$INSTRUMENT_TYPE_2]
  write.csv(df_filtered, paste0(adduct, "_filtered_pairs_df.csv"), row.names = FALSE)
  cat(adduct, ": Filtered =", nrow(df_filtered), "\n")
})
