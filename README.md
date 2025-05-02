# Adduct-Based MS/MS Analysis and Visualization Pipeline

This repository contains a complete workflow for analyzing and visualizing adduct-specific MS/MS spectral similarity using R scripts.

---

##  Workflow Overview

1. **Step1_Adduct_Sorting.R**  
   Organizes MS/MS data by adduct types from input spectral libraries.

2. **Step2_Finding_the_Adduct_Pairs.R**  
   Identifies all valid adduct pairings based on metadata and collision energy matching.

3. **Step3_Similarity_Calculation.R**  
   Calculates cosine similarity (and optionally modified versions) between MS/MS spectra across adduct pairs.

Once Steps 1–3 are complete, you can generate specific figures using the scripts below.

---

## 📊 Figure Generation Scripts

- `Plot_for_fig1.R`: Overview or summary plot
- `Plot_for_fig2.R`: Density distribution of spectral similarities
- `Plot_for_fig3.R`: Match type categorization
- `Plot_for_fig5A_same_CE.R`: Cosine similarity at matched CE values
- `Plot_for_fig5B.R`: Similarity distributions for various adducts
- `Plot_for_fig5C_different_CE.R`: Spectral similarity under CE mismatch
- `Plot_for_fig6A_Direct_searching.R`: Performance of direct spectral searching
- `Plot_for_fig6B_Molecular_networking.R`: Molecular networking connectivity
- `Plot_for_fig6C_1_Model_Training_without_opt.R`: Model training without optimization
- `Plot_for_fig6C_2_Model_Training_with_opt.R`: Optimized model training pipeline
- `Plot_for_fig6C_3_Model_Comparisons.R`: ML-based model comparisons across settings

---
