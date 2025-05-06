# Adduct-Based MS/MS Analysis and Visualization Pipeline

This repository contains a complete workflow for analyzing and visualizing adduct-specific MS/MS spectral similarity using R scripts.

---

##  Workflow Overview

1. **Step1_Adduct_Sorting.R**  
   Organizes MS/MS data by adduct types from input spectral libraries.

2. **Step2_Finding_the_Adduct_Pairs.R**  
   Identifies all valid adduct pairings based on identical chemical structures, the same or similar MS instrument type, and closely matched CEs (within 3 eV).

3. **Step3_Similarity_Calculation.R**  
   Calculates modified cosine similarity between MS/MS spectra across adduct pairs.

Once Steps 1–3 are complete, you can generate specific figures using the scripts below.

---

## 📊 Figure Generation Scripts

- `Plot_for_fig1.R`: Example MS/MS spectral comparison shows two key factors: adducts and CE
- `Plot_for_fig2.R`: Spectral similarity analysis: Density plot of modified cosine similarity of MS/MS comparison
- `Plot_for_fig3.R`: Matched fragment ratios and fragment matched types analysis
- `Plot_for_fig5A.R`: The trend between modified cosine similarity and CE
- `Plot_for_fig5B.R`: The trend between modified cosine similarity and ΔCE
- `Plot_for_fig5C.R`: The trend between mean geometric m/z center shift values and CE levels
- `Plot_for_fig6A.R`: Direct spectral searching annotation results impaired by mismatched adduct forms
- `Plot_for_fig6B.R`: Molecular networking based annotation results impaired by mismatched adduct forms
- `Plot_for_fig6C_Step1.R`: Effects of adduct forms and CE thresholds on detecting PO4-related structural indicators in MS/MS
- `Plot_for_fig6C_Step2.R`: TPR and FPR of ML-based annotation model performance comparisons with/without excluding alkali adduct and applying CE thresholds
- `Plot_for_fig6C_Step3.R`: ML annotation AUC curve comparisons with/without excluding alkali adduct and applying CE thresholds

---
