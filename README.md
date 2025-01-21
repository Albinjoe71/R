# R Mutational Signature Analysis Using Maftools  

This repository contains a comprehensive analysis of mutational signatures in genomic data using the **maftools** package in R. The workflow includes reading MAF files, generating sample and gene summaries, visualizing mutation data, performing transitions and transversions analysis, and estimating mutational signatures.  

---

## 🧬 Key Analysis Steps
1. **Load MAF Data**: Read MAF file using `read.maf()` function from `maftools`.  
2. **Sample and Gene Summaries**: Obtain sample and gene summaries using `getSampleSummary()` and `getGeneSummary()`.  
3. **Clinical Data**: Access and review associated clinical data using `getClinicalData()`.  
4. **MAF Summary**: Write out MAF summary with `write.mafSummary()`.  
5. **Visualization**:  
   - MAF summary visualization using `plotmafSummary()`.  
   - Oncoplot of top 10 mutated genes using `oncoplot()`.  
6. **Transition & Transversion (TiTv)**: Analyze and visualize transitions and transversions with `titv()`.  
7. **Somatic Interactions**: Identify exclusive or co-occurring events in top mutated genes using `somaticInteractions()`.  
8. **Mutational Signatures**:  
   - Estimate mutational signatures using `estimateSignatures()`.  
   - Compare with COSMIC signatures (legacy and updated).  
9. **Heatmaps**: Visualize cosine similarities of mutational signatures using `pheatmap()`.  
10. **Final Signature Plot**: Plot final mutational signatures.

---

## 📂 Files Included  
- **R Script**: Code implementing the analysis (`mutational_signature_analysis.R`).  
- **Example MAF File**: Sample MAF file to test the code (optional).  

---

## 🛠️ Getting Started  

1. Clone the repository:  
   ```bash
   git clone https://github.com/Albinjoe71/R-Mutational-Signature-Analysis.git
