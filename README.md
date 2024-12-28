# HER2-gene-expression-analysis
This is a gene expression analysis of HER2+ breast cancer.

Dataset required:  https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
Files: data_mrna_seq_v2_rsem.txt (RNA-seq gene expression data), data_clinical_patient.txt (Patient survival and clinical data), data_cna.txt (Copy number aberration data)

Packages required: "DESeq2", "clusterProfiler", "enrichplot", "org.Hs.eg.db", "ReactomePA", "pathview", "pheatmap", "glmnet", "survival", "stringr"

Libraries and Data loading: lines 1-37 
- Load the 3 files from the dataset
  
Data preprocessing: lines 39-90.
- Patient ID formatting and sort order of patients the same throughout 3 files
  
Statistical analysis for identifying differentially expressed genes: lines 91-116
- run DESeq analysis on data & identify 10 most differentially expressed genes based on fold change
  
Over and under expressed genes pathway enrichment analyses: lines 121-230
- Seperate genes for over and under expression and run enrichGo, enrichKEGG, enrichPathway and pairwise_termsim.
- Plot dotplots for each analysis
  
Variance stabilised gene expressions visualisation: lines 234-253
- get vst values and plot PCA plot and heatmap
   
Survival model: lines 255-290.
- calculate survival, run cross-validation on fitted Lasso Cox Model
- print sorted genes based on descending coefficent attributed 
- separate risk groups based on survival median and plot Kaplan-Meier line

