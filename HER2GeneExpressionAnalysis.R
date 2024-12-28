#All library imports
library(stringr)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(pathview)
library(pheatmap)
library(glmnet)
library(survival)

###2

#untar folder and extract files
folder = 'C:/Users/LENOVO/Downloads/brca_tcga_pan_can_atlas_2018.tar.gz'
untar(folder)

#set new directory
new_dir = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/" )
setwd(new_dir)
#print list of files in directory
list.files(new_dir)

###3
#read RNA-seq file
rna_seq = read.delim('data_mrna_seq_v2_rsem.txt', header = TRUE)

###4
#read patient data
patient_data = read.delim('data_clinical_patient.txt', skip=4, header=TRUE)

###5
#read copy number aberrations data
cna_data = read.delim('data_cna.txt', header=TRUE)

###6
#adjust column names so they all have same format of patient ID
colnames(rna_seq)[3:ncol(rna_seq)] = gsub("\\.", "-", colnames(rna_seq)[3:ncol(rna_seq)])
colnames(rna_seq)[3:ncol(rna_seq)] = gsub("-01$", "", colnames(rna_seq)[3:ncol(rna_seq)])

colnames(cna_data)[3:ncol(cna_data)] = gsub("\\.", "-", colnames(cna_data)[3:ncol(cna_data)])
colnames(cna_data)[3:ncol(cna_data)] = gsub("-01$", "", colnames(cna_data)[3:ncol(cna_data)])

rna_seq_patient_ids = colnames(rna_seq)[3:ncol(rna_seq)]
cna_patient_ids = colnames(cna_data)[3:ncol(cna_data)]
clinical_patient_ids = patient_data$PATIENT_ID

#find the common patient IDs across three datasets
common_patient_ids = Reduce(intersect, list(rna_seq_patient_ids, cna_patient_ids, clinical_patient_ids))

#see what columns are in common and reorder
rna_seq_columns = which(colnames(rna_seq) %in% common_patient_ids)
rna_seq = rna_seq[, c(1, 2, rna_seq_columns)]  # Keep the first two columns (Hugo_Symbol, Entrez_Gene_Id)

#do the same for cna and patient data
cna_columns = which(colnames(cna_data) %in% common_patient_ids)
cna_data = cna_data[, c(1, 2, cna_columns)]  # Keep the first two columns (Hugo_Symbol, Entrez_Gene_Id)
patient_data = patient_data[patient_data$PATIENT_ID %in% common_patient_ids, ]


###7
#make count matric of rna data so the dimensions can be used for the metadata & assign row names
assay = round(as.matrix(rna_seq[, -c(1,2)])) 
rownames(assay) = rna_seq[,1]

#find erbb2 row index
erbb2_idx = which(cna_data$Hugo_Symbol=="ERBB2")
#isolate that row in variable
erbb2_cna = cna_data[erbb2_idx, -c(1,2)]
#make function that assigns 1 for amplified and 0 for not amplified based on erbb2 variable
amplified = as.numeric(erbb2_cna>0)
#create metadata matrix
cna_metadata = matrix(0, dim(assay)[2], ncol=1)
cna_metadata[,1] = (amplified)
#name rows and columns
rownames(cna_metadata) = colnames(cna_data[erbb2_idx, -c(1,2)])
colnames(cna_metadata) = c("ERBB2_Amplified") 

###8
#impute nas in count matrix with 0s and negative values with 0s
assay[is.na(assay)] = 0  
assay[assay<0] = 0

#filter out genes with too many missing values
keep = rowSums(assay >= 10) >= 3
assay = assay[keep,]

#normalise
cna_metadata = as.data.frame(cna_metadata)
cna_metadata$ERBB2_Amplified = as.factor(cna_metadata$ERBB2_Amplified)
dds =  DESeqDataSetFromMatrix(countData = assay,
                              colData = cna_metadata,
                              design = ~ERBB2_Amplified)

dds = DESeq(dds)


###9

res = results(dds)
nrow(res)

#get genes that have a significant p value 
res_sig = res[which( res$padj < 0.05), ]

nrow(res_sig)


#10 most differentially expressed genes based on Fold Change
res_sig_10 = head(res_sig[order(abs(res_sig$log2FoldChange), decreasing = TRUE),], 10)

rownames(res_sig_10)
res_sig_10

###10


# separate into over and under expressed using log2foldchange
DE_over = rownames(res_sig[res_sig$log2FoldChange>0,])
DE_under = rownames(res_sig[res_sig$log2FoldChange<0,])

go_results_over = enrichGO(
  gene          = DE_over,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

print(head(go_results_over))
dotplot(go_results_over, showCategory=10) + ggtitle("Gene Ontology Enrichment Over Expressed")

go_results_under = enrichGO(
  gene          = DE_under,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# print and plot results
print(head(go_results_under))
dotplot(go_results_under, showCategory=10) + ggtitle("Gene Ontology Enrichment Under Expressed")

gene_entrez_over = bitr(
  DE_over,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

gene_entrez_under = bitr(
  DE_under,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

kegg_results_over =  enrichKEGG(
  gene          = gene_entrez_over[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

kegg_results_under =  enrichKEGG(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

print(head(kegg_results_over))

dotplot(kegg_results_over, showCategory=10) + ggtitle("Kegg Pathway Enrichment Over Expressed")

print(head(kegg_results_under))

dotplot(kegg_results_under, showCategory=10) + ggtitle("Kegg Pathway Enrichment Under Expressed")

###############################
#enrichment pathway analysis
reactome_results_over =  enrichPathway(
  gene          = gene_entrez_over[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)

reactome_results_under =  enrichPathway(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)


print(head(reactome_results_over))

dotplot(reactome_results_over, showCategory=10) + ggtitle("Reactome Pathway Enrichment Over Expressed")

print(head(reactome_results_under, 10))

dotplot(reactome_results_under, showCategory=10) + ggtitle("Reactome Pathway Enrichment Under Expressed")



go_results_over_pw = pairwise_termsim(go_results_over)
treeplot(go_results_over_pw)+ ggtitle("GO Enrichment Pairwise Comparison Over Expressed")

go_results_under_pw = pairwise_termsim(go_results_under)
treeplot(go_results_under_pw)+ ggtitle("GO Enrichment Pairwise Comparison Under Expressed")

kegg_results_over_pw = pairwise_termsim(kegg_results_over)
treeplot(kegg_results_over_pw)+ ggtitle("KEGG Enrichment Over Expressed")

kegg_results_under_pw = pairwise_termsim(kegg_results_under)
treeplot(kegg_results_under_pw)+ ggtitle("KEGG Enrichment Under Expressed")
####################################################

###11

vst = vst(dds)

###12
###13
#PCA graph
plotPCA(vst, intgroup=c("ERBB2_Amplified")) + ggtitle("ERBB2 Amplified PCA Plot")

#subset dataset on differentially expressed genes 
top_DE = order(res$padj)

vst_DE = assay(vst)[top_DE[1:10],]

pheatmap(
  vst_DE,
  cluster_rows = TRUE,      
  cluster_cols = TRUE,  
  scale = 'row',
  show_colnames = FALSE,
  show_rownames = TRUE) + ggtitle("ERBB2 Amplified Heatmap")


###14
deg_genes = rownames(res_sig)  
vst_sig_values = assay(vst(dds))[deg_genes, ] 

survival = data.frame(time=patient_data$OS_MONTHS, status=as.integer(substr(patient_data$OS_STATUS, 1, 1)))


#replace 0 values with small positive value so cox model works
survival$time[which(survival$time <=0)] = 1e-10

#split x=significant values and y=survival status & time
x = t(vst_sig_values)
y = Surv(survival$time, survival$status)


#fit cross-validation Lasso Cox Model to find optimal lambda
cv_fit = cv.glmnet(x, y, family="cox", nfolds = 10)
best_lambda = cv_fit$lambda.min

#get selected genes
coef = as.matrix(coef(cv_fit, s = "lambda.min"))
#display them in descending order based on coefficient and exclude genes with coefficient 0
selected_genes = coef[coef!=0,,drop=FALSE]
selected_genes = selected_genes[order(abs(selected_genes), decreasing = TRUE), , drop = FALSE]
print(selected_genes)

#get risk scores & separate risk groups into High and Low
risk_scores = predict(cv_fit, newx = x, s = "lambda.min", type = "link")
risk_group = ifelse(risk_scores > median(risk_scores), "High", "Low")

#fit km line
km_fit = survfit(Surv(survival$time, survival$status) ~ risk_group)
#plot survival line
plot(km_fit, col = c("red", "blue"), xlab = "Time (months)", ylab = "Survival Probability")
legend("bottomleft", legend = c("High Risk", "Low Risk"), col = c("red", "blue"), lty = 1)

