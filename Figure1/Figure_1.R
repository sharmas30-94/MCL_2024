# Main Figure 1 ###

#Figure 1a
# All mutations have been included based on these four initial criteria:

# 1: Variants marked as "." in dbSNP (non-flagged) AND called by Varscan = "Y" 
#    AND also detected (or called) in Mutect.

# 2: Variants in TP53 marked as non-flagged in dbSNP AND called by Varscan = "Y" 
#    AND detected in Mutect.

# 3: Variants marked as "." in dbSNP (non-flagged) AND called by Mutect = "Y" 
#    AND detected in Varscan, AND have "PASS" in the VCFFILTER.Mutect2 field.

# 4: Variants with WES.non_cancer_AF_popmax < 0.01 (1%)
#    AND VAF.Mutect2.Percentage > 5 
#    AND VAF.Varscan2.Percentage > 5.

# Initial setup
setwd("/Users/sharmas30/Desktop/code for MCL/Github_code/Figure1/Data")

# Install and load required packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("maftools", update = FALSE)
library(maftools)
library("survival")
library("survminer")
library("dplyr")

# Step 0: Apply next set of filters
data1 <- read.table("WES_130_cases.txt", header = TRUE, sep = "\t")
data1 <- subset(
   data1,
   data1$WES.non_cancer_AF_popmax < 0.01 &
     (data1$VAF.Mutect2.Percenatge > 5 & data1$VAF.Varscan2.Percentage > 5))
 write.table(data1, file="WES_130_cases_filtered_mutations.txt", sep="\t", quote=FALSE, row.names=FALSE)


# Step 1: Define input files for NAMC and BC datasets
namc_var_file <- "WES_130_cases_filtered_mutations.txt"
namc_clin_file <- "WES_130_clinical.tsv"

bc_var_file <- "BC_23_cases.txt"
bc_clin_file <- "BC_23_clinical.tsv"

# Step 2: Convert ANNOVAR data to MAF format and read clinical data
namc_maf <- annovarToMaf(
  annovar = namc_var_file, 
  Center = "CSI-NUS", 
  refBuild = "hg38", 
  tsbCol = "Tumor_Sample_Barcode", 
  table = "refGene"
)

namc <- read.maf(maf = namc_maf, clinicalData = namc_clin_file)

bc_maf <- annovarToMaf(
  annovar = bc_var_file, 
  Center = "CSI-NUS", 
  refBuild = "hg38", 
  tsbCol = "Tumor_Sample_Barcode", 
  table = "refGene"
)

bc <- read.maf(maf = bc_maf, clinicalData = bc_clin_file)

# Step 3: Merge NAMC and BC  for Cohort-1
merged_maf <- merge_mafs(maf = c(namc, bc), verbose = TRUE)
#Remove all Flags shown to be highly mutated across cancers
#https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-014-0064-y

merged_maf<- filterMaf(merged_maf, genes = c('TTN', 'MUC3A', 'GXYLT1', 'MUC16', 'TTC14', 'MUC6', 'CFTR', 'MUC5AC', 'APOB', 'ADGRV1', 'USH2A', 'MYO7A', 'MUC17', 'MUC4', 'USH2A', 'ZNF343', 'RYR2', 'ZNF429', 'OBSCN', 'PALLD',  'TTCN14', 'ZFYVE26', 'ZNF91', 'DIS3', 'MTO1', 'SRCAP', 'ABCA4', 'DNAH12', 'DYNC111', 'TTCN14', 'FLNA', 'RBP3', 'AR', 'BAG6', 'CCDC168', 'ZNF107', 'ZNF708', 'LAMB1', 'PLEC', 'RAI1', 'TSFM', 'APOE', 'CDK5RAP2', 'CDH7', 'CNTNAP5', 'DST', 'FLRT2', 'GJB2', 'ITPR1', 'KCNT1', 'MYO15A', 'NRXN2', 'PCDH15', 'PDSS1', 'PRKN', 'RNF17', 'SLC22A5', 'TENM4', 'TRPM4', 'ZNF254', 'ZNF493', 'PRR12', 'ASB10', 'GOLGB1', 'CACNA1E', 'ANKRD31', 'CDH23', 'HMCN1', 'CSMD3', 'DNAH17', ''))

# Step 4: Define color palettes for clinical annotations

#Genes retained after removing flagged genes (Gene Mutation frequency <= 100 in population database) , flagged gene included only if present in previous MCL studies. 
k<-read.table("genes.of.interest.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)


fabcolors1 = c("blue", "yellow", "grey")
names(fabcolors1) = c("Positive", "Negative", "N/A")
fabcolors1 = list(SOX11_IHC = fabcolors1)


fabcolors2 = c("red", "cyan", "grey")
names(fabcolors2) = c("Positive", "Negative", "N/A")
fabcolors2 = list(TP53_IHC = fabcolors2)

fabcolors3 = c("green", "purple", "grey")
names(fabcolors3) = c("Positive", "Negative", "N/A")
fabcolors3 = list(Ki67_Counts = fabcolors3)


fabcolors4 = c("plum", "tan3")
names(fabcolors4) = c("FF", "FFPE")
fabcolors4 = list(BIOPSY_TYPE = fabcolors4)


fabcolors5 = c("deepskyblue", "pink", "grey")
names(fabcolors5) = c("Male", "Female", "N/A")
fabcolors5 = list(SEX = fabcolors5)


fabcolors6 = c("darkorange1", "gold1", "grey")
names(fabcolors6) = c("Yes", "No", "N/A")
fabcolors6 = list(RITUXIMAB_ADMINISTERED = fabcolors6)


fabcolors7 = c("yellow4", "deeppink4", "grey")
names(fabcolors7) = c("Early", "Late", "N/A")
fabcolors7 = list(AAS = fabcolors7)

fabcolors8 = c("lightcoral", "lightgreen", "grey")
names(fabcolors8) = c(">=60", "<60", "N/A")
fabcolors8 = list(AGE_Status = fabcolors8)

# Step 5: Generate an Oncoplot
pdf("Figure_1A.pdf", height=20, width=15)
oncoplot(
  maf = merged_maf,
  fontSize = 0.4,
  sepwd_genes = 1,
  pathways = k,
  clinicalFeatures = c('BIOPSY_TYPE', 'SOX11_IHC', 'TP53_IHC', 'Ki67_Counts', 'RITUXIMAB_ADMINISTERED', 'AGE_Status', 'AAS', 'SEX'), 
  annotationColor = c(fabcolors4, fabcolors1, fabcolors2, fabcolors3, fabcolors6, fabcolors7, fabcolors8, fabcolors5),
  anno_height = 2,
  gene_mar = 9,
  drawBox = TRUE,
  sepwd_samples = 0.8,
  legend_height = 2,
  logColBar =FALSE, 
  
)
dev.off()

################ Figure 1B ###############
namc_var_file <- "WES_130_cases_filtered_mutations.txt"
namc_clin_file <- "WES_130_clinical.tsv"

bc_var_file <- "BC_23_cases.txt"
bc_clin_file <- "BC_23_clinical.tsv"

# Step 2: Convert ANNOVAR data to MAF format and read clinical data
namc_maf <- annovarToMaf(
  annovar = namc_var_file, 
  Center = "CSI-NUS", 
  refBuild = "hg38", 
  tsbCol = "Tumor_Sample_Barcode", 
  table = "refGene"
)

namc <- read.maf(maf = namc_maf, clinicalData = namc_clin_file)

bc_maf <- annovarToMaf(
  annovar = bc_var_file, 
  Center = "CSI-NUS", 
  refBuild = "hg38", 
  tsbCol = "Tumor_Sample_Barcode", 
  table = "refGene"
)

bc <- read.maf(maf = bc_maf, clinicalData = bc_clin_file)

# Step 3: Merge NAMC and BC  for Cohort-1
merged_maf <- merge_mafs(maf = c(namc, bc), verbose = TRUE)
#Remove all Flags shown to be highly mutated across cancers
#https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-014-0064-y

merged_maf<- filterMaf(merged_maf, genes = c('TTN', 'MUC3A', 'GXYLT1', 'MUC16', 'TTC14', 'MUC6', 'CFTR', 'MUC5AC', 'APOB', 'ADGRV1', 'USH2A', 'MYO7A', 'MUC17', 'MUC4', 'USH2A', 'ZNF343', 'RYR2', 'ZNF429', 'OBSCN', 'PALLD',  'TTCN14', 'ZFYVE26', 'ZNF91', 'DIS3', 'MTO1', 'SRCAP', 'ABCA4', 'DNAH12', 'DYNC111', 'TTCN14', 'FLNA', 'RBP3', 'AR', 'BAG6', 'CCDC168', 'ZNF107', 'ZNF708', 'LAMB1', 'PLEC', 'RAI1', 'TSFM', 'APOE', 'CDK5RAP2', 'CDH7', 'CNTNAP5', 'DST', 'FLRT2', 'GJB2', 'ITPR1', 'KCNT1', 'MYO15A', 'NRXN2', 'PCDH15', 'PDSS1', 'PRKN', 'RNF17', 'SLC22A5', 'TENM4', 'TRPM4', 'ZNF254', 'ZNF493', 'PRR12', 'ASB10', 'GOLGB1', 'CACNA1E', 'ANKRD31', 'CDH23', 'HMCN1', 'CSMD3', 'DNAH17', ''))


WES<-mutCountMatrix(
  merged_maf,
  includeSyn = FALSE,
  countOnly = NULL,
  removeNonMutated = FALSE
)
write.csv(WES, file="Cohort-1-mutation-matrix.csv")

tar_var_file <- "Targeted_137_cases_filtered_mutations.txt"
tar_clin_file <- "Targeted_137_clinical.tsv"

# Step 2: Convert ANNOVAR data to MAF format and read clinical data
tar_maf <- annovarToMaf(
  annovar = tar_var_file, 
  Center = "CSI-NUS", 
  refBuild = "hg38", 
  tsbCol = "Tumor_Sample_Barcode", 
  table = "refGene"
)

tar_maf <- read.maf(maf = tar_maf, clinicalData = tar_clin_file)

merged_maf1<- filterMaf(tar_maf, genes = c('TTN', 'MUC3A', 'GXYLT1', 'MUC16', 'TTC14', 'MUC6', 'CFTR', 'MUC5AC', 'APOB', 'ADGRV1', 'USH2A', 'MYO7A', 'MUC17', 'MUC4', 'USH2A', 'ZNF343', 'RYR2', 'ZNF429', 'OBSCN', 'PALLD',  'TTCN14', 'ZFYVE26', 'ZNF91', 'DIS3', 'MTO1', 'SRCAP', 'ABCA4', 'DNAH12', 'DYNC111', 'TTCN14', 'FLNA', 'RBP3', 'AR', 'BAG6', 'CCDC168', 'ZNF107', 'ZNF708', 'LAMB1', 'PLEC', 'RAI1', 'TSFM', 'APOE', 'CDK5RAP2', 'CDH7', 'CNTNAP5', 'DST', 'FLRT2', 'GJB2', 'ITPR1', 'KCNT1', 'MYO15A', 'NRXN2', 'PCDH15', 'PDSS1', 'PRKN', 'RNF17', 'SLC22A5', 'TENM4', 'TRPM4', 'ZNF254', 'ZNF493', 'PRR12', 'ASB10', 'GOLGB1', 'CACNA1E', 'ANKRD31', 'CDH23', 'HMCN1', 'CSMD3', 'DNAH17'))


Targeted<-mutCountMatrix(
  merged_maf1,
  includeSyn = FALSE,
  countOnly = NULL,
  removeNonMutated = FALSE
)

write.csv(Targeted, file="Cohort-2-mutation-matrix.csv")

# Step 1: Read data files
d1 <- read.csv("Cohort-1-mutation-matrix.csv", header = TRUE, sep = ",", check.names = FALSE)
d2 <- read.csv("Cohort-2-mutation-matrix.csv", header = TRUE, sep = ",", check.names = FALSE)
colnames(d1)[1] <- "Genes"
colnames(d2)[1] <- "Genes"

# Step 2: Merge the two datasets by the "Genes" column
merged_data <- merge(d1, d2, by = "Genes", all = TRUE)

# Step 3: Remove rows with NA values
merged_data <- na.omit(merged_data)

# Step 4: Set the "Genes" column as row names for transposing
rownames(merged_data) <- merged_data$Genes
merged_data$Genes <- NULL  # Remove the "Genes" column as it is now part of row names

# Step 5: Transpose the data
transposed_data <- as.data.frame(t(merged_data))

# Step 6: Print the transposed data
print("Transposed Data:")
print(transposed_data)

# Step 7: Combine CLinical Information
j1<-read.table("WES_130_clinical.tsv", header=TRUE, sep="\t")
j2<-read.table("BC_23_clinical.tsv", header=TRUE, sep="\t")
j3<-read.table("Targeted_137_clinical.tsv", header=TRUE, sep="\t")
com<-rbind(j1, j2, j3)

# Step 8: Add all clinical information into transposed matrix

m<-match(rownames(transposed_data), com[,1])
com<-com[m,]
f<-cbind(transposed_data, com)


slots<-list()
i=14 #ATM
df <- select(f, colnames(f[i]), OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0, RITUXIMAB_ADMINISTERED)
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'MI'
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'MC'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'WTI'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'WTC'
df[df == "N/A"] <- NA
df <- na.omit(df)
df$OVERALL_SURVIVAL_MONTHS <- as.numeric(as.character(df$OVERALL_SURVIVAL_MONTHS))
df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0 <- as.numeric(as.character(df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0))
fit <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ Mut_Status, data=df)
slots[[1]]<-ggsurvplot(fit, pval = TRUE, data = df, risk.table = TRUE, palette = c("blue", "red", "#8DCFEC", "#FF8A8A"), title = names(f[i]), surv.median.line = "hv", risk.table.legend = FALSE, risk.table.col = "strata")
pairwise_survdiff(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ Mut_Status, data=df)
ggsurvplot(fit, pval = TRUE, data = df, risk.table = TRUE, palette = c("blue", "red", "#8DCFEC", "#FF8A8A"), title = names(f[i]), surv.median.line = "hv", risk.table.legend = FALSE, risk.table.col = "strata")


i=188 #TP53
df <- select(f, colnames(f[i]), OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0, RITUXIMAB_ADMINISTERED)
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'MI'
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'MC'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'WTI'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'WTC'
df[df == "N/A"] <- NA
df <- na.omit(df)
df$OVERALL_SURVIVAL_MONTHS <- as.numeric(as.character(df$OVERALL_SURVIVAL_MONTHS))
df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0 <- as.numeric(as.character(df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0))
fit <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ Mut_Status, data=df)
slots[[2]]<-ggsurvplot(fit, pval = TRUE, data = df, risk.table = TRUE, palette = c("blue", "red", "#8DCFEC", "#FF8A8A"), title = names(f[i]), surv.median.line = "hv", risk.table.legend = FALSE, risk.table.col = "strata")
pairwise_survdiff(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ Mut_Status, data=df)

i=94 #KMT2C
df <- select(f, colnames(f[i]), OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0, RITUXIMAB_ADMINISTERED)
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'MI'
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'MC'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'WTI'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'WTC'
df[df == "N/A"] <- NA
df <- na.omit(df)
df$OVERALL_SURVIVAL_MONTHS <- as.numeric(as.character(df$OVERALL_SURVIVAL_MONTHS))
df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0 <- as.numeric(as.character(df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0))
fit <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ Mut_Status, data=df)
slots[[3]]<-ggsurvplot(fit, pval = TRUE, data = df, risk.table = TRUE, palette = c("blue", "red", "#8DCFEC", "#FF8A8A"), title = names(f[i]), surv.median.line = "hv", risk.table.legend = FALSE, risk.table.col = "strata")
pairwise_survdiff(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ Mut_Status, data=df)


i=29 #CCND1
df <- select(f, colnames(f[i]), OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0, RITUXIMAB_ADMINISTERED)
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'MI'
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'MC'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'WTI'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'WTC'
df[df == "N/A"] <- NA
df <- na.omit(df)
df$OVERALL_SURVIVAL_MONTHS <- as.numeric(as.character(df$OVERALL_SURVIVAL_MONTHS))
df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0 <- as.numeric(as.character(df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0))
fit <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ Mut_Status, data=df)
slots[[4]]<-ggsurvplot(fit, pval = TRUE, data = df, risk.table = TRUE, palette = c("blue", "red", "#8DCFEC", "#FF8A8A"), title = names(f[i]), surv.median.line = "hv", risk.table.legend = FALSE, risk.table.col = "strata")
pairwise_survdiff(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ Mut_Status, data=df)

pdf("Figure_1B.pdf", height = 18, width = 30)
arrange_ggsurvplots(slots, print = TRUE,
                    ncol = 2, nrow = 2, risk.table.height = 0.4)
dev.off()

#PFS for selected genes
slots<-list()
i=14 #ATM
df <- select(f, colnames(f[i]), PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0, RITUXIMAB_ADMINISTERED)
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'MI'
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'MC'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'WTI'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'WTC'
df[df == "N/A"] <- NA
df <- na.omit(df)
df$PROGRESSION_FREE_SURVIVAL_MONTHS <- as.numeric(as.character(df$PROGRESSION_FREE_SURVIVAL_MONTHS))
df$PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0 <- as.numeric(as.character(df$PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0))
fit <- survfit(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ Mut_Status, data=df)
slots[[1]]<-ggsurvplot(fit, pval = TRUE, data = df, risk.table = TRUE, palette = c("blue", "red", "#8DCFEC", "#FF8A8A"), title = names(f[i]), surv.median.line = "hv", risk.table.legend = FALSE, risk.table.col = "strata")
pairwise_survdiff(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ Mut_Status, data=df)



i=188 #TP53
df <- select(f, colnames(f[i]), PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0, RITUXIMAB_ADMINISTERED)
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'MI'
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'MC'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'WTI'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'WTC'
df[df == "N/A"] <- NA
df <- na.omit(df)
df$PROGRESSION_FREE_SURVIVAL_MONTHS <- as.numeric(as.character(df$PROGRESSION_FREE_SURVIVAL_MONTHS))
df$PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0 <- as.numeric(as.character(df$PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0))
fit <- survfit(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ Mut_Status, data=df)
slots[[2]]<-ggsurvplot(fit, pval = TRUE, data = df, risk.table = TRUE, palette = c("blue", "red", "#8DCFEC", "#FF8A8A"), title = names(f[i]), surv.median.line = "hv", risk.table.legend = FALSE, risk.table.col = "strata")
pairwise_survdiff(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ Mut_Status, data=df)


i=94 #KMT2C
df <- select(f, colnames(f[i]), PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0, RITUXIMAB_ADMINISTERED)
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'MI'
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'MC'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'WTI'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'WTC'
df[df == "N/A"] <- NA
df <- na.omit(df)
df$PROGRESSION_FREE_SURVIVAL_MONTHS <- as.numeric(as.character(df$PROGRESSION_FREE_SURVIVAL_MONTHS))
df$PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0 <- as.numeric(as.character(df$PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0))
fit <- survfit(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ Mut_Status, data=df)
slots[[3]]<-ggsurvplot(fit, pval = TRUE, data = df, risk.table = TRUE, palette = c("blue", "red", "#8DCFEC", "#FF8A8A"), title = names(f[i]), surv.median.line = "hv", risk.table.legend = FALSE, risk.table.col = "strata")
pairwise_survdiff(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ Mut_Status, data=df)



i=29 #CCND1
df <- select(f, colnames(f[i]), PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0, RITUXIMAB_ADMINISTERED)
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'MI'
df$Mut_Status[df[ , colnames(df)[1]] == 1 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'MC'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'Yes'] <- 'WTI'
df$Mut_Status[df[ , colnames(df)[1]] == 0 & df$RITUXIMAB_ADMINISTERED == 'No'] <- 'WTC'
df[df == "N/A"] <- NA
df <- na.omit(df)
df$PROGRESSION_FREE_SURVIVAL_MONTHS <- as.numeric(as.character(df$PROGRESSION_FREE_SURVIVAL_MONTHS))
df$PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0 <- as.numeric(as.character(df$PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0))
fit <- survfit(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ Mut_Status, data=df)
slots[[4]]<-ggsurvplot(fit, pval = TRUE, data = df, risk.table = TRUE, palette = c("blue", "red", "#8DCFEC", "#FF8A8A"), title = names(f[i]), surv.median.line = "hv", risk.table.legend = FALSE, risk.table.col = "strata")
pairwise_survdiff(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ Mut_Status, data=df)


pdf("PFS_selected_genes.pdf", height = 18, width = 20)
arrange_ggsurvplots(slots, print = TRUE,
                    ncol = 2, nrow = 2, risk.table.height = 0.4)
dev.off()



####Figure 1C #####
#Multivariate analysis of mutations with clinical factors which show significant correlation with survival
d1 <- read.csv("Cohort-1-mutation-matrix.csv", header = TRUE, sep = ",", check.names = FALSE)
d2 <- read.csv("Cohort-2-mutation-matrix.csv", header = TRUE, sep = ",", check.names = FALSE)
colnames(d1)[1] <- "Genes"
colnames(d2)[1] <- "Genes"
# Step 2: Merge the two datasets by the "Genes" column
merged_data <- merge(d1, d2, by = "Genes", all = TRUE)

# Step 3: Remove rows with NA values
merged_data <- na.omit(merged_data)

# Step 4: Set the "Genes" column as row names for transposing
rownames(merged_data) <- merged_data$Genes
merged_data$Genes <- NULL  # Remove the "Genes" column as it is now part of row names

# Step 5: Transpose the data
transposed_data <- as.data.frame(t(merged_data))

# Step 6: Print the transposed data
print("Transposed Data:")
print(transposed_data)

# Step 7: Combine CLinical Information
j1<-read.table("WES_130_clinical.tsv", header=TRUE, sep="\t")
j2<-read.table("BC_23_clinical.tsv", header=TRUE, sep="\t")
j3<-read.table("Targeted_137_clinical.tsv", header=TRUE, sep="\t")
com<-rbind(j1, j2, j3)

# Step 8: Add all clinical information into transposed matrix

m<-match(rownames(transposed_data), com[,1])
com<-com[m,]
f<-cbind(transposed_data, com)
write.csv(f, file="inter.csv")

#### Forest Plot ####
f <- read.csv("inter_final.csv", header = TRUE, sep = ",", row.names = 1)
i = 14
# Select relevant columns
df <- select(f, colnames(f[i]), OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0, Transplant, B_Symptom, RITUXIMAB_ADMINISTERED, COHORT_TYPE, AGE_Status)

# Filter rows for "Cohort-1"
df <- subset(df, df$COHORT_TYPE == "Cohort-1")

# Convert the selected column to "Mutated" or "Non-Mutated"
df[, 1] <- ifelse(df[, 1] > 0, "Mutated", "Non-Mutated")

# Replace "N/A" with NA and remove rows with NA values
df[df == "N/A"] <- NA
df <- na.omit(df)

# Convert survival columns to numeric
df$OVERALL_SURVIVAL_MONTHS <- as.numeric(as.character(df$OVERALL_SURVIVAL_MONTHS))
df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0 <- as.numeric(as.character(df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0))

# Cox proportional hazards models
mut_therapy <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:RITUXIMAB_ADMINISTERED, data = df)
mut_age <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:AGE_Status, data = df)
mut_transplant <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:Transplant, data = df)
mut_Bsymptom <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:B_Symptom, data = df)

# Define a function to process and save results for each model
# Update the process_model function to include confidence intervals
process_model <- function(model, filename, gene_name) {
  # Extract coefficients and their summary
  coefficients <- summary(model)$coefficients
  conf_intervals <- confint(model)
  
  # Filter rows where the row names start with "df[, 1]"
  df1_values <- coefficients[grepl("^df\\[, 1\\]", rownames(coefficients)), ]
  df1_conf <- conf_intervals[grepl("^df\\[, 1\\]", rownames(conf_intervals)), ]
  
  # Convert coefficients and confidence intervals to data frames
  df1_values_df <- as.data.frame(df1_values)
  df1_conf_df <- as.data.frame(df1_conf)
  
  # Add confidence interval columns (lower and upper bounds)
  df1_values_df$CI_Lower <- df1_conf_df[, 1]
  df1_values_df$CI_Upper <- df1_conf_df[, 2]
  
  # Add a column for the gene name
  df1_values_df$Gene <- gene_name
  
  # Save to CSV
  write.csv(df1_values_df, filename, row.names = TRUE)
  
  # Return the data frame for verification if needed
  return(df1_values_df)
}

# Process each model and save the results
therapy_results <- process_model(mut_therapy, "mut_therapy_results.csv", colnames(f)[i])
age_results <- process_model(mut_age, "mut_age_results.csv", colnames(f)[i])
transplant_results <- process_model(mut_transplant, "mut_transplant_results.csv", colnames(f)[i])
Bsymptom_results <- process_model(mut_Bsymptom, "mut_Bsymptom_results.csv", colnames(f)[i])

# Combine results into one data frame
atm <- rbind(therapy_results, age_results, transplant_results, Bsymptom_results)

# Update row names to include the gene name
rownames(atm) <- paste(atm$Gene, rownames(atm), sep = ": ")

# View the final data frame with confidence intervals
print(atm)

####
i = 188

# Select relevant columns
df <- select(f, colnames(f[i]), OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0, Transplant, B_Symptom, RITUXIMAB_ADMINISTERED, COHORT_TYPE, AGE_Status)

# Filter rows for "Cohort-1"
df <- subset(df, df$COHORT_TYPE == "Cohort-1")

# Convert the selected column to "Mutated" or "Non-Mutated"
df[, 1] <- ifelse(df[, 1] > 0, "Mutated", "Non-Mutated")

# Replace "N/A" with NA and remove rows with NA values
df[df == "N/A"] <- NA
df <- na.omit(df)

# Convert survival columns to numeric
df$OVERALL_SURVIVAL_MONTHS <- as.numeric(as.character(df$OVERALL_SURVIVAL_MONTHS))
df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0 <- as.numeric(as.character(df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0))

# Cox proportional hazards models
mut_therapy <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:RITUXIMAB_ADMINISTERED, data = df)
mut_age <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:AGE_Status, data = df)
mut_transplant <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:Transplant, data = df)
mut_Bsymptom <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:B_Symptom, data = df)

# Define a function to process and save results for each model
# Update the process_model function to include confidence intervals
process_model <- function(model, filename, gene_name) {
  # Extract coefficients and their summary
  coefficients <- summary(model)$coefficients
  conf_intervals <- confint(model)
  
  # Filter rows where the row names start with "df[, 1]"
  df1_values <- coefficients[grepl("^df\\[, 1\\]", rownames(coefficients)), ]
  df1_conf <- conf_intervals[grepl("^df\\[, 1\\]", rownames(conf_intervals)), ]
  
  # Convert coefficients and confidence intervals to data frames
  df1_values_df <- as.data.frame(df1_values)
  df1_conf_df <- as.data.frame(df1_conf)
  
  # Add confidence interval columns (lower and upper bounds)
  df1_values_df$CI_Lower <- df1_conf_df[, 1]
  df1_values_df$CI_Upper <- df1_conf_df[, 2]
  
  # Add a column for the gene name
  df1_values_df$Gene <- gene_name
  
  # Save to CSV
  write.csv(df1_values_df, filename, row.names = TRUE)
  
  # Return the data frame for verification if needed
  return(df1_values_df)
}

# Process each model and save the results
therapy_results <- process_model(mut_therapy, "mut_therapy_results.csv", colnames(f)[i])
age_results <- process_model(mut_age, "mut_age_results.csv", colnames(f)[i])
transplant_results <- process_model(mut_transplant, "mut_transplant_results.csv", colnames(f)[i])
Bsymptom_results <- process_model(mut_Bsymptom, "mut_Bsymptom_results.csv", colnames(f)[i])

# Combine results into one data frame
tp53 <- rbind(therapy_results, age_results, transplant_results, Bsymptom_results)

# Update row names to include the gene name
rownames(tp53) <- paste(tp53$Gene, rownames(tp53), sep = ": ")

# View the final tp53 data frame with confidence intervals
print(tp53)

#####
i = 29

# Select relevant columns
df <- select(f, colnames(f[i]), OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0, Transplant, B_Symptom, RITUXIMAB_ADMINISTERED, COHORT_TYPE, AGE_Status)

# Filter rows for "Cohort-1"
df <- subset(df, df$COHORT_TYPE == "Cohort-1")

# Convert the selected column to "Mutated" or "Non-Mutated"
df[, 1] <- ifelse(df[, 1] > 0, "Mutated", "Non-Mutated")

# Replace "N/A" with NA and remove rows with NA values
df[df == "N/A"] <- NA
df <- na.omit(df)

# Convert survival columns to numeric
df$OVERALL_SURVIVAL_MONTHS <- as.numeric(as.character(df$OVERALL_SURVIVAL_MONTHS))
df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0 <- as.numeric(as.character(df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0))

# Cox proportional hazards models
mut_therapy <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:RITUXIMAB_ADMINISTERED, data = df)
mut_age <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:AGE_Status, data = df)
mut_transplant <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:Transplant, data = df)
mut_Bsymptom <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:B_Symptom, data = df)

# Define a function to process and save results for each model
# Update the process_model function to include confidence intervals
process_model <- function(model, filename, gene_name) {
  # Extract coefficients and their summary
  coefficients <- summary(model)$coefficients
  conf_intervals <- confint(model)
  
  # Filter rows where the row names start with "df[, 1]"
  df1_values <- coefficients[grepl("^df\\[, 1\\]", rownames(coefficients)), ]
  df1_conf <- conf_intervals[grepl("^df\\[, 1\\]", rownames(conf_intervals)), ]
  
  # Convert coefficients and confidence intervals to data frames
  df1_values_df <- as.data.frame(df1_values)
  df1_conf_df <- as.data.frame(df1_conf)
  
  # Add confidence interval columns (lower and upper bounds)
  df1_values_df$CI_Lower <- df1_conf_df[, 1]
  df1_values_df$CI_Upper <- df1_conf_df[, 2]
  
  # Add a column for the gene name
  df1_values_df$Gene <- gene_name
  
  # Save to CSV
  write.csv(df1_values_df, filename, row.names = TRUE)
  
  # Return the data frame for verification if needed
  return(df1_values_df)
}

# Process each model and save the results
therapy_results <- process_model(mut_therapy, "mut_therapy_results.csv", colnames(f)[i])
age_results <- process_model(mut_age, "mut_age_results.csv", colnames(f)[i])
transplant_results <- process_model(mut_transplant, "mut_transplant_results.csv", colnames(f)[i])
Bsymptom_results <- process_model(mut_Bsymptom, "mut_Bsymptom_results.csv", colnames(f)[i])

# Combine results into one data frame
ccnd1 <- rbind(therapy_results, age_results, transplant_results, Bsymptom_results)

# Update row names to include the gene name
rownames(ccnd1) <- paste(ccnd1$Gene, rownames(ccnd1), sep = ": ")

# View the final ccnd1 data frame with confidence intervals
print(ccnd1)



####
i = 94
# Select relevant columns
df <- select(f, colnames(f[i]), OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0, Transplant, B_Symptom, RITUXIMAB_ADMINISTERED, COHORT_TYPE, AGE_Status)

# Filter rows for "Cohort-1"
df <- subset(df, df$COHORT_TYPE == "Cohort-1")

# Convert the selected column to "Mutated" or "Non-Mutated"
df[, 1] <- ifelse(df[, 1] > 0, "Mutated", "Non-Mutated")

# Replace "N/A" with NA and remove rows with NA values
df[df == "N/A"] <- NA
df <- na.omit(df)

# Convert survival columns to numeric
df$OVERALL_SURVIVAL_MONTHS <- as.numeric(as.character(df$OVERALL_SURVIVAL_MONTHS))
df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0 <- as.numeric(as.character(df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0))

# Cox proportional hazards models
mut_therapy <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:RITUXIMAB_ADMINISTERED, data = df)
mut_age <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:AGE_Status, data = df)
mut_transplant <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:Transplant, data = df)
mut_Bsymptom <- coxph(Surv(df$OVERALL_SURVIVAL_MONTHS, df$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ df[, 1]:B_Symptom, data = df)

# Define a function to process and save results for each model
# Update the process_model function to include confidence intervals
process_model <- function(model, filename, gene_name) {
  # Extract coefficients and their summary
  coefficients <- summary(model)$coefficients
  conf_intervals <- confint(model)
  
  # Filter rows where the row names start with "df[, 1]"
  df1_values <- coefficients[grepl("^df\\[, 1\\]", rownames(coefficients)), ]
  df1_conf <- conf_intervals[grepl("^df\\[, 1\\]", rownames(conf_intervals)), ]
  
  # Convert coefficients and confidence intervals to data frames
  df1_values_df <- as.data.frame(df1_values)
  df1_conf_df <- as.data.frame(df1_conf)
  
  # Add confidence interval columns (lower and upper bounds)
  df1_values_df$CI_Lower <- df1_conf_df[, 1]
  df1_values_df$CI_Upper <- df1_conf_df[, 2]
  
  # Add a column for the gene name
  df1_values_df$Gene <- gene_name
  
  # Save to CSV
  write.csv(df1_values_df, filename, row.names = TRUE)
  
  # Return the data frame for verification if needed
  return(df1_values_df)
}

# Process each model and save the results
therapy_results <- process_model(mut_therapy, "mut_therapy_results.csv", colnames(f)[i])
age_results <- process_model(mut_age, "mut_age_results.csv", colnames(f)[i])
transplant_results <- process_model(mut_transplant, "mut_transplant_results.csv", colnames(f)[i])
Bsymptom_results <- process_model(mut_Bsymptom, "mut_Bsymptom_results.csv", colnames(f)[i])

# Combine results into one data frame
kmt2c <- rbind(therapy_results, age_results, transplant_results, Bsymptom_results)

# Update row names to include the gene name
rownames(kmt2c) <- paste(kmt2c$Gene, rownames(kmt2c), sep = ": ")

# View the final kmt2c data frame with confidence intervals
print(kmt2c)

final<-rbind(atm, tp53, kmt2c, ccnd1)
mutated_rows <- final[grepl("Mutated", rownames(final)) & !grepl("Non-Mutated", rownames(final)), ]
mutated_rows <- na.omit(mutated_rows)
rownames(mutated_rows) <- gsub("df\\[, 1\\]", "", rownames(mutated_rows))
rownames(mutated_rows) <- gsub(": Mutated:", " MUT:", rownames(mutated_rows))
rownames(mutated_rows) <- gsub("B_Symptom0", "B Symptom Yes", rownames(mutated_rows))
rownames(mutated_rows) <- gsub("B_Symptom1", "B Symptom No", rownames(mutated_rows))
rownames(mutated_rows) <- gsub("Transplant0", "Transplant Yes", rownames(mutated_rows))
rownames(mutated_rows) <- gsub("Transplant1", "Transplant No", rownames(mutated_rows))
data<-mutated_rows


write.csv(data, file="forest_plot.csv")


# Add significance levels as stars
data$Significance <- cut(
  data$`Pr(>|z|)`,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns"),
  right = FALSE
)

#Remove all insignificant interactions 
data<-subset(data, data$Significance != "ns")
data<-data[, c(2,6,7,9)]
data$CI_Lower <- exp(data$CI_Lower)
data$CI_Upper <- exp(data$CI_Upper)
colnames(data)<-c("HR", "CI_Lower", "CI_Upper", "Significance")
pdf("Figure_1C.pdf", height=10, width=15)
ggplot(data, aes(x = HR, y = rownames(data))) +
  geom_point(aes(color = Significance), size = 3) +  # Points for hazard ratios
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +  # Horizontal error bars for CI
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +  # Reference line at HR = 1
  scale_x_log10() +  # Log scale for hazard ratio
  labs(
    title = "Forest Plot Showing Interaction of Covariates",
    x = "Hazard Ratio (HR)",
    y = "Intercation Terms"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  ) + theme_classic()

dev.off()



