#Consensus clustering
library("dplyr")
library("survival")
library("survminer")
library("ggplot2")
library("gridExtra")
library("ComplexHeatmap")
library('pdftools')
library("stringr")
library("impute")

# Read the data file
data <- read.table("~/Figure3/Data/DLBCL_Nat_Med_April_2018-master/consensus_clustering/input_data/MCL_Significant_events.txt", header=TRUE, sep="\t", check.names=FALSE)

# Define the output folder
input_folder <- "input_data"
dir.create(input_folder, showWarnings = FALSE)

# Number of rows to select in each iteration
rows_to_select <- 35
iterations <- 300

# List of required rows
required_rows <- c("ATM", "Deletion_11q22.3", "TP53", "Deletion_TP53", "Amplification-8q23.2-8q24.21", "Deletion_9p21.1", "Amplification-3q22.1-3q22.3")

# Loop through each iteration
for (i in 1:iterations) {
  
  # Ensure required rows are included
  required_subset <- data[data$gene %in% required_rows, ]  # Adjust Column_Name to match actual column
  
  # Sample additional rows excluding required rows
  remaining_rows <- data[!data$gene %in% required_rows, ]
  sampled_rows <- remaining_rows %>% sample_n(rows_to_select - nrow(required_subset))
  
  # Combine required and sampled rows
  final_sample <- rbind(required_subset, sampled_rows)
  
  # Create a subfolder for this iteration
  iteration_folder <- file.path(input_folder, paste0("input_set_", i))
  dir.create(iteration_folder, showWarnings = FALSE)
  
  # Define the filename inside the iteration folder
  input_filename <- file.path(iteration_folder, paste0("set_", i, ".txt"))
  
  # Save the sampled rows to the iteration folder
  write.table(final_sample, file = input_filename, row.names = FALSE, sep="\t", quote=FALSE)
}

iteration_dirs <- list.files("~/Figure3/Data/DLBCL_Nat_Med_April_2018-master/consensus_clustering/input_data")

# Sort directories numerically
iteration_dirs <- iteration_dirs[order(as.numeric(gsub("\\D", "", basename(iteration_dirs))))]

# Print sorted directories
print(iteration_dirs)

setwd("/Users/sharmas30/Desktop/code for MCL/Figure3/DLBCL_Nat_Med_April_2018-master/consensus_clustering/")

for (iteration_dir in 1:length(iteration_dirs)) {
  
  # Run Topgenes_v1.R
  source('src/GDAC_TopgenesforCluster/Topgenes_v1.R')
  result <- main("-s./src/GDAC_TopgenesforCluster/", paste0("-m", "input_data/", iteration_dirs[iteration_dir], "/set_", paste0(iteration_dir), ".txt"), "-uALL","-oMCL")
  source("src/GDAC_NmfConsensusClustering/GDAC_CNMF.R", echo = TRUE)
  result <- main("-ssrc/GDAC_NmfConsensusClustering/","-mMCL.expclu.gct" ,"-oMCL",'-u4','-v10') 
  source("src/GDAC_selectBestcluster/select_best_cluster_chc.R")
  result <- main("-uMCL.expclu.gct","-mPearson","-cMCL.cophenetic.coefficient.txt","-wMCL.membership.txt","-voutput_dir/MCL", paste0("-p", "input_data/", iteration_dirs[iteration_dir], "/set_", paste0(iteration_dir), ".txt"))
  
  
  ## Run Kaplan Meyer Curve on each set ####
  library("survival")
  library("survminer")
  data1 <- read.table("MCL.membership.txt", header = TRUE, sep = "\t", row.names = 1)
  # List of membership columns to process (replace these with your actual column names)
  membership_columns <- c("membership", "membership.1", "membership.2", "membership.3", "membership.4", "membership.5", "membership.6")
  
 
  # Loop through each membership column
  data2 <- read.table("SURVIVAL.csv", header=TRUE, sep=",")
  
  # Match rownames from data1 with the first column of data2
  m <- match(rownames(data1), data2[, 1])
  dat_original <- data2[m, ]
  
  # Open a multi-page PDF
  pdf("Survival_Analysis_Report.pdf", width=8, height=12)
  
  # Loop over all columns of data1 (clusters)
  for (k in membership_columns) {
    dat <- dat_original
    dat$Cluster <- data1[, k]  # Assign cluster from the current column
    
    # Select necessary columns
    dat <- dat[, c("SAMPLE_ID", "OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0", 
                   "OVERALL_SURVIVAL_MONTHS", "Cluster", 
                   "RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0")]
    
    # Subset for Rituximab-treated patients
    dat <- subset(dat, dat$RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0 == 1)
    
    # Convert columns to numeric
    dat$OVERALL_SURVIVAL_MONTHS <- as.numeric(as.character(dat$OVERALL_SURVIVAL_MONTHS))
    dat$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0 <- as.numeric(as.character(dat$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0))
    
    # Remove NA values
    dat <- na.omit(dat)
    
    # Fit survival model
    fit <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ Cluster, data=dat)
    
    # Generate survival plot
    surv_plot <- ggsurvplot(fit, pval=TRUE, data=dat, risk.table=TRUE, title=paste("Survival Curve -", k))
    
    # Print plot to PDF
    print(surv_plot)
    
  }
  
  # Close PDF
  dev.off()
  
  #Move files to directories
  source_dir <- "/Users/sharmas30/Desktop/code for MCL/Figure3/DLBCL_Nat_Med_April_2018-master/consensus_clustering"
  dest_dir <- paste0("/Users/sharmas30/Desktop/code for MCL/Figure3/DLBCL_Nat_Med_April_2018-master/consensus_clustering/input_data/", iteration_dirs[[iteration_dir]])
  file_patterns <- c("*.png", "*.txt", "*.gct", "*.pdf")
  files_to_move <- list.files(source_dir, pattern = "\\.(png|txt|gct|pdf)$", full.names = TRUE)
  file.rename(files_to_move, file.path(dest_dir, basename(files_to_move)))
  source_dir1 <- "/Users/sharmas30/Desktop/code for MCL/Figure3/DLBCL_Nat_Med_April_2018-master/consensus_clustering/output_dir"
  file_patterns <- c("*.png", "*.txt")
  dest_dir1 <- paste0("/Users/sharmas30/Desktop/code for MCL/Figure3/DLBCL_Nat_Med_April_2018-master/consensus_clustering/input_data/", iteration_dirs[[iteration_dir]])
  files_to_move1 <- list.files(source_dir1, pattern = "\\.(png|txt)$", full.names = TRUE)
  file.rename(files_to_move1, file.path(dest_dir1, basename(files_to_move1)))
  }




#Selected Set

selected_set <- "input_set_86"
memb<- "membership.2"
setwd(paste0("/Users/sharmas30/Desktop/code for MCL/Figure3/DLBCL_Nat_Med_April_2018-master/consensus_clustering/input_data/", selected_set))
print(getwd())
 
data1<-read.table("MCL.membership.txt", header=TRUE, sep="\t", row.names=1)
data1 <- data1[order(data1[, memb]), ]
data2 <- read.table("/Users/sharmas30/Desktop/code for MCL/Figure3/DLBCL_Nat_Med_April_2018-master/consensus_clustering/input_data/input_set_86/set_86.txt", header=TRUE, sep="\t", row.names=1)
rownames(data1)<-gsub("-", ".", rownames(data1))
m<-match(rownames(data1), colnames(data2))
data3<-data2[,m]


#create an output directory. Folder created by membership column
data5<-data.frame()
################################################################
class<-unique(data1[, memb])
################################################################
for(i in 1:length(class)){
 l<-c()
 m<-c()
 n<-c()
 o<-c()
 names<-c()
 ################################################################
 a<-which(data1[,memb] == class[i])
 ################################################################ ?
 start<-a[1]
 end<-a[length(a)]
 for(j in 1:nrow(data3)){
    row1<-data3[j, start:end]
    row2<-data3[j, -c(start:end)]
    altClust1<- row1 %>% select_if(function(col) col > 0)
    notaltClust1<- row1 %>% select_if(function(col) col == 0)
    altClustother<-row2 %>% select_if(function(col) col > 0)
    notaltClustother<- row2 %>% select_if(function(col) col == 0)
    l<-c(l, ncol(altClust1))
    m<-c(m, ncol(notaltClust1))
    n<-c(n, ncol(altClustother))
    o<-c(o, ncol(notaltClustother))
    names<-c(names, rownames(data3)[j])

}
data4<-data.frame(names, l,m, n, o)
write.csv(data4, file=paste0(i, "_frequency.csv"), row.names=FALSE)
rnames<-c()
pvalue<-c()
odds.ratio<-c()
for(k in 1:nrow(data4)){
   cont_table=matrix(c(data4[k,2], data4[k,3], data4[k,4], data4[k,5]), nrow=2)
   test<-fisher.test(cont_table)
   pvalue<-c(pvalue, test$p.value)
   odds.ratio<-c(odds.ratio, test$estimate)
   rnames<-c(rnames, data4[k,1])
}
 result<-data.frame(rnames, pvalue, odds.ratio)
 result<-subset(result, result$pvalue <= 0.1)
 result<-result[order(result$pvalue),]
 result$cluster=c(rep(i, times = nrow(result)))
 write.csv(result, file=paste0(i, "_significant_events.csv"), row.names=FALSE)
 data5<-rbind(data5, result)
 write.csv(data5, file="All.significant.events.csv", row.names=FALSE)
 
}

#retain only one of the duplicate events with highest odds.ratio across different clusters
data<-read.csv("All.significant.events.csv", header=TRUE, sep=",")
data_ordered<- data[order(data$odds.ratio, decreasing = TRUE), ]
data_highest<-data_ordered[!duplicated(data_ordered$rnames),]
data_highest
data_highest<-data_highest[order(data_highest$cluster),]
write.csv(data_highest, file="unique_across_groups.csv", row.names=FALSE)

#Arrange samples and events for heatmap
match_row<-match(data_highest[,1], rownames(data3))
data6<-data3[match_row,]
write.csv(data6, file="Heatmap.csv", row.names=TRUE)


# Load the dataset
df <- read.csv("Heatmap.csv", header = TRUE, sep = ",", row.names = 1, check.names = FALSE)

# --- Handle Amplification Rows ---
amplification_rows <- grepl("Amplification", rownames(df), ignore.case = TRUE)
df[amplification_rows, ] <- lapply(df[amplification_rows, ], function(x) {
  x[x == 2] <- "Broad_Amp;"
  x[x == 1] <- "Low_Amp;"
  x[x == 0] <- ""
  return(x)
})

# --- Handle Deletion Rows ---
deletion_rows <- grepl("Deletion", rownames(df), ignore.case = TRUE)
df[deletion_rows, ] <- lapply(df[deletion_rows, ], function(x) {
  x[x == 1] <- "Loss;"
  x[x == 0] <- "";
  x[x == 2] <- "Gain;"
  return(x)
})

rows_to_modify <- which(!(grepl("Amplification|Deletion", rownames(df), ignore.case = TRUE)))

df[rows_to_modify, ] <- lapply(df[rows_to_modify, ], function(x) {
  x[x == 0] <- ""
  x[x > 0] <- "MUT;"
  return(x)
})

# Save the modified dataframe
write.csv(df, file = "check.csv")


alter_fun_list = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "white", col = NA))
  },
  Low_Amp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "salmon1", col = NA))
  },
  Broad_Amp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "red", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "black", col = NA))
  },
  Loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  Gain = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "maroon", col = NA))
  }
)

col = c(
  "Low_Amp" = "salmon1",
  "Broad_Amp" = "red",
  "MUT" = "black",
  "Loss" = "blue",
  "Gain" = "maroon"
)

# Read the data from CSV files
data3 <- read.csv("check.csv", header = TRUE, sep = ",", check.names = FALSE, row.names = 1)
data4 <- read.csv("unique_across_groups.csv", header = TRUE, sep = ",", check.names = FALSE)

# Remove leading 'X' and replace '-' with '.' in column names
colnames(data3) <- gsub("^X", "", colnames(data3))
colnames(data3) <- gsub("-", ".", colnames(data3))

# Initialize vectors to store final column and row orders
final_cols <- c()
final_rows <- c()

# Loop through clusters
data1$clust <- data1[, memb]
k <- unique(data1$clust)

for (i in 1:length(k)) {
  # Subset data based on cluster
  subset_data <- data3[, data1$clust == k[i]]
  ids_to_extract <- which(data4$cluster == k[i])
  subset_data <- subset_data[ids_to_extract, ]
  
  # Create oncoPrint plot
  h <- oncoPrint(
    subset_data,
    get_type = function(x) strsplit(x, ";")[[1]],
    alter_fun = alter_fun_list,
    col = col,
    remove_empty_columns = FALSE,
    row_names_gp = gpar(fontsize = 5),
    alter_fun_is_vectorized = FALSE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 4)
  )
  
  # Save oncoPrint plot as a PDF
  pdf_file <- paste0(k[i], ".pdf")
  pdf(pdf_file, height = 6, width = 6)
  print(h) # Print the plot
  dev.off()
  
}

# Read the data again
data3 <- read.csv("arranged.csv", header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
pdf("new10.pdf", height = 15, width = 20)
oncoPrint(
  data3,
  get_type = function(x) strsplit(x, ";")[[1]],
  alter_fun = alter_fun_list,
  col = col,
  remove_empty_columns = FALSE,
  remove_empty_rows = FALSE,
  row_names_gp = gpar(fontsize = 15),
  alter_fun_is_vectorized = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8),
  column_order = 1:ncol(data3), 
  row_order = 1:nrow(data3) 
)
dev.off()

#### Clinical cluster analysis ####

library(dplyr)
library(tidyr)

# Load your data (update the path if needed)
setwd("")
df <- read.csv("Cohort-1.csv")

# List of binary clinical variables to test
clinical_vars <- c(
  "SOX11_IHC_GREATER_THAN_EQUAL_10_1_LESS_THAN_10_0",
  "MTP53_GREATER_THAN_EQUAL_90_1_ELSE_0",
  "Ki67_COUNTS_GREATER_THAN_30_1_ELSE_0",
  "SUBTYPES_CLASSICAL_100_0_B.P_1",
  "GROWTH_PATTERN_OTHER_0_DIFFUSE_100_1",
  "GENDER_M_1_F_0",
  "RACE_CAUCASIAN_1_OTHER_0",
  "NODAL_1_EXTRANODAL_2_NODAL_WITH_EXTRANODAL_3",
  "PERFORMANCE_ABMULATORY_0_ELSE_1",
  "AAS_III_IV_1_II_I_0",
  "WHETHER_TRANSPLANT_YES_0_NO_1",
  "AGE_GREATER_THAN_EQUAL_60_1_LESS_THAN_60_0",
  "RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0"
)

# Initialize results list
results <- list()

# Loop through each feature
for (var in clinical_vars) {
  sub_df <- df %>%
    select(Cluster, all_of(var)) %>%
    filter(!is.na(Cluster), !is.na(.data[[var]]))
  
  if (nrow(sub_df) == 0) next
  
  # Create contingency table
  tbl <- table(sub_df$Cluster, sub_df[[var]])
  
  # Run chi-square test
  test <- chisq.test(tbl)
  p_value <- test$p.value
  
  # Compute proportions of '1' per cluster
  proportions <- prop.table(tbl, margin = 1)
  
  # Identify enrichment and depletion
  if ("1" %in% colnames(proportions)) {
    enriched_cluster <- names(which.max(proportions[, "1"]))
    depleted_cluster <- names(which.min(proportions[, "1"]))
    enriched_category <- 1
  } else if ("0" %in% colnames(proportions)) {
    enriched_cluster <- names(which.max(proportions[, "0"]))
    depleted_cluster <- names(which.min(proportions[, "0"]))
    enriched_category <- 0
  } else {
    enriched_cluster <- NA
    depleted_cluster <- NA
    enriched_category <- NA
  }
  
  results[[var]] <- data.frame(
    Feature = var,
    P_value = p_value,
    Significant = p_value < 0.05,
    Enriched_Category = enriched_category,
    Enrichment = if (!is.na(enriched_cluster)) paste0("Higher in Cluster ", enriched_cluster) else NA,
    Depletion = if (!is.na(depleted_cluster)) paste0("Lower in Cluster ", depleted_cluster) else NA
  )
}

# Combine and print results
enrichment_results <- do.call(rbind, results)
print(enrichment_results)

# Optionally, write to CSV
write.csv(enrichment_results, "cluster_enrichment_results.csv", row.names = FALSE)
 


##Arange Clinical data by the order of the clusters
 column_order <- c(
   "GROUP6_1", "UNMC_44", "47_3", "MCL_28", "GROUP7_15", "MCL_38", "MCL_98", "UNMC_28", "MCL_128",
   "MCL_35", "UNMC_23", "MCL_71", "GROUP5_7", "MCL_24", "UNMC_37", "IA_02", "MCL_109", "SC01_10",
   "M47_025", "IA4203", "UNMC_79_C1", "UNMC_47", "GROUP5_10", "UNMC_81", "UNMC_59", "UNMC_71_A5",
   "GROUP5_9", "GROUP4_11", "WCM032", "GROUP4_8", "UNMC_8", "GROUP7_5", "GROUP5_3", "HUMC_9",
   "1_33", "GROUP7_9", "1_87", "IA2084", "GROUP7_20", "HUMC4", "S92_5", "HUMC_7", "IA4227",
   "MCL_40", "UNMC_48", "GROUP7_22", "SPECS_3164", "IA4289", "GROUP4_3", "S03_8", "MCL_146", "HUMC_8",
   "UNMC_18", "HUMC_6", "S03_3", "HUMC3", "GROUP7_13", "1_90", "GROUP4_9", "GROUP7_1", "1_16",
   "WCM062", "UNMC_5", "47_027", "IA1857", "GROUP4_2", "GROUP5_6", "GROUP7_17", "22_51", "MCL_37",
   "IA2268", "GROUP4_1", "GROUP4_4", "UNMC_20", "GROUP6_3", "MCL_39", "GROUP7_6", "1_97",
   "GROUP7_3", "1_31", "GROUP4_5", "1_39", "1_64", "1_61", "GROUP4_6", "HUMC_12",
   "1_88", "GROUP5_2", "1_6", "1_86", "1_85", "SC08_2", "IA2745", "HUMC5", "GROUP5_5",
   "IA834", "GROUP5_8", "MCL_126", "MCL_14", "IA04", "IA2852", "GROUP5_4", "MCL_61", "HUMC1",
   "IA2334", "MCL_122", "UNMC_7", "UNMC_68_A2", "GROUP5_1", "S87_4", "SPECS_3163", "SPECS_3155",
   "GROUP7_11", "SC02_11", "S06_166", "S89_9", "HUMC_10", "WCM041", "MCL_42", "UNMC_55", "GROUP7_14",
   "47_022", "IA4885", "S08_1", "1_84", "SC86_6", "IA68", "47_023", "1_51", "GROUP7_2",
   "47_024", "HUMC_15", "BC04", "47_029", "1_24", "UNMC_69", "SPECS_3157", "GROUP7_19",
   "MCL_138", "MCL_33", "SPECS_3158", "GROUP6_4", "GROUP6_2", "IA1147", "GROUP4_7", "S94_7",
   "WCM045", "GROUP4_10", "UNMC_62", "SPECS_1917", "MCL_75", "MCL_45", "MCL_25"
 )
 
 
 pdf("Figure_3_clinical_annotation.pdf", height=20, width=15)
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
   sampleOrder = column_order 
 )
 dev.off()

 
 
 ## Figure 3B
 ##Left Panel
 data<-read.csv("~/Figure3/Data/analysis_nmf.csv", header=TRUE, sep=",")
 data<-data[, c(1, 2, 40, 41, 59)]
 data<-na.omit(data)
 data<-subset(data, data$RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0 == "1")
 surv_object <- Surv(time = data$OVERALL_SURVIVAL_MONTHS, event = data$OVERALL_SURVIVAL_STATUS)
 km_fit <- survfit(surv_object ~ NMF.membership.2, data = data)
 pdf("surv_left.pdf", height = 8, width=10)
 ggsurvplot(km_fit, data = data,
            title = "Kaplan-Meier Survival Curve",
            xlab = "Survival Time (Months)",
            ylab = "Survival Probability",
            surv.median.line = "hv",   # Add median survival line
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_minimal(), pval = TRUE)  # 
 dev.off()
 
##Right Panel
 data<-read.csv("~/Figure3/Data/analysis_nmf.csv", header=TRUE, sep=",")
 data<-data[, c(1, 2, 40, 41, 59)]
 data<-na.omit(data)
 data<-subset(data, data$RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0 == "1")
 surv_object <- Surv(time = data$OVERALL_SURVIVAL_MONTHS, event = data$OVERALL_SURVIVAL_STATUS)
 km_fit <- survfit(surv_object ~ NMF.membership.2, data = data)
 pdf("surv_left.pdf", height = 8, width=10)
 ggsurvplot(km_fit, data = data,
            title = "Kaplan-Meier Survival Curve",
            xlab = "Survival Time (Months)",
            ylab = "Survival Probability",
            surv.median.line = "hv",   # Add median survival line
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_minimal(), pval = TRUE)  # 
 dev.off()
 
 # Load and clean the data
 data <- read.csv("~/Figure3/Data/hj.csv", header = TRUE, sep = ",")
 data <- data[, c(1, 2, 3, 32, 33, 51)]
 data <- na.omit(data)
 
 # Subset only Rituximab-treated patients
 data <- subset(data, RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0 == 1)
 
 # Ensure survival status is numeric (0 = alive, 1 = dead)
 data$OVERALL_SURVIVAL_STATUS <- as.numeric(as.character(data$OVERALL_SURVIVAL_STATUS))
 
 # Ensure cluster is a factor
 data$Merged.cluster <- as.factor(data$Merged.cluster)
 
 # Build the survival object and fit
 surv_object <- Surv(time = data$OVERALL_SURVIVAL_MONTHS, event = data$OVERALL_SURVIVAL_STATUS)
 km_fit <- survfit(surv_object ~ Merged.cluster, data = data)
 
 # Plot Kaplan-Meier curve
 ggsurvplot(km_fit, data = data,
            title = "Kaplan-Meier Survival Curve",
            xlab = "Survival Time (Months)",
            ylab = "Survival Probability",
            surv.median.line = "hv",
            risk.table = TRUE,
            pval = TRUE,
            palette = c("black", "grey", "red"),   # Correct spelling
            ggtheme = theme_minimal())

 
 ## Figure 3C analysis performed using GSEA
 data<-read.table("~/Figure3/Data/RNA_seq_Heatmap.txt", header=TRUE, sep="\t")
 if (!require("pheatmap")) install.packages("pheatmap")
 library(pheatmap)
 
 # Set rownames as pathway names
 rownames(data) <- data$NAME
 
 # Select only numeric expression columns
 mat <- as.matrix(data[, c("Average.poor", "Average_Int", "Average_Good")])
 
 # Generate the heatmap
 pheatmap(mat,
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          scale = "row",  # normalize each row (optional)
          color = colorRampPalette(c("navy", "white", "firebrick3"))(300),
          main = "Hallmark Signature Activity Across Prognostic Groups")
 
 ## Figure 3D
 ### Cibersort Analysis in three groups
 if (!require("FSA")) install.packages("FSA", dependencies=TRUE)
 install.packages("rcompanion")
 
 library(FSA)
 library("rcompanion")
 
 
 
 data<-read.csv("~/Figure3/Data/nmf_analysis_corrected.csv", header=TRUE, sep=",")
 data$Merged.cluster <- NA
 
 data$Merged.cluster[data$NMF.membership.2 %in% c('2', '3', '4')] <- "Moderate OS"
 data$Merged.cluster[data$NMF.membership.2 == '5'] <- "Poor OS"
 data$Merged.cluster[data$NMF.membership.2 %in% c('1', '6')] <- "Good OS"
 
 # Remove rows where Merged.cluster is NA
 data <- na.omit(data)
 
 # Filter dataset for RNA.Seq.Status = TRUE
 data <- data[data$RNA.Seq.Status == TRUE, ]
 
 # Convert Merged.cluster to a factor (categorical variable)
 data$Merged.cluster <- as.factor(data$Merged.cluster)
 
 # Define selected immune-related columns
 selected_columns <- c("T.cells.follicular.helper", "Macrophages.M0")
 
 # Create an empty list to store plots
 plots <- list()
 
 # Loop through selected immune-related columns to generate violin-box plots
 for (i in seq_along(selected_columns)) {
   col <- selected_columns[i]
   
   # Perform Dunn's test
   dunn_result <- tryCatch({
     dunnTest(data[[col]] ~ data$Merged.cluster, method = "bonferroni")
   }, error = function(e) NULL)
   
   if (!is.null(dunn_result)) {
     dunn_pvals <- as.data.frame(dunn_result$res)
     colnames(dunn_pvals) <- c("group1", "group2", "Statistic", "p_value")
     
     # Format p-values for labeling
     dunn_pvals$label <- ifelse(dunn_pvals$p_value < 0.001, "***", 
                                ifelse(dunn_pvals$p_value < 0.01, "**", 
                                       ifelse(dunn_pvals$p_value < 0.05, "*", "ns")))
   } else {
     dunn_pvals <- NULL
   }
   
   # Generate the plot
   p <- ggplot(data, aes(x = Merged.cluster, y = .data[[col]], fill = Merged.cluster)) +
     geom_violin(alpha = 0.6, trim = FALSE) +
     geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA) +
     stat_compare_means(method = "kruskal.test", label = "p.signif") +
     theme_minimal() +
     labs(title = paste(col, "Distribution Across Merged Clusters"), x = "Merged Cluster", y = col) +
     theme(legend.position = "none")
   
   if (!is.null(dunn_pvals) && nrow(dunn_pvals) > 0) {
     p <- p + stat_pvalue_manual(dunn_pvals, label = "label", tip.length = 0.01, y.position = max(data[[col]], na.rm = TRUE) * 1.1)
   }
   
   plots[[i]] <- p
 }
 
 # Save plots as a PDF with a layout
 pdf("selected_immune_violin_boxplot.pdf", height = length(selected_columns) * 3, width = 10)
 do.call(grid.arrange, c(plots, ncol = 1))
 dev.off()
 
 
 
