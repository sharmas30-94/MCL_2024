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
data <- read.table("/Users/sharmas30/Desktop/code for MCL/Figure3/DLBCL_Nat_Med_April_2018-master/consensus_clustering/input_data/MCL_Significant_events.txt", header=TRUE, sep="\t", check.names=FALSE)

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

iteration_dirs <- list.files("/Users/sharmas30/Desktop/code for MCL/Figure3/DLBCL_Nat_Med_April_2018-master/consensus_clustering/input_data")

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

selected_set <- "input_set_94"
memb<- "membership.2"
setwd(paste0("/Users/sharmas30/Desktop/code for MCL/Figure3/DLBCL_Nat_Med_April_2018-master/consensus_clustering/input_data/", selected_set))
print(getwd())
 
data1<-read.table("MCL.membership.txt", header=TRUE, sep="\t", row.names=1)
data1 <- data1[order(data1[, memb]), ]
data2 <- read.table("/Users/sharmas30/Desktop/code for MCL/Figure3/DLBCL_Nat_Med_April_2018-master/consensus_clustering/input_data/input_set_94/set_94.txt", header=TRUE, sep="\t", row.names=1)
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
  
  # Extract text content from the PDF
  pdf_text_content <- pdf_text(pdf_file)
  cleaned_text <- gsub("\\s+", " ", pdf_text_content)
  
  # Extract column and row names from the cleaned text
  cols <- colnames(subset_data) # Extract columns specific to this subset
  rows <- rownames(subset_data) # Extract rows specific to this subset
  
  # Create regex pattern with word boundaries for accurate matching
  cols_pattern <- paste0("\\b(", paste(cols, collapse = "|"), ")\\b")
  rows_pattern <- paste0("\\b(", paste(rows, collapse = "|"), ")\\b")
  
  # Extract matched column names and row names
  extracted_cols <- unique(unlist(str_extract_all(cleaned_text, cols_pattern)))
  extracted_rows <- unique(unlist(str_extract_all(cleaned_text, rows_pattern)))
  
  # Append to final lists
  final_cols <- unique(c(final_cols, extracted_cols))
  final_rows <- unique(c(final_rows, extracted_rows))
}

# Print or store final extracted column and row names
print(final_cols)
print(final_rows)


# Read the data again
data3 <- read.csv("check.csv", header = TRUE, sep = ",", row.names = 1, check.names = FALSE)


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
  column_order = sample_ids, 
  row_order = genomic_order  
)
dev.off()






# Load required libraries
library(pdftools)
library(stringr)
library(dplyr)

# Define the PDF file path
pdf_file <- "Survival_Analysis_Report.pdf"

# Extract text from the PDF
pdf_text <- pdf_text(pdf_file)

# Combine text into a single string
pdf_text_combined <- paste(pdf_text, collapse = " ")

# Extract p-values using regex
p_values <- str_extract_all(pdf_text_combined, "p\\s*=\\s*[0-9\\.]+")[[1]]

# Clean extracted p-values
p_values <- as.numeric(str_extract(p_values, "[0-9\\.]+"))

# Create a data frame with the extracted p-values
df <- data.frame(PDF_Name = pdf_file, t(p_values))
colnames(df)[-1] <- paste0("membership", ifelse(seq_along(p_values) == 1, "", paste0(".", seq_along(p_values) - 1)))

# Print the data frame
print(df)

# Save as CSV (optional)
write.csv(df, "extracted_p_values.csv", row.names = FALSE)




# Load necessary libraries
library(dplyr)
library(tidyr)

# Load your data (update the path if needed)
setwd("/Users/sharmas30/Downloads/Fig4")
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

