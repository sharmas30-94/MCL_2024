setwd('/Users/ssharma/Downloads/DLBCL_Nat_Med_April_2018-master 2/consensus_clustering')
data<-read.table("Targeted_149_nmf_new.txt", header=TRUE, sep="\t", check.names=FALSE)
colnames(data)<-gsub(".ATM.MUT", "", colnames(data))
colnames(data)<-gsub(".TP53.MUT", "", colnames(data))
colnames(data)<-gsub(".TP53.CN", "", colnames(data))
colnames(data)<-gsub(".ATM.CN", "", colnames(data))
colnames(data)<-gsub(".WT", "", colnames(data))
colnames(data)<-gsub(".ATM.TP53.BOTH", "", colnames(data))
colnames(data)<-gsub(".ATM.DL", "", colnames(data))
colnames(data)<-gsub(".TP53.DL", "", colnames(data))

library(dplyr)
library(purrr)
library(pdftools)
library(stringr)

output_folder <- "input_data"
dir.create(output_folder, showWarnings = FALSE)

# Number of rows to select in each iteration
rows_to_select <- 60

# Number of iterations
iterations <- 200

# Loop through each iteration
for (i in 1:iterations) {
  # Randomly select 5 rows from the dataset
  sampled_rows <- data %>%
    sample_n(rows_to_select)
  
  # Create a filename for the current iteration
  output_filename <- file.path(output_folder, paste0("iteration_", i, ".txt"))
  
  # Save the sampled rows to a CSV file
  write.table(sampled_rows, file = output_filename, row.names = FALSE, sep="\t", quote=FALSE)
  
}
############
##################
##################################
###############################################
setwd('/Users/ssharma/Downloads/DLBCL_Nat_Med_April_2018-master 2/consensus_clustering')
iteration_files <- list.files('/Users/ssharma/Downloads/DLBCL_Nat_Med_April_2018-master 2/consensus_clustering/input_data/')
#iteration_file <- which(iteration_files == "iteration_69.txt")


# Loop through each iteration file
for (iteration_file in iteration_files) {
  
  # Run Topgenes_v1.R (ensure the correct paths are used within the script)
  source('src/GDAC_TopgenesforCluster/Topgenes_v1.R')
  result <- main("-s./src/GDAC_TopgenesforCluster/", paste0("-minput_data/", iteration_file), "-uALL", "-oDLBCL")
  
  # Run GDAC_CNMF.R (ensure the correct paths are used within the script)
  source("src/GDAC_NmfConsensusClustering/GDAC_CNMF.R", echo = TRUE)
  result <- main("-ssrc/GDAC_NmfConsensusClustering/", "-mDLBCL.expclu.gct", "-oDLBCL", '-u4', '-v10')
  
  # Run select_best_cluster_chc.R (ensure the correct paths are used within the script)
  source("src/GDAC_selectBestcluster/select_best_cluster_chc.R")
  
  # Load necessary libraries
  library("survival")
  library("survminer")
  
  # Read the membership data
  data1 <- read.table("DLBCL.membership.txt", header = TRUE, sep = "\t", row.names = 1)
  #rownames(data1) <- gsub("^X", "", rownames(data1))
  rownames(data1) <- gsub("-", ".", rownames(data1))
  
  # List of membership columns to process (replace these with your actual column names)
  membership_columns <- c("membership", "membership.1", "membership.2", "membership.3", "membership.4", "membership.5", "membership.6")
  
  # Open a PDF device to save all plots in one PDF file
  out_filename <- paste0("iteration-", iteration_file, ".pdf")
  pdf(out_filename)
  k=1
  # Loop through each membership column
  for (k in membership_columns) {
    # Sort data1 based on the current membership column (k)
    data1 <- data1[order(data1[, k]), ]
    
    # Read clinical data
    data2 <- read.table("Targeted_clinical.txt", header = TRUE, sep = "\t")
    
    # Match data based on sample ID
    m <- match(rownames(data1), data2[, 1])
    dat <- data2[m, ]
    dat$Cluster <- data1[, k]
    
    dat <- dat[, c("Case", "PFS_status", "PFS_years", "OS.censor", "OS.YEARS", "Cluster")]
    
    # OS analysis
    os <- dat[, c("Case", "OS.censor", "OS.YEARS", "Cluster")]
    os <- na.omit(os)
    
    # Create survival plot
    fit1 <- survfit(Surv(OS.YEARS, OS.censor) ~ Cluster, data = os)
    
    # Generate and save the survival plot
    print(ggsurvplot(fit1, pval = TRUE, data = os, risk.table = TRUE))
  }
  
  # Close the PDF device
  dev.off()
  
}


memb<-"membership.3"

library("dplyr")
setwd("/Users/ssharma/Downloads/DLBCL_Nat_Med_April_2018-master 2/consensus_clustering")

# Remove X and - in sample names
data1<-read.table("DLBCL.membership.txt", header=TRUE, sep="\t", row.names=1)
rownames(data1)<-gsub("^X", "", rownames(data1))
rownames(data1)<-gsub("-", ".", rownames(data1))
#rownames(data1)<-gsub(" ", "", rownames(data1))

#arrange according to membership column

##########################################
data1 <- data1[order(data1[, memb]), ]
##########################################
#iteration_file <- iteration_files[[iteration_file]]  # Replace with the actual dynamic file name

# Construct the full file path
file_path <- paste0("/Users/ssharma/Downloads/DLBCL_Nat_Med_April_2018-master 2/consensus_clustering/input_data/", iteration_file)

# Read the table
data2 <- read.table(file_path, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
m<-match(rownames(data1), colnames(data2))
data3<-data2[,m]



#Create a folder for storing results
###############################################################
name<-iteration_file
dir.create(paste0("/Users/ssharma/Downloads/DLBCL_Nat_Med_April_2018-master 2/consensus_clustering/", name, sep="_", memb))
setwd(paste0("/Users/ssharma/Downloads/DLBCL_Nat_Med_April_2018-master 2/consensus_clustering/", name, sep="_", memb))
###############################################################

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
  ################################################################  
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

df<-read.csv("Heatmap.csv", header=TRUE, sep=",", row.names=1, check.names=FALSE)


amplification_rows <- grepl("Amplification", rownames(df), ignore.case = TRUE)
df[amplification_rows, ] <- lapply(df[amplification_rows, ], function(x) {
  x[x == 2] <- "Broad_Amp;"
  x[x == 1] <- "Low_Amp;"
  x[x == 0] <- ""
  return(x)
})

#deletion
deletion_rows <- grepl("DeletionPeak", rownames(df), ignore.case = TRUE)
df[deletion_rows, ] <- lapply(df[deletion_rows, ], function(x) {
  x[x == 1] <- "Loss;"
  x[x == 0] <- ""
  return(x)
})

#Broad deletion
del <- grepl("Broad_Deletion", rownames(df), ignore.case = TRUE)
df[del, ] <- lapply(df[del, ], function(x) {
  x[x == 1] <- "Broad_Deletion;"
  x[x == 0] <- ""
  return(x)
})

rows_to_modify <- which(!(grepl("amplification|Deletion|Broad-Deletion|Broad_Gain", rownames(df), ignore.case = TRUE)))

# Replace values based on conditions
df[rows_to_modify, ] <- lapply(df[rows_to_modify, ], function(x) {
  x[x == 0] <- ""
  x[x == 1] <- "Mut;"
  return(x)
})

##Broad Gain
bg <- grepl("Broad_Gain", rownames(df), ignore.case = TRUE)
df[bg, ] <- lapply(df[bg, ], function(x) {
  x[x == 1] <- "Broad_Gain;"
  x[x == 0] <- ""
  return(x)
})

write.csv(df, file="check.csv")

library("ComplexHeatmap")

# Read the data from CSV files
data3 <- read.csv("check.csv", header = TRUE, sep = ",", check.names = FALSE, row.names = 1)
data4 <- read.csv("unique_across_groups.csv", header = TRUE, sep = ",", check.names = FALSE)

# Remove leading 'X' and replace '-' with '.' in row names
colnames(data3) <- gsub("^X", "", colnames(data3))
colnames(data3) <- gsub("-", ".", colnames(data3))

# Initialize vectors to store final column and row orders
final_cols <- c()
final_rows <- c()

# Loop through clusters
data1$clust<-data1[,memb]
k<-unique(data1$clust)
for (i in 1:length(k)) {
  # Subset data based on cluster
  subset_data <- data3[, data1$clust == k[i]]
  ids_to_extract <- which(data4$cluster == k[i])
  subset_data <- subset_data[ids_to_extract, ]
  
  # Define alter_fun_list and col as you did in the original code
  
  # Create oncoPrint plot
  h <- oncoPrint(subset_data, get_type = function(x) strsplit(x, ";")[[1]], alter_fun = alter_fun_list, col = col, remove_empty_columns = FALSE, row_names_gp = gpar(fontsize = 5), alter_fun_is_vectorized = FALSE, show_column_names = TRUE, column_names_gp = gpar(fontsize = 4))
  
  # Save oncoPrint plot as a PDF
  pdf_file <- paste0(k[i], ".pdf")
  pdf(pdf_file, height = 6, width = 6)
  print(h)  # Print the plot
  dev.off()
  
  # Extract text content from the PDF
  pdf_text_content <- pdf_text(pdf_file)
  cleaned_text <- gsub("\\s+", " ", pdf_text_content)
  
  # Extract columns and rows from the cleaned text and append to final_cols and final_rows
  cols <- colnames(data3)
  rows <- rownames(data3)
  extracted_cols <- unlist(str_extract_all(cleaned_text, paste(cols, collapse = "|")))
  extracted_rows <- unlist(str_extract_all(cleaned_text, paste(rows, collapse = "|")))
  final_cols <- c(final_cols, extracted_cols)
  final_rows <- c(final_rows, extracted_rows)
}

# Convert final_cols and final_rows to data frames
final_cols_df <- data.frame(final_cols)
final_rows_df <- data.frame(final_rows)

# Create the consolidated oncoPrint plot
pdf("new.pdf", height = 15, width = 20)
oncoPrint(data3, get_type = function(x) strsplit(x, ";")[[1]], alter_fun = alter_fun_list, col = col, remove_empty_columns = FALSE, remove_empty_rows = FALSE, row_names_gp = gpar(fontsize = 15), alter_fun_is_vectorized = FALSE, show_column_names = TRUE, column_names_gp = gpar(fontsize = 8), column_order = final_cols_df[, 1], row_order = final_rows_df[, 1])
dev.off()

#####
#########
############ Targeted Cohort ########
setwd("/Users/ssharma/Downloads/DLBCL_Nat_Med_April_2018-master 2/consensus_clustering")

data1<-read.csv("DLBCL.membership.txt", header=TRUE, sep="\t", row.names=1)
rownames(data1)<-gsub("-", ".", rownames(data1))
rownames(data1)<-gsub(".ATM.CN", "", rownames(data1))
rownames(data1)<-gsub(".ATM.TP53.BOTH", "", rownames(data1))
rownames(data1)<-gsub(".WT", "", rownames(data1))
rownames(data1)<-gsub(".ATM.MUT", "", rownames(data1))
rownames(data1)<-gsub(".ATM.DL", "", rownames(data1))
rownames(data1)<-gsub(".TP53.DL", "", rownames(data1))
rownames(data1)<-gsub(".TP53.MUT", "", rownames(data1))
rownames(data1)<-gsub(".TP53.CN", "", rownames(data1))


data2<-read.csv("Targeted_clinical.txt", header=TRUE, sep="\t")
# Match data based on sample ID
m <- match(rownames(data1), data2[, 1])
dat <- data2[m, ]
dat$Cluster <- data1$membership.4

dat <- dat[, c("Case",  "OS.YEARS", "OS.censor", "Cluster")]

# OS analysis

os <- na.omit(dat)

# Create survival plot
fit1 <- survfit(Surv(OS.YEARS, OS.censor) ~ Cluster, data = os)

# Generate and save the survival plot
pdf("member8.pdf")
print(ggsurvplot(fit1, pval = TRUE, data = os, risk.table = TRUE))
dev.off()




library("ComplexHeatmap")

# Read the data from CSV files
data3 <- read.csv("check.csv", header = TRUE, sep = ",", check.names = FALSE, row.names = 1)
colnames(data3)<-gsub("-", ".", colnames(data3))
colnames(data3)<-gsub(".ATM.CN", "", colnames(data3))
colnames(data3)<-gsub(".ATM.TP53.BOTH", "", colnames(data3))
colnames(data3)<-gsub(".WT", "", colnames(data3))
colnames(data3)<-gsub(".ATM.MUT", "", colnames(data3))
colnames(data3)<-gsub(".ATM.DL", "", colnames(data3))
colnames(data3)<-gsub(".TP53.DL", "", colnames(data3))
colnames(data3)<-gsub(".TP53.MUT", "", colnames(data3))
colnames(data3)<-gsub(".TP53.CN", "", colnames(data3))

data4 <- read.csv("unique_across_groups.csv", header = TRUE, sep = ",", check.names = FALSE)

# Remove leading 'X' and replace '-' with '.' in row names
rownames(data3) <- gsub("^X", "", rownames(data3))
rownames(data3) <- gsub("-", ".", rownames(data3))

# Initialize vectors to store final column and row orders
final_cols <- c()
final_rows <- c()

# Loop through clusters
data1$clust<-data1[,memb]
k<-unique(data1$clust)
for (i in 1:length(k)) {
  # Subset data based on cluster
  subset_data <- data3[, data1$clust == k[i]]
  ids_to_extract <- which(data4$cluster == k[i])
  subset_data <- subset_data[ids_to_extract, ]
  
  # Define alter_fun_list and col as you did in the original code
  
  # Create oncoPrint plot
  h <- oncoPrint(subset_data, get_type = function(x) strsplit(x, ";")[[1]], alter_fun = alter_fun_list, col = col, remove_empty_columns = FALSE, row_names_gp = gpar(fontsize = 5), alter_fun_is_vectorized = FALSE, show_column_names = TRUE, column_names_gp = gpar(fontsize = 4))
  
  # Save oncoPrint plot as a PDF
  pdf_file <- paste0(k[i], ".pdf")
  pdf(pdf_file, height = 6, width = 6)
  print(h)  # Print the plot
  dev.off()
  
  # Extract text content from the PDF
  pdf_text_content <- pdf_text(pdf_file)
  cleaned_text <- gsub("\\s+", " ", pdf_text_content)
  
  # Extract columns and rows from the cleaned text and append to final_cols and final_rows
  cols <- colnames(data3)
  rows <- rownames(data3)
  extracted_cols <- unlist(str_extract_all(cleaned_text, paste(cols, collapse = "|")))
  extracted_rows <- unlist(str_extract_all(cleaned_text, paste(rows, collapse = "|")))
  final_cols <- c(final_cols, extracted_cols)
  final_rows <- c(final_rows, extracted_rows)
}

# Convert final_cols and final_rows to data frames
final_cols_df <- data.frame(final_cols)
final_rows_df <- data.frame(final_rows)

# Create the consolidated oncoPrint plot
pdf("new.pdf", height = 15, width = 20)
oncoPrint(data3, get_type = function(x) strsplit(x, ";")[[1]], alter_fun = alter_fun_list, col = col, remove_empty_columns = FALSE, remove_empty_rows = FALSE, row_names_gp = gpar(fontsize = 15), alter_fun_is_vectorized = FALSE, show_column_names = TRUE, column_names_gp = gpar(fontsize = 8), column_order = final_cols_df[, 1], row_order = final_rows_df[, 1])
dev.off()
