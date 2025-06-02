#Figure 2
setwd("/Users/sharmas30/Desktop/code for MCL/Github_code/MCL_2024/Figure2/Data")
library(ggplot2)
finalDat <-read.delim("WES.Targeted.Final.Freq.Plot.Data.052825.txt")
chr=seq(1,22,by=1)
chr.ends<-c("248956422", "491149951", "689445510", "879660065", "1061198324", "1232004303", "1391350276", "1536488912", "1674883629", "1808681051", "1943767673", "2077042982", "2191407310", "2298451028", "2400442217", "2490780562", "2574038003", "2654411288", "2713028904", "2777473071", "2824183054", "2875001522")
chr.ends <-as.numeric(chr.ends)
midpt<- c("124478211", "370053186.5", "590297730.5", "784552787.5", "970429194.5", "1146601313.5", "1311677289.5", "1463919594", "1605686270.5", "1741782340", "1876224362", "2010405327.5", "2134225146", "2244929169", "2349446622.5", "2445611389.5", "2532409282.5", "2614224645.5", "2683720096", "2745250987.5", "2800828062.5", "2849592288")
midpt<-as.numeric(midpt)
plotWEStargeted<-ggplot()+geom_polygon(data= finalDat,aes(x=pos,y= WES.gain.freq),fill="#FF8785",alpha=.5)+ geom_polygon(data= finalDat,aes(x=pos,y= WES.loss.freq),fill="#60E4FF",alpha=.5) +geom_line(data= finalDat,aes(x=pos,y= targeted.gain.freq),color="black")+ geom_line(data= finalDat,aes(x=pos,y= targeted.loss.freq),color="black") +geom_vline(xintercept=chr.ends, linewidth=.65)+theme_bw()+ labs(y ="% of cases with abnormality",x="") + scale_x_continuous(expand = c(0, 0) )+geom_text(aes(x=midpt,y=50,label=chr))+ggtitle("") +coord_cartesian(ylim=c(-60,60)) +theme_classic()+ theme(axis.title.x=element_blank(),  axis.text.x=element_blank(), axis.ticks.x=element_blank())
pdf(file=paste("Final.Plot.Fig2A",Sys.Date(),"pdf",sep="."),width=9.5,height=3.5,useDingbats=F,onefile=F)
plotWEStargeted
dev.off()


### Figure 2b ####
setwd("/Users/sharmas30/Desktop/code for MCL")
data1<-read.csv("Gain_driven.csv", header=TRUE, sep=",")
pdf("Gain_Driven_OS.pdf")
pval<-c()
gene<-c()

for(i in 7:ncol(data1)){
  df1<-data1[, c(1,2,3,6)]
  df2<-data1[,i]
  k<-colnames(data1)[i]
  one<-cbind(df1, df2)
  one<-na.omit(one)
  one<-subset(one, one$RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0 == 1)
  one$group <- ifelse(one[, "df2"] > 0, "Gain", "Other")
  fit <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ group, data = one)
  val<- surv_pvalue(fit)$pval
  pval<-c(pval, val)
  gene<-c(gene, k)
  print(ggsurvplot(fit, pval = TRUE, data = one, risk.table = TRUE, title = colnames(data1)[i]))
}

df <- data.frame(gene, pval)
write.csv(df, file="Gain_gene_pval_OS.csv")

dev.off()


data1<-read.csv("Gain_driven.csv", header=TRUE, sep=",")
pdf("Gain_Driven_PFS.pdf")
pval<-c()
gene<-c()

for(i in 7:ncol(data1)){
  df1<-data1[, c(1,4,5,6)]
  df2<-data1[,i]
  k<-colnames(data1)[i]
  one<-cbind(df1, df2)
  one<-na.omit(one)
  one<-subset(one, one$RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0 == 1)
  one$group <- ifelse(one[, "df2"] > 0, "Gain", "Other")
  fit <- survfit(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ group, data = one)
  val<- surv_pvalue(fit)$pval
  pval<-c(pval, val)
  gene<-c(gene, k)
  print(ggsurvplot(fit, pval = TRUE, data = one, risk.table = TRUE, title = colnames(data1)[i]))
}

df <- data.frame(gene, pval)
write.csv(df, file="Gain_gene_pval_pfS.csv")

dev.off()

data1<-read.csv("Loss_driven.csv", header=TRUE, sep=",")
pdf("Loss_Driven_OS.pdf")
pval<-c()
gene<-c()

for(i in 7:ncol(data1)){
  df1<-data1[, c(1,2,3,6)]
  df2<-data1[,i]
  k<-colnames(data1)[i]
  one<-cbind(df1, df2)
  one<-na.omit(one)
  one<-subset(one, one$RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0 == 1)
  one$group <- ifelse(one[, "df2"] < 0, "Loss", "Other")
  fit <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ group, data = one)
  val<- surv_pvalue(fit)$pval
  pval<-c(pval, val)
  gene<-c(gene, k)
  print(ggsurvplot(fit, pval = TRUE, data = one, risk.table = TRUE, title = colnames(data1)[i]))
}

df <- data.frame(gene, pval)
write.csv(df, file="Loss_gene_pval_OS.csv")

dev.off()


data1<-read.csv("Loss_driven.csv", header=TRUE, sep=",")
pdf("Loss_Driven_PFS.pdf")
pval<-c()
gene<-c()

for(i in 7:ncol(data1)){
  df1<-data1[, c(1,4,5,6)]
  df2<-data1[,i]
  k<-colnames(data1)[i]
  one<-cbind(df1, df2)
  one<-na.omit(one)
  one<-subset(one, one$RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0 == 1)
  one$group <- ifelse(one[, "df2"] < 0, "Loss", "Other")
  fit <- survfit(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ group, data = one)
  val<- surv_pvalue(fit)$pval
  pval<-c(pval, val)
  gene<-c(gene, k)
  print(ggsurvplot(fit, pval = TRUE, data = one, risk.table = TRUE, title = colnames(data1)[i]))
}

df <- data.frame(gene, pval)
write.csv(df, file="Loss_gene_pval_PFS.csv")

dev.off()


### Figure 2c ####
### GSEA Run on GenePatterns for pathway analysis fo cases with 15q21.2 amps and cases with no 15q21.2 amps
### BCL2L10 OS and PFS analysis ####



### Figure 2d ###
data<-read.csv("loss_median_centered_wes.csv", header=TRUE, sep=",")
data<-na.omti(data)
data$ATM.status<-factor(data$ATM.status, levels=c("Mut_CN", "WT", "CN", "Mut"))
ATM<-ggplot(data, aes(x=reorder(GENE, ATM.median.centered), y=ATM.median.centered, label="", hjust=0.8)) + geom_bar(stat="identity", position="identity", aes(fill=ATM.status)) + scale_fill_manual(values=c("grey", "blue", "forestgreen", "purple"))+ theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("Patients") + ylab("Median Centered GEP")

data$TP53.status<-factor(data$TP53.status, levels=c("Mut_CN", "WT", "CN", "Mut"))
TP53<-ggplot(data, aes(x=reorder(GENE, TP53.median.centered), y=TP53.median.centered, label="", hjust=0.8)) + geom_bar(stat="identity", position="identity", aes(fill=TP53.status)) + scale_fill_manual(values=c("grey", "blue", "forestgreen", "purple"))+ theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("Patients") + ylab("Median Centered GEP")

data<-read.csv("Gain_median_centered_cn_values.csv", header=TRUE, sep=",")
data$MYC.cn<-factor(data$MYC.cn, levels=c(0,1))
MYC<-ggplot(data, aes(x=reorder(GENE, MYC.gep.median.centered), y=MYC.gep.median.centered, label="", hjust=0.8)) + geom_bar(stat="identity", position="identity", aes(fill=MYC.cn)) + scale_fill_manual(values=c("grey", "firebrick1"))+ theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("Patients") + ylab("Median Centered GEP")

data$DVL3.cn<-factor(data$DVL3.cn, levels=c(0,1))
DVL3<-ggplot(data, aes(x=reorder(GENE, DVL3.gep.median.centered), y=DVL3.gep.median.centered, label="", hjust=0.8)) + geom_bar(stat="identity", position="identity", aes(fill=DVL3.cn)) + scale_fill_manual(values=c("grey", "firebrick1"))+ theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("Patients") + ylab("Median Centered GEP")
pdf("Fig2d.pdf", height=10, width=10)
ggarrange(ATM, TP53, MYC, DVL3, ncol=1, row=4)
dev.off()


####. Figure 2e ######
setwd("")
# Load required libraries
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# Read data
data1 <- read.csv("CN_lesion.csv", header = TRUE, sep = ",")
data2 <- read.csv("clinical.csv", header = TRUE, sep = ",")

# Merge clinical data into lesion dataset
data1$COHORT_TYPE <- data2$COHORT_TYPE[match(data1[,1], data2[,1])]
data1$OVERALL_SURVIVAL_MONTHS <- data2$OVERALL_SURVIVAL_MONTHS[match(data1[,1], data2[,1])]
data1$OVERALL_SURVIVAL_STATUS <- data2$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0[match(data1[,1], data2[,1])]
data1$PROGRESSION_FREE_SURVIVAL_MONTHS <- data2$PROGRESSION_FREE_SURVIVAL_MONTHS[match(data1[,1], data2[,1])]
data1$PROGRESSION_FREE_SURVIVAL_STATUS <- data2$PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0[match(data1[,1], data2[,1])]
data1$RITUXIMAB <- data2$RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0[match(data1[,1], data2[,1])]

# Convert survival status to numeric
data1$OVERALL_SURVIVAL_STATUS <- as.numeric(data1$OVERALL_SURVIVAL_STATUS)
data1$PROGRESSION_FREE_SURVIVAL_STATUS <- as.numeric(data1$PROGRESSION_FREE_SURVIVAL_STATUS)
data1$RITUXIMAB <- as.numeric(data1$RITUXIMAB)

# Filter only Rituximab positive cases
data1_filtered <- subset(data1, RITUXIMAB == 1)

# Select lesion columns
lesion_columns <- grep("^Amplification_|^Deletion_", colnames(data1_filtered), value = TRUE)

# Keep only relevant columns
data1_filtered <- data1_filtered[, c("OVERALL_SURVIVAL_MONTHS", "OVERALL_SURVIVAL_STATUS",
                                     "PROGRESSION_FREE_SURVIVAL_MONTHS", "PROGRESSION_FREE_SURVIVAL_STATUS", lesion_columns)]

# Remove rows with NA values in survival columns
data1_filtered <- na.omit(data1_filtered)

# Set output directory for plots
output_dir <- "Survival_Analysis_Plots_Rituximab_Positive"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Initialize data frame for storing results
significant_results <- data.frame(Lesion = character(),
                                  HR = numeric(),
                                  Lower_CI = numeric(),
                                  Upper_CI = numeric(),
                                  P_Value = numeric(),
                                  Type = character(),
                                  stringsAsFactors = FALSE)

# Function to perform survival analysis
perform_survival_analysis <- function(time_col, status_col, analysis_type) {
  for (lesion in lesion_columns) {
    
    # Ensure lesion is a factor with two levels (Present/Absent)
    data1_filtered[[lesion]] <- as.factor(ifelse(data1_filtered[[lesion]] > 0, "Present", "Absent"))
    
    # Cox proportional hazards model
    cox_model <- coxph(Surv(data1_filtered[[time_col]], data1_filtered[[status_col]]) ~ get(lesion), data = data1_filtered)
    cox_summary <- summary(cox_model)
    
    # Extract log HR and standard error
    log_hr <- cox_summary$coefficients[,"coef"]
    se <- cox_summary$coefficients[,"se(coef)"]
    
    # Compute HR and its 95% confidence intervals using the correct formula
    hr <- exp(log_hr)
    lower_ci <- exp(log_hr - 1.96 * se)
    upper_ci <- exp(log_hr + 1.96 * se)
    p_value <- cox_summary$coefficients[,"Pr(>|z|)"]
    
    # Ensure HR is between Lower_CI and Upper_CI (detect potential issues)
    if (!(lower_ci < hr & hr < upper_ci)) {
      print(paste("Warning: HR is outside CI bounds for lesion:", lesion))
      print(paste("HR:", hr, "Lower_CI:", lower_ci, "Upper_CI:", upper_ci))
    }
    
    # Save result
    significant_results <<- rbind(significant_results, data.frame(Lesion = lesion, HR = hr, 
                                                                  Lower_CI = lower_ci, Upper_CI = upper_ci, 
                                                                  P_Value = p_value, Type = analysis_type))
  }
}

# Perform OS analysis
perform_survival_analysis("OVERALL_SURVIVAL_MONTHS", "OVERALL_SURVIVAL_STATUS", "OS")

# Perform PFS analysis
perform_survival_analysis("PROGRESSION_FREE_SURVIVAL_MONTHS", "PROGRESSION_FREE_SURVIVAL_STATUS", "PFS")

# Print the results table
print("Cox Model Outputs (HR, Lower CI, Upper CI, P-Value):")
print(significant_results)

# Filter for significant values (P < 0.05)
filtered_results <- significant_results %>% filter(P_Value < 0.05)
values <- filtered_results$Lesion


matched_results <- significant_results %>%
  filter(Lesion %in% values) %>%
  arrange(Lesion, Type)


matched_results <- significant_results %>%
  filter(Lesion %in% k) %>%
  arrange(Lesion, Type)


# Define colors for OS and PFS
color_map <- c("OS" = "Magenta", "PFS" = "Green")

# Ensure data is sorted correctly: OS appears first, then PFS for each lesion
matched_results <- matched_results %>%
  arrange(Lesion, Type)

# Assign y-axis positions ensuring OS is plotted first, then PFS immediately after
matched_results <- matched_results %>%
  mutate(y_position = as.numeric(factor(Lesion, levels = rev(unique(Lesion)))) * 2 - ifelse(Type == "OS", 0.6, 0.2))

# Define significance annotation (* for p < 0.05, "ns" otherwise)
matched_results <- matched_results %>%
  mutate(Significance = ifelse(P_Value < 0.05, "*", "ns"))

# Create the forest plot with significance annotations
ggplot(matched_results, aes(x = HR, y = y_position, color = Type)) +
  geom_point(size = 3) +  # Plot points for HR
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +  # Add error bars
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  # Reference line at HR = 1
  geom_text(aes(label = Significance), hjust = -0.5, vjust = 0.5, size = 4) +  # Add significance labels
  scale_color_manual(values = color_map) +
  labs(title = "Rituximab Administered",
       x = "Hazard Ratio (HR)",
       y = "Lesion") +
  scale_y_continuous(breaks = seq(1, max(matched_results$y_position), by = 2),
                     labels = unique(matched_results$Lesion)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")


#Rituximab naive cases
data1 <- read.csv("CN_lesion.csv", header = TRUE, sep = ",")
data2 <- read.csv("clinical.csv", header = TRUE, sep = ",")

# Merge clinical data into lesion dataset
data1$COHORT_TYPE <- data2$COHORT_TYPE[match(data1[,1], data2[,1])]
data1$OVERALL_SURVIVAL_MONTHS <- data2$OVERALL_SURVIVAL_MONTHS[match(data1[,1], data2[,1])]
data1$OVERALL_SURVIVAL_STATUS <- data2$OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0[match(data1[,1], data2[,1])]
data1$PROGRESSION_FREE_SURVIVAL_MONTHS <- data2$PROGRESSION_FREE_SURVIVAL_MONTHS[match(data1[,1], data2[,1])]
data1$PROGRESSION_FREE_SURVIVAL_STATUS <- data2$PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0[match(data1[,1], data2[,1])]
data1$RITUXIMAB <- data2$RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0[match(data1[,1], data2[,1])]

# Convert survival status to numeric
data1$OVERALL_SURVIVAL_STATUS <- as.numeric(data1$OVERALL_SURVIVAL_STATUS)
data1$PROGRESSION_FREE_SURVIVAL_STATUS <- as.numeric(data1$PROGRESSION_FREE_SURVIVAL_STATUS)
data1$RITUXIMAB <- as.numeric(data1$RITUXIMAB)

# Filter only Rituximab positive cases
data1_filtered <- subset(data1, RITUXIMAB == 0)

# Select lesion columns
lesion_columns <- grep("^Amplification_|^Deletion_", colnames(data1_filtered), value = TRUE)

# Keep only relevant columns
data1_filtered <- data1_filtered[, c("OVERALL_SURVIVAL_MONTHS", "OVERALL_SURVIVAL_STATUS",
                                     "PROGRESSION_FREE_SURVIVAL_MONTHS", "PROGRESSION_FREE_SURVIVAL_STATUS", lesion_columns)]

# Remove rows with NA values in survival columns
data1_filtered <- na.omit(data1_filtered)

# Set output directory for plots
output_dir <- "Survival_Analysis_Plots_Rituximab_Positive"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Initialize data frame for storing results
significant_results <- data.frame(Lesion = character(),
                                  HR = numeric(),
                                  Lower_CI = numeric(),
                                  Upper_CI = numeric(),
                                  P_Value = numeric(),
                                  Type = character(),
                                  stringsAsFactors = FALSE)

# Function to perform survival analysis
perform_survival_analysis <- function(time_col, status_col, analysis_type) {
  for (lesion in lesion_columns) {
    
    # Ensure lesion is a factor with two levels (Present/Absent)
    data1_filtered[[lesion]] <- as.factor(ifelse(data1_filtered[[lesion]] > 0, "Present", "Absent"))
    
    # Cox proportional hazards model
    cox_model <- coxph(Surv(data1_filtered[[time_col]], data1_filtered[[status_col]]) ~ get(lesion), data = data1_filtered)
    cox_summary <- summary(cox_model)
    
    # Extract log HR and standard error
    log_hr <- cox_summary$coefficients[,"coef"]
    se <- cox_summary$coefficients[,"se(coef)"]
    
    # Compute HR and its 95% confidence intervals using the correct formula
    hr <- exp(log_hr)
    lower_ci <- exp(log_hr - 1.96 * se)
    upper_ci <- exp(log_hr + 1.96 * se)
    p_value <- cox_summary$coefficients[,"Pr(>|z|)"]
    
    # Ensure HR is between Lower_CI and Upper_CI (detect potential issues)
    if (!(lower_ci < hr & hr < upper_ci)) {
      print(paste("Warning: HR is outside CI bounds for lesion:", lesion))
      print(paste("HR:", hr, "Lower_CI:", lower_ci, "Upper_CI:", upper_ci))
    }
    
    # Save result
    significant_results <<- rbind(significant_results, data.frame(Lesion = lesion, HR = hr, 
                                                                  Lower_CI = lower_ci, Upper_CI = upper_ci, 
                                                                  P_Value = p_value, Type = analysis_type))
  }
}

# Perform OS analysis
perform_survival_analysis("OVERALL_SURVIVAL_MONTHS", "OVERALL_SURVIVAL_STATUS", "OS")

# Perform PFS analysis
perform_survival_analysis("PROGRESSION_FREE_SURVIVAL_MONTHS", "PROGRESSION_FREE_SURVIVAL_STATUS", "PFS")

# Print the results table
print("Cox Model Outputs (HR, Lower CI, Upper CI, P-Value):")
print(significant_results)

matched_resuls<-significant_results


# Ensure PFS appears first
matched_results$Type <- factor(matched_results$Type, levels = c("PFS", "OS"))

# **Custom lesion order**
custom_order <- c("Deletion 9p21.3", "Deletion 9p21.1", 
                  "Deletion 1p22.1", "Deletion 1p21.2", 
                  "Amplification 18q21.33", "Amplification 14q32.33")

# Assign ordered factor levels
matched_results$Lesion <- factor(matched_results$Lesion, levels = rev(custom_order))

# Create y-axis positioning: PFS appears first, then OS for each lesion
matched_results <- matched_results %>%
  arrange(Lesion, Type) %>%  # Apply the custom order
  mutate(y_position = as.numeric(factor(Lesion, levels = rev(custom_order))) * 2 - as.numeric(Type == "OS"))

# Define colors: OS (magenta) and PFS (green)
color_mapping <- c("OS" = "Magenta", "PFS" = "Green")

# Create the forest plot with custom order and theme_bw()
ggplot(matched_results, aes(x = HR, y = y_position, color = Type)) +
  geom_point(size = 3) +  # Points for HR
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.15) +  # Error bars
  scale_x_log10() +  # Log scale for HR
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  # Reference line at HR = 1
  scale_color_manual(values = color_mapping) +  # Apply custom color scheme
  labs(title = "Rituximab Administered Cases",
       x = "Hazard Ratio (HR)", y = "Lesion", color = "Type") +
  theme_bw() +  # Change to black and white theme
  scale_y_continuous(breaks = seq(1, max(matched_results$y_position), 2), 
                     labels = rev(custom_order)) +  # Use custom order for y-axis
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 14))


# Figure 2f ####
# Load necessary libraries
library(survival)
library(survminer)

# Read the dataset
data1 <- read.csv("plot_signature.csv", header = TRUE, sep = ",")

# Initialize vectors to store Cox regression results
gene_cox <- c()
cox_pval <- c()
cox_hr <- c()  # Hazard Ratio
cox_lower_ci <- c()  # Lower CI for HR
cox_upper_ci <- c()  # Upper CI for HR

# Loop through each gene column (assuming gene expression data starts from column 7)

for (i in 7:ncol(data1)) {
  
  # Extract the gene name (column name)
  gene_name <- colnames(data1)[i]
  
  # Subset the data for the current gene expression column
  one_gene_data <- data1[, c("OVERALL_SURVIVAL_MONTHS", "OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0", "RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0", gene_name)]
  
  # Remove rows with missing values
  one_gene_data <- na.omit(one_gene_data)
  
  # Subset data for Rituximab cases only
  one_gene_data <- subset(one_gene_data, one_gene_data$RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0 == 1)
  
  # Create a new group variable based on gene expression (Loss or Other)
  one_gene_data$group <- ifelse(one_gene_data[, gene_name] < 0, "Loss", "Other")
  
  # Ensure "Other" is the reference category for the Cox regression
  one_gene_data$group <- factor(one_gene_data$group, levels = c("Other", "Loss"))
  
  # Fit the Cox model
  cox_model <- coxph(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ group, data = one_gene_data)
  
  # Extract Cox regression summary
  cox_summary <- summary(cox_model)
  
  # Extract necessary statistics from the Cox model
  
  pval_cox <- cox_summary$coefficients[1, "Pr(>|z|)"]
  hr_cox <- cox_summary$coefficients[1, "exp(coef)"]
  lower_ci_cox <- cox_summary$conf.int[1, "lower .95"]
  upper_ci_cox <- cox_summary$conf.int[1, "upper .95"]
  
  # Store the results
  gene_cox <- c(gene_cox, gene_name)
  cox_pval <- c(cox_pval, pval_cox)
  cox_hr <- c(cox_hr, hr_cox)
  cox_lower_ci <- c(cox_lower_ci, lower_ci_cox)
  cox_upper_ci <- c(cox_upper_ci, upper_ci_cox)
}

# Combine the results into a data frame
cox_results <- data.frame(
  gene = gene_cox,
  pval = cox_pval,
  hazard_ratio = cox_hr,
  lower_ci = cox_lower_ci,
  upper_ci = cox_upper_ci
)

# View the results of the Cox regression
print(cox_results)




library(ggplot2)
library(dplyr)

# Assuming cox_results is already generated, as per the previous code

# Modify color column based on p-value
cox_results$color <- ifelse(cox_results$pval < 0.05, "blue", "grey")

# Add stars or dots based on p-value
cox_results$significance <- ifelse(cox_results$pval < 0.05, "*", ".")

# Create a forest plot with distance for the stars, bold text, and a dotted line at HR = 1

pdf("Fig_2f_OS.pdf")
ggplot(cox_results, aes(x = hazard_ratio, y = gene, xmin = lower_ci, xmax = upper_ci, color = color)) +
  geom_point(size = 4) + # Points for the HR, colored based on p-value
  geom_errorbarh(height = 0.2, color = "black") + # Horizontal error bars for CI (black color)
  scale_color_identity() + # Use the color as specified in the data (no scale mapping)
  geom_text(aes(label = significance), 
            hjust = -0.3, 
            nudge_y = 0.2,   # Nudge the text labels vertically to create distance
            nudge_x = 0.1,   # Nudge the text horizontally to create more space from points
            size = 5, 
            color = "black",  # Color of the stars/dots
            fontface = "bold") + # Make the stars and dots bold
  geom_vline(xintercept = 1, linetype = "dotted", color = "black", size = 1) + # Dotted line at HR = 1
  theme_minimal() + # Minimal theme for the plot
  labs(x = "Hazard Ratio (HR)", y = "Prognostic Loss Driven Genes") +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(legend.position = "none") + theme_bw() # Remove legend for clarity

dev.off()

### PFS 
library(survival)
library(survminer)

# Read the dataset
data1 <- read.csv("plot_signature.csv", header = TRUE, sep = ",")

# Initialize vectors to store Cox regression results
gene_cox <- c()
cox_pval <- c()
cox_hr <- c()  # Hazard Ratio
cox_lower_ci <- c()  # Lower CI for HR
cox_upper_ci <- c()  # Upper CI for HR

# Loop through each gene column (assuming gene expression data starts from column 7)

for (i in 7:ncol(data1)) {
  
  # Extract the gene name (column name)
  gene_name <- colnames(data1)[i]
  
  # Subset the data for the current gene expression column
  one_gene_data <- data1[, c("PROGRESSION_FREE_SURVIVAL_MONTHS", "PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0", "RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0", gene_name)]
  
  # Remove rows with missing values
  one_gene_data <- na.omit(one_gene_data)
  
  # Subset data for Rituximab cases only
  one_gene_data <- subset(one_gene_data, one_gene_data$RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0 == 1)
  
  # Create a new group variable based on gene expression (Loss or Other)
  one_gene_data$group <- ifelse(one_gene_data[, gene_name] < 0, "Loss", "Other")
  
  # Ensure "Other" is the reference category for the Cox regression
  one_gene_data$group <- factor(one_gene_data$group, levels = c("Other", "Loss"))
  
  # Fit the Cox model
  cox_model <- coxph(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ group, data = one_gene_data)
  
  # Extract Cox regression summary
  cox_summary <- summary(cox_model)
  
  # Extract necessary statistics from the Cox model
  
  pval_cox <- cox_summary$coefficients[1, "Pr(>|z|)"]
  hr_cox <- cox_summary$coefficients[1, "exp(coef)"]
  lower_ci_cox <- cox_summary$conf.int[1, "lower .95"]
  upper_ci_cox <- cox_summary$conf.int[1, "upper .95"]
  
  # Store the results
  gene_cox <- c(gene_cox, gene_name)
  cox_pval <- c(cox_pval, pval_cox)
  cox_hr <- c(cox_hr, hr_cox)
  cox_lower_ci <- c(cox_lower_ci, lower_ci_cox)
  cox_upper_ci <- c(cox_upper_ci, upper_ci_cox)
}

# Combine the results into a data frame
cox_results <- data.frame(
  gene = gene_cox,
  pval = cox_pval,
  hazard_ratio = cox_hr,
  lower_ci = cox_lower_ci,
  upper_ci = cox_upper_ci
)

# View the results of the Cox regression
print(cox_results)




library(ggplot2)
library(dplyr)

# Assuming cox_results is already generated, as per the previous code

# Modify color column based on p-value
cox_results$color <- ifelse(cox_results$pval < 0.05, "blue", "grey")

# Modify significance labels based on p-value
cox_results_combined$significance <- ifelse(cox_results_combined$pval < 0.05, "*", 
                                            ifelse(cox_results_combined$pval < 0.1, ".", "NS"))

pdf("Fig_2f_PFS.pdf")
ggplot(cox_results_combined, aes(x = hazard_ratio, y = gene, xmin = lower_ci, xmax = upper_ci)) +
  geom_point(size = 4, color = "blue") +  # Keep all circles blue
  geom_errorbarh(height = 0.2, color = "black") +  # Horizontal error bars for CI (black color)
  geom_text(aes(label = significance, color = significance), 
            hjust = -0.3, 
            nudge_y = 0.2,   # Nudge the text labels vertically to create distance
            nudge_x = 0.1,   # Nudge the text horizontally to create more space from points
            size = 5, 
            fontface = "bold") + # Make the stars and dots bold
  geom_vline(xintercept = 1, linetype = "dotted", color = "red", size = 1) +  # Dotted red line at HR = 1
  facet_wrap(~ type, scales = "free_y") +  # Separate the OS and PFS into two panels
  scale_color_manual(values = c("*" = "blue", "." = "grey", "NS" = "black")) +  # Set custom colors for significance
  labs(x = "Hazard Ratio (HR)", y = "Gene", title = "Forest Plot") +
  theme_minimal() +  # Minimal theme for the plot
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "bottom") +  # Move legend to the bottom
  theme(legend.title = element_blank(),  # Remove legend title
        legend.key = element_blank(),  # Remove legend key
        legend.text = element_text(size = 12)) +  # Make legend text a bit larger
  theme(legend.position = "bottom") + theme_bw()# Show legend at the bottom


dev.off()


### Figure 2g #####


