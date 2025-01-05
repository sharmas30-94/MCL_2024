#Supplementary Figure 1a ####
setwd("/Users/sharmas30/Desktop/code for MCL")
library(survival)
library(forestplot)

# Load clinical data
data <- read.csv("clinical.csv", header = TRUE, sep = ",")

# Subset data for "Cohort-1"
data <- subset(data, data$COHORT_TYPE == "Cohort-1")

# Clinical parameters to analyze (short names used for visualization)
clinical_parameters <- c(
  "SOX11_IHC_GREATER_THAN_EQUAL_10_1_LESS_THAN_10_0",
  "MTP53_GREATER_THAN_EQUAL_90_1_ELSE_0",
  "Ki67_COUNTS_GREATER_THAN_30_1_ELSE_0",
  "SUBTYPES_CLASSICAL_100_0_B.P_1",
  "GROWTH_PATTERN_OTHER_0_DIFFUSE_100_1",
  "NODAL_1_EXTRANODAL_2_NODAL_WITH_EXTRANODAL_3",
  "NO_B_SYMPTOM_1_B_SYMPTOM_0",
  "WHETHER_TRANSPLANT_YES_0_NO_1",
  "MIPI",
  "GENDER_M_1_F_0",
  "AGE_GREATER_THAN_EQUAL_60_1_LESS_THAN_60_0",
  "RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0"
)

# Short names for visualization
short_names <- c(
  "SOX11 IHC",
  "TP53 IHC",
  "Ki-67 IHC",
  "Pathology Subtype",
  "Growth Pattern",
  "Nodal Status",
  "B Symptoms",
  "Transplant",
  "MIPI",
  "Sex",
  "Age",
  "Rituximab Administered"
)

# Initialize results dataframes
results_os <- data.frame(
  Parameter = character(),
  HR = numeric(),
  Lower95CI = numeric(),
  Upper95CI = numeric(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

results_pfs <- data.frame(
  Parameter = character(),
  HR = numeric(),
  Lower95CI = numeric(),
  Upper95CI = numeric(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

# Perform Cox proportional hazards regression for each parameter (OS and PFS)
for (i in seq_along(clinical_parameters)) {
  param <- clinical_parameters[i]
  short_name <- short_names[i]
  
  # OS Cox model
  formula_os <- as.formula(paste(
    "Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~", 
    param
  ))
  cox_os <- coxph(formula_os, data = data)
  summary_os <- summary(cox_os)
  
  # Extract OS results
  HR_os <- summary_os$coefficients[,"exp(coef)"]
  Lower95CI_os <- summary_os$conf.int[,"lower .95"]
  Upper95CI_os <- summary_os$conf.int[,"upper .95"]
  PValue_os <- summary_os$coefficients[,"Pr(>|z|)"]
  
  results_os <- rbind(
    results_os, 
    data.frame(
      Parameter = short_name,
      HR = HR_os,
      Lower95CI = Lower95CI_os,
      Upper95CI = Upper95CI_os,
      PValue = PValue_os
    )
  )
  
  # PFS Cox model
  formula_pfs <- as.formula(paste(
    "Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~", 
    param
  ))
  cox_pfs <- coxph(formula_pfs, data = data)
  summary_pfs <- summary(cox_pfs)
  
  # Extract PFS results
  HR_pfs <- summary_pfs$coefficients[,"exp(coef)"]
  Lower95CI_pfs <- summary_pfs$conf.int[,"lower .95"]
  Upper95CI_pfs <- summary_pfs$conf.int[,"upper .95"]
  PValue_pfs <- summary_pfs$coefficients[,"Pr(>|z|)"]
  
  results_pfs <- rbind(
    results_pfs, 
    data.frame(
      Parameter = short_name,
      HR = HR_pfs,
      Lower95CI = Lower95CI_pfs,
      Upper95CI = Upper95CI_pfs,
      PValue = PValue_pfs
    )
  )
}

# Remove duplicates (e.g., MIPI appears twice)
results_os <- results_os[!duplicated(results_os$Parameter), ]
results_pfs <- results_pfs[!duplicated(results_pfs$Parameter), ]

# Combine OS and PFS results
combined_results <- data.frame(
  Parameter = results_os$Parameter,
  OS_HR = sprintf("%.2f (%.2f-%.2f)", results_os$HR, results_os$Lower95CI, results_os$Upper95CI),
  OS_PValue = sprintf("%.3f", results_os$PValue),
  PFS_HR = sprintf("%.2f (%.2f-%.2f)", results_pfs$HR, results_pfs$Lower95CI, results_pfs$Upper95CI),
  PFS_PValue = sprintf("%.3f", results_pfs$PValue)
)

# Prepare forest plot table
forest_table <- rbind(
  c("Parameter", "OS HR (95% CI)", "OS P-Value", "PFS HR (95% CI)", "PFS P-Value"),
  cbind(
    combined_results$Parameter, 
    combined_results$OS_HR, 
    combined_results$OS_PValue, 
    combined_results$PFS_HR, 
    combined_results$PFS_PValue
  )
)

# Create the forest plot
pdf("Supllementary_Figure_1.pdf", height=15, width=20)
forestplot(
  labeltext = forest_table,
  mean = rbind(
    c(NA, NA), 
    cbind(results_os$HR, results_pfs$HR)
  ),
  lower = rbind(
    c(NA, NA), 
    cbind(results_os$Lower95CI, results_pfs$Lower95CI)
  ),
  upper = rbind(
    c(NA, NA), 
    cbind(results_os$Upper95CI, results_pfs$Upper95CI)
  ),
  is.summary = c(TRUE, rep(FALSE, nrow(results_os))),
  xlab = "Hazard Ratio",
  zero = 1,
  boxsize = 0.2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.1,
  col = fpColors(
    box = c("red", "darkblue"),
    line = "darkblue",
    summary = "royalblue"
  )
)
dev.off()


###Supplementary Figure 1b #####
setwd("/Users/sharmas30/Desktop/code for MCL")
library("survival")
library("survminer")
data<-read.csv("clinical.csv", header=TRUE, sep=",")

splots<-list()
d3<-subset(data, data$COHORT_TYPE == "Cohort-1")
d3<-d3[, c(1, 17, 12, 13)]
d3<-na.omit(d3)
fit <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ MTP53_GREATER_THAN_EQUAL_90_1_ELSE_0, data = d3)
splots[[1]] <- ggsurvplot(fit, 
                          pval = TRUE, 
                          data = d3, 
                          risk.table.col = "strata", 
                          surv.median.line = "hv", 
                          palette = c("purple", "green"), 
                          risk.table = TRUE, 
                          legend.title = "TP53 IHC",
                          legend.labs = c("Negative", "Positive"),
                          risk.table.height = 0.25,
                          title="OS",
                          risk.table.y.text = FALSE)

d23<-subset(data, data$COHORT_TYPE == "Cohort-1")
d23<-d23[, c(1, 32, 12, 13)]
d23<-na.omit(d23)
fit <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ RITUXIMAB_ADMINISTERED_PRIMARY_REGIMEN_YES_1_NO_0, data = d23)
splots[[2]] <- ggsurvplot(fit, 
                          pval = TRUE, 
                          data = d23, 
                          risk.table.col = "strata", 
                          surv.median.line = "hv", 
                          palette = c("#E4717A", "#C19A6B"), 
                          risk.table = TRUE, 
                          legend.title = "Rituximab adminstered (Primary Regimen: Cohort-1)",
                          legend.labs=c("No", "Yes"),
                          risk.table.height = 0.25,
                          title = "OS",
                          risk.table.y.text = FALSE)

pdf("Supplementary_Figure_1b.pdf", height=10, width=7)
arrange_ggsurvplots(splots, print = TRUE,
                    ncol = 1, nrow = 2, risk.table.height = 0.4)
dev.off()

### Supplementary Figure 1c ######
d15<-subset(data, data$COHORT_TYPE == "Cohort-1")
d15<-d15[, c(1, 24, 12, 13)]
d15<-na.omit(d15)
fit <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ NO_B_SYMPTOM_1_B_SYMPTOM_0, data = d15)
splots[[1]] <- ggsurvplot(fit, 
                           pval = TRUE, 
                           data = d15, 
                           risk.table.col = "strata", 
                           surv.median.line = "hv", 
                           palette = c("#E4D00A", "#C95A49"), 
                           risk.table = TRUE, 
                           legend.title = "B Symptom",
                           legend.labs=c("Yes", "No"),
                           risk.table.height = 0.25,
                           title = "OS",
                           risk.table.y.text = FALSE)

####Univariate analysis with PFS with NO_B_SYMPTOM_1_B_SYMPTOM_0"
d16<-subset(data, data$COHORT_TYPE == "Cohort-1")
d16<-d16[, c(1, 24, 14, 15)]
d16<-na.omit(d16)
fit <- survfit(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ NO_B_SYMPTOM_1_B_SYMPTOM_0, data = d16)
splots[[2]] <- ggsurvplot(fit, 
                           pval = TRUE, 
                           data = d16, 
                           risk.table.col = "strata", 
                           surv.median.line = "hv", 
                           palette = c("#E4D00A", "#C95A49"), 
                           risk.table = TRUE, 
                           legend.title = "B Symptom",
                           legend.labs=c("Yes", "No"),
                           risk.table.height = 0.25,
                           title = "PFS",
                           risk.table.y.text = FALSE)

### Univariate analysis Translant status with OS ###
d17<-subset(data, data$COHORT_TYPE == "Cohort-1")
d17<-d17[, c(1, 27, 12, 13)]
d17<-na.omit(d17)
fit <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ WHETHER_TRANSPLANT_YES_0_NO_1, data = d17)
splots[[3]] <- ggsurvplot(fit, 
                           pval = TRUE, 
                           data = d17, 
                           risk.table.col = "strata", 
                           surv.median.line = "hv", 
                           palette = c("#C88A65", "#B0BF1A"), 
                           risk.table = TRUE, 
                           legend.title = "Transplant Status",
                           legend.labs=c("Yes", "No"),
                           risk.table.height = 0.25,
                           title="OS",
                           risk.table.y.text = FALSE)


## Univariate analysis Transplant with PFS ###
d18<-subset(data, data$COHORT_TYPE == "Cohort-1")
d18<-d18[, c(1, 27, 14, 15)]
d18<-na.omit(d18)
fit <- survfit(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ WHETHER_TRANSPLANT_YES_0_NO_1, data = d18)
splots[[4]] <- ggsurvplot(fit, 
                           pval = TRUE, 
                           data = d18, 
                           risk.table.col = "strata", 
                           surv.median.line = "hv", 
                           palette = c("#C88A65", "#B0BF1A"), 
                           risk.table = TRUE, 
                           legend.title = "Transplant Status",
                           legend.labs=c("Yes", "No"),
                           risk.table.height = 0.25,
                           title="PFS",
                           risk.table.y.text = FALSE)



d21<-subset(data, data$COHORT_TYPE == "Cohort-1")
d21<-d21[, c(1, 31, 12, 13)]
d21<-na.omit(d21)
fit <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_STATUS_Dead.1_Alive_0) ~ AGE_GREATER_THAN_EQUAL_60_1_LESS_THAN_60_0, data = d21)
splots[[5]] <- ggsurvplot(fit, 
                           pval = TRUE, 
                           data = d21, 
                           risk.table.col = "strata", 
                           surv.median.line = "hv", 
                           palette = c("#79443B", "#D891EF"), 
                           risk.table = TRUE, 
                           legend.title = "Age (years)",
                           legend.labs=c("< 60", ">= 60"),
                           risk.table.height = 0.25,
                           title="OS",
                           risk.table.y.text = FALSE)

## Univariate analysis of age with PFS ##
d22<-subset(data, data$COHORT_TYPE == "Cohort-1")
d22<-d22[, c(1, 31, 14, 15)]
d22<-na.omit(d22)
fit <- survfit(Surv(PROGRESSION_FREE_SURVIVAL_MONTHS, PROGRESSION_FREE_SURVIVAL_STATUS_YES_1_NO_0) ~ AGE_GREATER_THAN_EQUAL_60_1_LESS_THAN_60_0, data = d22)
splots[[6]] <- ggsurvplot(fit, 
                           pval = TRUE, 
                           data = d22, 
                           risk.table.col = "strata", 
                           surv.median.line = "hv", 
                           palette = c("#79443B", "#D891EF"), 
                           risk.table = TRUE, 
                           legend.title = "Age (years)",
                           legend.labs=c("< 60", ">=60"),
                           risk.table.height = 0.25,
                           title="PFS",
                           risk.table.y.text = FALSE)


pdf("Supplementary_Figure_1c.pdf", height=20, width=15)
arrange_ggsurvplots(splots, print = TRUE,
                    ncol = 2, nrow = 3, risk.table.height = 0.4)
dev.off()