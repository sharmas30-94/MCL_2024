# === Figure 4a: Predicting NMF Clusters Using Random Forest ===

# Set working directory
setwd("/Users/sharmas30/Downloads")

# Load required package
if (!require("randomForest")) install.packages("randomForest")
library(randomForest)
library("caret")

# === Load Data ===

# Cohort-1: WES features + NMF cluster
wes <- read.csv("WES_random_forest.csv", header = TRUE)
wes$Cluster <- as.factor(wes$Cluster)

# Cohort-2: Targeted panel features
cohort2 <- read.csv("Targeted_randomForest_new.csv", header = TRUE)
names(cohort2)[1] <- "Sample"  # rename first column to 'Sample' if needed

# === Define Feature Set Used in Best Model ===
common_features <- c(
  "Deletion_6q21", "SP140", "Deletion_8p23.3", "NSD2", "Deletion_TP53",
  "Amplification_13q31.3", "Amplification_7p22.3", "KMT2D", "Deletion_19p13.3",
  "Deletion_12p13.2", "Deletion_11q22.3", "Deletion_13q12.11", "Deletion_1p13.1",
  "TET2", "Deletion_9p21.1", "ARID1A", "SMARCA4", "SPEN", "Deletion_13q14.2",
  "CARD11", "Deletion_6q25.3", "Deletion_9p21.3", "Amplification-8q23.2-8q24.21",
  "Deletion_17p13.3", "TP53", "ATM", "Amplification-3q22.1-3q22.3",
  "Amplification_17q24.3"
)

# === Check which features are present in WES ===
available_features <- intersect(common_features, names(wes))

# === Prepare Training Data (Cohort-1) ===
X_train <- wes[, available_features]
y_train <- wes$Cluster

# === Prepare Test Data (Cohort-2) ===

# Add missing features to Cohort-2 and fill with 0
missing_in_cohort2 <- setdiff(available_features, names(cohort2))
for (col in missing_in_cohort2) {
  cohort2[[col]] <- 0
}

# Align columns to training features
X_test <- cohort2[, available_features]

# === Train Random Forest Model ===
set.seed(42)
rf_model <- randomForest(x = X_train, y = y_train, ntree = 500)

# === Predict Clusters for Cohort-2 ===
predicted_clusters <- predict(rf_model, newdata = X_test)
cohort2$Predicted_Cluster <- predicted_clusters

# === Output Results ===
write.csv(cohort2, "Cohort2_with_Predicted_Clusters.csv", row.names = FALSE)

# === Preview Output ===
head(cohort2[, c("Sample", "Predicted_Cluster")])

library(caret)

# true vs predicted
conf <- confusionMatrix(data = predicted, reference = actual)

# View overall metrics
print(conf)

# To extract per-class stats:
as.data.frame(conf$byClass)

