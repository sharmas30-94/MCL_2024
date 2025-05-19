# IMPORT DATA & PERFORM NORMALIZATION -------------------------

# Install the remotes package if it is not already installed
if (!require(remotes)) {
  install.packages("remotes")
  library(remotes)
}

# Install the here package if it is not already installed
if (!require(here)) {
  install.packages("here")
  library(here)
}

# Load or install & load required packages with specific versions
install_and_load <- function(pkg, version) {
  if (!require(pkg, character.only = TRUE)) {
    remotes::install_version(pkg, version = version)
    library(pkg, character.only = TRUE)
  }
}

install_and_load("cyCombine", "0.2.15")
install_and_load("dplyr", "1.1.3")
install_and_load("tidyr", "1.3.0")

current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
data_master <- file.path(current_dir, "dataset")
meta_master <- file.path(current_dir, "dataset meta")

panel_file <- file.path(meta_master, "panel.csv")
data_dir <- data_master
meta_dir <- meta_master
meta_file <- file.path(meta_dir, "metadata.csv")


markers <- read.csv(panel_file) %>%
  filter(Type != "none") %>%
  pull(Marker)


uncorrected <- prepare_data(
  data_dir = data_dir,
  metadata = meta_file,
  filename_col = "Filename",
  batch_ids = "batch",
  condition = "condition",
  sample_ids = "Patient_id",
  transform = FALSE,
  down_sample = FALSE,
  # sample_size = 30000,
  # sampling_type = "sample_ids",
  panel = panel_file,
  panel_channel = "Channel",
  panel_antigen = "Marker"
  # markers = markers
)


detect_batch_effect_express(
  uncorrected,
  downsample = 40000,
  out_dir = meta_dir )


labels <- uncorrected %>%
  normalize(markers = markers,
            norm_method = "scale",
            ties.method = "average") %>%
  create_som(markers = markers,
             rlen = 30,
             seed = 353,
             xdim = 8,
             ydim = 8 )


corrected <- correct_data(
  df = uncorrected,
  label = labels,
  covar = NULL,
  markers = markers,
  parametric = TRUE )


saveRDS(corrected, file.path(data_dir, "corrected.RDS"))


# Re-run clustering on corrected data
labels <- corrected %>%
  create_som(markers = markers,
             rlen = 30,
             seed = 719,
             xdim = 8,
             ydim = 8 )
uncorrected$label <- corrected$label <- labels


# Evaluate EMD
emd <- evaluate_emd(uncorrected, corrected, cell_col = "label")

# # Reduction
# emd$reduction

# Violin & Scatter plots
viol_plot <- emd$violin
scat_plot <- emd$scatter
plot_save_two(viol_plot, scat_plot, file.path(meta_dir, "comparison_plots.png"))


meta_dir_corrected <- file.path(meta_dir, "corrected")
if (!file.exists(meta_dir_corrected)) {
  dir.create(meta_dir_corrected)
}

# # Create UMAPs
# sam <- sample(1:nrow(uncorrected), 40000)
# plot1 <- plot_dimred(uncorrected[sam, ], "Uncorrected", type = "umap", plot = "batch", markers = markers)
# plot2 <- plot_dimred(corrected[sam, ], "Corrected", type = "umap", plot = "batch", markers = markers)
# plot_save_two(plot1, plot2, file.path(meta_dir,"umap.png"))

detect_batch_effect_express(
  corrected,
  downsample = 40000,
  out_dir = file.path(meta_dir, "corrected") )
