library(rhdf5)
library(ggplot2)
library(dplyr)

file <- "fulldata.h5ad"

# Get raw data of umap coordinate
umap_data <- h5read(file, "/obsm/X_umap")

age        <- h5read(file, "/obs/age_or_mean_of_age_range")
condition   <- h5read(file, "/obs/lung_condition")
condition   <- condition$categories[condition$codes + 1]
ann_level_1 <- h5read(file, "/obs/ann_level_1")
ann_level_1 <- ann_level_1$categories[ann_level_1$codes + 1]

# Build Dataframe
umap_dataframe <- data.frame(
  umap_xvar      = umap_data[1, ],
  umap_yvar      = umap_data[2, ],
  age        = age,
  condition  = condition,
  ann_level_1 = ann_level_1
)

# Classify based on age
umap_dataframe$age_group <- cut(
  umap_dataframe$age,
  breaks = c(5, 40, 65, Inf),
  labels = c("Young (5-40)", "Middle-aged (41-65)", "Older (>65)"),
  right  = TRUE
)

# Randomly choose 100k samples for plotting
set.seed(1234)
umap_sample <- umap_dataframe %>%
  filter(!is.na(age_group)) %>%
  sample_n(100000)

age_colors <- c("Young (5-40)"="#90BE6D","Middle-aged (41-65)"="#F9C74F","Older (>65)"="#F94144")
# Plotting based on age (UMAP Plot)
umap_plot_age <- ggplot(umap_sample, aes(x = umap_xvar, y = umap_yvar, color = age_group)) +
  geom_point(size = 0.05) +
  scale_color_manual(values = age_colors) +
  guides(color = guide_legend(override.aes = list(size=2, alpha=1))) +
  labs(
    title = "UMAP: Cell Distribution by Age Group",
    subtitle = "Subsampled to 100,000 cells for visualization",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Age Group"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "grey40"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold")
  )

print(umap_plot_age)
ggsave("Visualization_umap_age5.png", plot = umap_plot_age, width = 8, height = 6, dpi = 300, bg = "#F2F0F0")


cell_type_colors <- c("Endothelial"="#F54927","Epithelial"="#CEF525","Immune"="#18F25C","Stroma"="#1E6DF7","Unknown"="#F71EB6")
centroids_celltype <- umap_sample %>%
  group_by(ann_level_1) %>%
  summarise(umap_xvar = median(umap_xvar),umap_yvar = median(umap_yvar),.groups = "drop"
)

# UMAP colored by cell type
umap_plot_celltype <- ggplot(umap_sample, aes(x = umap_xvar, y = umap_yvar, color = ann_level_1)) +
  geom_point(size = 0.05) +
  scale_color_manual(values = cell_type_colors) +
  guides(color = guide_legend(override.aes = list(size=2, alpha=1))) +
  labs(
    title = "UMAP: Cell Distribution by Cell Type",
    subtitle = "Subsampled to 100,000 cells for visualization",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cell Type"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "grey40"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold")
  )

print(umap_plot_celltype)
ggsave("Visualization_umap_celltype5.png", plot = umap_plot_celltype, width = 8, height = 6, dpi = 300, bg = "#F2F0F0")