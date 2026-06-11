library(rhdf5)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(patchwork)

# Read through the file
file <- "fulldata.h5ad"

age         <- h5read(file, "/obs/age_or_mean_of_age_range")
condition   <- h5read(file, "/obs/lung_condition")
condition   <- condition$categories[condition$codes + 1]
ann_level_1 <- h5read(file, "/obs/ann_level_1")
ann_level_1 <- ann_level_1$categories[ann_level_1$codes + 1]
cell_type   <- h5read(file, "/obs/cell_type")
cell_type   <- cell_type$categories[cell_type$codes + 1]
tissue      <- h5read(file, "/obs/tissue")
tissue      <- tissue$categories[tissue$codes + 1]
sex         <- h5read(file, "/obs/sex")
sex         <- sex$categories[sex$codes + 1]

# Build Dataframe
meta <- data.frame(
  age         = age,
  condition   = condition,
  ann_level_1 = ann_level_1,
  cell_type   = cell_type,
  tissue      = tissue,
  sex         = sex,
  stringsAsFactors = FALSE
)

# Classify based on age
meta$age_group <- cut(
  meta$age,
  breaks = c(5, 40, 65, Inf),
  labels = c("Young (5-40)", "Middle-aged (41-65)", "Older (>65)"),
  right  = TRUE
)

# Check the counts of each group
cat("Age groups\n")
print(table(meta$age_group, useNA = "always"))
cat("Age groups X Ann_level1\n")
# Clean the data with no age
meta_aged <- meta %>% filter(!is.na(age_group))
print(table(meta_aged$ann_level_1, meta_aged$age_group))
cat("Age groups X Lung Condition\n")
print(table(meta_aged$condition, meta_aged$age_group))


# Color Setting
age_colors <- c("Young (5-40)"="#90BE6D","Middle-aged (41-65)"="#F9C74F","Older (>65)"="#F94144")
cell_colors <- c("Endothelial"="#4E9BCD","Epithelial"="#F28E2B","Immune"="#E15759","Stroma"="#76B7B2")
cell_colors_with_unknown <- c("Endothelial"="#4E9BCD","Epithelial"="#F28E2B","Immune"="#E15759","Stroma"="#76B7B2", "Unknown"="#B0B0B0")

# Data Visualization

# Plot1: Cell type composition for each age
# Data
plot1_data <- meta_aged %>%
  group_by(age_group, ann_level_1) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(age_group) %>%
  mutate(percentage = n / sum(n) * 100)
# Plotting
plot1 <- ggplot(plot1_data, aes(x = age_group, y = percentage, fill = ann_level_1)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(
    aes(label = ifelse(percentage >= 1, sprintf("%.1f%%", percentage), "")),
    position = position_stack(vjust = 0.5),
    size = 5, fontface = "bold", color = "white"
  ) +
  scale_fill_manual(values = cell_colors_with_unknown) +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  labs(
    title = "Plot 1: Composition of Lung Cell Type Across Age Groups",
    subtitle = "Proportion of Ann Level 1 cell shifts with age",
    x = "Age Group",
    y = "Percentage of Cells (%)",
    fill = "Cell Lineage"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "grey60"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold")
  )

print(plot1)
ggsave("Visualization1_new.png", plot = plot1, width = 8, height = 6, dpi = 300, bg = "#F2F0F0")


# Plot2: Age Data Coverage
# Data
all_conditions_cell_with_age <- meta %>%
  group_by(condition) %>%
  summarise(total_cells = n(), cells_with_age = sum(!is.na(age)), cells_without_age = sum(is.na(age)), percentage_has_age = round(sum(!is.na(age))/n()* 100, 1), .groups = "drop") %>%
  mutate(has_age = ifelse(cells_with_age > 0, "YES", "NO")) %>%
  arrange(has_age, desc(total_cells))
print(all_conditions_cell_with_age, n = Inf)

conditions_with_age <- c("Healthy", "COVID-19", "Healthy (tumor adjacent)",
                         "Chronic rhinitis", "Systemic sclerosis-associated ILD",
                         "Lymphangioleiomyomatosis", "End-stage lung fibrosis, unknown etiology")

plot2_data <- all_conditions_cell_with_age %>%
  arrange(desc(condition)) %>%
  mutate(condition = factor(condition, levels = condition))
plot2_axis_colors <- ifelse(levels(plot2_data$condition) %in% plot2_data$condition[plot2_data$has_age == "YES"],"black", "#FF0000")

# Plotting
plot2 <- ggplot(plot2_data, aes(x = condition, y = percentage_has_age, fill = has_age)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = paste0(percentage_has_age, "%")), hjust = -0.15, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("YES" = "#90BE6D", "NO" = "#FF0000")) +
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  coord_flip() +
  labs(
    title = "Plot2: Age Data Percentage by Lung Condition",
    subtitle = "Black labels = Has Age Data,   Red labels = No Age Data",
    x = "Lung Condition",
    y = "% of Cells with Age Data",
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "grey60"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10, color = plot2_axis_colors),
    axis.title = element_text(size = 14),
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.5, linetype = "dashed")
  )
print(plot2)
ggsave("Visualization2_new.png", plot = plot2, width = 12, height = 8, dpi = 300, bg = "#F2F0F0")



# Plot3: Cell counts for each condition for each group
# Data
plot3_data <- meta_aged %>%
  # Drop all conditions with no age data
  filter(condition %in% conditions_with_age) %>%
  group_by(condition, age_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  # Sort by total cell counts
  group_by(condition) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  arrange(desc(total)) %>%
  mutate(condition = factor(condition, levels = unique(condition)))

# Plotting
plot3 <- ggplot(plot3_data, aes(x = condition, y = n, fill = age_group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  geom_text(aes(label = comma(n)), position = position_dodge(width = 0.75), vjust = -0.4, size = 3, fontface = "bold") +
  scale_fill_manual(values = age_colors) +
  scale_y_continuous(label = comma, limits = c(0,550000), expand = c(0,0)) +
  labs(
    title = "Plot 3: Cell Counts by Lung Condition and Age Group",
    subtitle = "Shown conditions with age data from Plot 2",
    x = "Lung Condition",
    y = "Number of Cells",
    fill = "Age Group"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "grey60"),
    axis.text.x = element_text(angle = 20, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5, linetype = "dashed")
  )

print(plot3)
ggsave("Visualization3_new.png", plot = plot3, width = 12, height = 8, dpi = 300, bg = "#F2F0F0")



# Plot4: Cell composition for each condition for each group
# Data
plot4_data <- meta_aged %>%
  filter(condition %in% conditions_with_age, ann_level_1 != "Unknown") %>%
  group_by(condition, age_group, ann_level_1) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition, age_group) %>%
  mutate(percentage = n / sum(n) * 100)

plot4_conditions_list <- unique(plot4_data$condition)

# Plot a single graph function
make_single_plot <- function(cond) {
  current_data <- plot4_data %>% filter(condition == cond)
  ggplot(current_data, aes(x = age_group, y = percentage, fill = ann_level_1)) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_text(aes(label = ifelse(percentage>=5, sprintf("%.1f%%", percentage),"")), position = position_stack(vjust = 0.5), size = 3, fontface = "bold", color = "white") +
    scale_fill_manual(values = cell_colors) +
    scale_y_continuous(limits = c(0,105), expand = c(0,0)) +
    labs(title = cond, x = "Age Group", y = "Percentage of Cells (%)", fill  = "Cell Lineage") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "none"
    )
}

# Plot the last graph with legend
make_legend_plot <- function(cond) {
  current_data <- plot4_data %>% filter(condition == cond)
  ggplot(current_data, aes(x = age_group, y = percentage, fill = ann_level_1)) +
    geom_bar(stat = "identity", width = 0.75) +
    geom_text(aes(label = ifelse(percentage>=5, sprintf("%.1f%%", percentage),"")), position = position_stack(vjust = 0.5), size = 3, fontface = "bold", color = "white") +
    scale_fill_manual(values = cell_colors) +
    scale_y_continuous(limits = c(0,105), expand = c(0,0)) +
    labs(title = cond, x = "Age Group", y = "Percentage of Cells (%)", fill  = "Cell Lineage") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "right",
      legend.text     = element_text(size = 9),
      legend.title    = element_text(size = 11, face = "bold")
    )
}

# Combine 7 graphs
plots <- lapply(plot4_conditions_list, make_single_plot)
plots[[7]] <- make_legend_plot(plot4_conditions_list[7])

plot4 <- (plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) / (plots[[5]] | plots[[6]] | plots[[7]] | plot_spacer()) +
  plot_annotation(
    title = "Plot 4: Bar chart of Cell Lineage Proportions",
    subtitle = "Comparing cell lineage proportions across all conditions",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "grey60")
    )
  )

print(plot4)
ggsave("Visualization4_new.png", plot = plot4, width = 12, height = 8, dpi = 300, bg = "#F2F0F0")


# Plot4.2 Heat map
# Data
plot4_data_2 <- meta_aged %>%
  filter(condition %in% conditions_with_age,ann_level_1 != "Unknown") %>%
  group_by(condition, age_group, ann_level_1) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition, age_group) %>%
  mutate(pct = n / sum(n) * 100) %>%
  unite("group", condition, age_group, sep = "\n")

plot4_2 <- ggplot(plot4_data_2, aes(x = group, y = ann_level_1, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), size = 4, fontface = "bold") +
  scale_fill_gradient(low = "#FFFFFF", high = "#FF0000") +
  labs(
    title = "Plot 4: Heatmap of Cell Lineage Proportions",
    subtitle = "Percentage of each cell lineage per condition and age group",
    x = "Condition / Age Group",
    y = "Cell Lineage",
    fill = "% of Cells"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "grey60"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold")
  )

print(plot4_2)
ggsave("Visualization4_Version2_new.png", plot = plot4_2, width = 12, height = 8, dpi = 300, bg = "#F2F0F0")
