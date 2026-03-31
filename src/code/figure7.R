library(ggplot2)
library(fastcluster)
library(dplyr)

data <- readRDS("data/1KGP_100PC.Rda")
super_pop <- read.table(
  "data/20131219.populations.tsv",
  sep = "\t",
  header = TRUE
)
super_pop <- head(super_pop, -3)
meta <- super_pop %>% 
  select(
    pop = Population.Code,
    superpop = Super.Population
  ) %>% 
  distinct()

data_hcl <- data %>% 
  select(-pop, -ID)
data_hcl <- as.matrix(data_hcl)

len <- length(unique(data$pop))
hcl <- hclust(
  dist(data_hcl, method = "euclidean"),
  method = "ward.D"
)
clusters <- cutree(hcl, k = len)

tab <- table(data$pop, clusters)
tab_prop <- prop.table(tab, margin = 1)
best_cluster <- apply(tab_prop, 1, which.max)
pop_order <- names(sort(best_cluster))

df <- as.data.frame(tab_prop)
colnames(df) <- c("Population", "cluster", "perc")

label_df <- data.frame(Population = pop_order) %>%
  left_join(meta, by = c("Population" = "pop")) %>%
  mutate(pop_label = Population)

df <- df %>%
  left_join(label_df, by = "Population")

df$pop_label <- factor(df$pop_label, levels = label_df$pop_label)
df$cluster <- factor(df$cluster, levels = sort(unique(clusters)))

superpop_colors <- c(
  "AFR" = "#ff9f1c",
  "AMR" = "#e71d36",
  "EAS" = "#2ec4b6",
  "EUR" = "#011627",
  "SAS" = "#984EA3"
)

label_colors <- label_df %>%
  mutate(color = superpop_colors[superpop]) %>%
  pull(color)

heatmap <- ggplot(
  df,
  aes(x = pop_label, y = cluster, fill = perc)
) +
  geom_tile(color = "black", linewidth = 0.3) +
  geom_point(aes(color = superpop), alpha = 0, size = 0) +
  labs(x = "Population", y = "HAC cluster") +
  scale_fill_gradient(
    low = "white",
    high = "#1b263b",
    name = "Proportion \nof 1KGP \npopulation"
  ) +
  scale_color_manual(
    name = "Super-population",
    values = superpop_colors,
    guide = guide_legend(override.aes = list(alpha = 1, size = 3))
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 7,
      colour = label_colors
    )
  )

heatmap

ggsave(
  filename = "heatmap_clust_pop.png",
  plot = heatmap,
  dpi = 300,
  width = 8,
  height = 5,
  units = "in"
)

