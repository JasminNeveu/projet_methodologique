library(ggplot2)
library(fastcluster)
library(dplyr)

data <- readRDS("data/1KGP_100PC.Rda")

data_hcl <- data %>% select(-pop,-ID)
data_hcl <- as.matrix(data_hcl)

len <- length(unique(data$pop))

hcl <- hclust(
  dist(data_hcl, method = "euclidian"),
  method = "ward.D"
)

clusters <- cutree(hcl, k = len)

tab <- table(data$pop, clusters)
tab_prop <- prop.table(tab, margin = 1)

best_cluster <- apply(tab_prop, 1, which.max)
pop_order <- names(sort(best_cluster))

df <- as.data.frame(tab_prop)
colnames(df) <- c("Population", "cluster", "perc")

df$Population <- factor(df$Population, levels = pop_order)
df$cluster <- factor(df$cluster, levels = sort(unique(clusters)))

heatmap<-ggplot(df, aes(x = Population, y = cluster, fill = perc)) +
  geom_tile(color = "black", linewidth = 0.3) +
  labs(x = "Population",
       y = "HAC cluster") +
  scale_fill_gradient(
    low = "white",
    high = "#1b263b",
    name = "Proportion \nof 1KGP \npopulation"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
heatmap
ggsave(filename="heatmap_clust_pop.png",plot=heatmap,dpi=300,width=6,height=4,units="in")

