library(ggplot2)
library(dplyr)
library(fastcluster)
require(clusterpval)
library(reshape2)
devtools::install_github("lucylgao/clusterpval")
data <- readRDS("data/1KGP_100PC.Rda")
df <- readRDS("data/df.rds")

# HAC
data_hcl <- data %>% select(-pop,-ID)
data_hcl <- as.matrix(data_hcl)
len <- length(unique(data$pop))
hcl <- hclust(dist(data_hcl, method="euclidean"), method="ward.D") 
labels <- cutree(hcl, k = len)

pvals <- matrix(NA_real_, nrow = len, ncol = len)
n_tests <- len * (len - 1) / 2
pvec   <- numeric(n_tests)
idx <- 1
for(i in 1:(len-1)){
  for(j in (i+1):len){
    p <- test_hier_clusters_exact(data_hcl[,1:10], link="ward.D", K=len, k1=i, k2=j, hcl=hcl)$pval
    pvals[i,j] <- p
    pvec[idx] <-p
    idx <- idx + 1
    print(paste(i, "-", j))
  }
}

pvals
res <- pvals
res[lower.tri(res, diag = TRUE)] <- NA
df <- melt(res, varnames = c("k1", "k2"), value.name = "pval")
df <- na.omit(df)
df$k1 <- factor(df$k1, levels = 1:len)
df$k2 <- factor(df$k2, levels = 1:len)



pop_list <- list(
  EUR = c(1, 2,10),
  EAS = c(3, 4, 6, 22),
  AMR = c(5, 7, 9, 11,24),
  SAS = c(12, 16, 17, 18, 26),
  AFR = c(8, 13, 14, 15, 19, 20, 21, 23,25)
)

cluster_mapping <- data.frame(
  cluster = unlist(pop_list),
  super_pop = rep(names(pop_list), lengths(pop_list))
)
ordered_clusters <- unlist(lapply(pop_list, sort))

cluster_mapping <- data.frame(
  cluster = ordered_clusters,
  super_pop = rep(names(pop_list), lengths(pop_list))
)
df_ordered <- bind_rows(
  df, 
  data.frame(k1 = df$k2, k2 = df$k1, pval = df$pval)
) %>%
  distinct(k1, k2, .keep_all = TRUE) %>%
  mutate(
    k1 = factor(k1, levels = ordered_clusters),
    k2 = factor(k2, levels = ordered_clusters)
  ) %>%
  filter(as.integer(k1) >= as.integer(k2))

pop_bounds <- cluster_mapping %>%
  mutate(position = row_number()) %>%
  group_by(super_pop) %>%
  summarise(start = min(position) - 0.5, end = max(position) + 0.5, mid = mean(position))



heatmap_pval <- ggplot(df_ordered, aes(x = k1, y = k2, fill = pval)) +
  geom_tile(color = "black", linewidth = 0.3) +
  scale_fill_gradient(low = "#ffffff", high = "#1b263b", name = "p-value", na.value = "white") +
  coord_equal(clip = "off") + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust =0.5),
    plot.margin = margin(t = 10, r = 10, b = 50, l = 50) 
  )

for(i in 1:nrow(pop_bounds)) {
  heatmap_pval <- heatmap_pval +
    annotate("rect", xmin = pop_bounds$start[i], xmax = pop_bounds$end[i], 
             ymin = -2, ymax = -0.7, fill = "#fff", color = "black") +
    annotate("text", x = pop_bounds$mid[i], y = -1.35, 
             label = pop_bounds$super_pop[i], size = 3) +
    annotate("rect", ymin = pop_bounds$start[i], ymax = pop_bounds$end[i], 
             xmin = -2, xmax = -0.7, fill = "#fff", color = "black") +
    annotate("text", y = pop_bounds$mid[i], x = -1.35, 
             label = pop_bounds$super_pop[i], size = 3, angle = 90)
}


ggsave(filename="heatmap_pval.png",plot=heatmap_pval,dpi=300,width=10,height=6,units="in")

