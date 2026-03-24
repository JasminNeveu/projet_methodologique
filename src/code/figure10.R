library(dplyr)
library(fastcluster)
require(clusterpval)
super_pop <- read.table("data/20131219.populations.tsv",sep="\t",header = TRUE)
super_pop <- head(super_pop,-3)
unique(super_pop$Super.Population)

data <- readRDS("data/1KGP_100PC.Rda")
unique(data$pop)


data_hcl <- data %>% select(-pop,-ID)
data_hcl <- as.matrix(data_hcl)
len <- length(unique(data$pop))
hcl <- hclust(dist(data_hcl, method="euclidean"), method="ward.D") 
labels <- cutree(hcl, k = len)

labels

# curse of dimensionality sous H0

p_list <- c(2,5,15,20,25,50,100)
results <- data.frame(n_pcs = p_list, pval = NA)
for(i in seq_along(p_list)){
  p <- p_list[i]
  pval <- test_hier_clusters_exact(data_hcl[,1:p], link="ward.D", K=len, k1=12, k=26, hcl=hcl)$pval
  results$pval[i] <-pval
  print(p)
}


plot_pval_12_26<-ggplot(results, aes(x = p_list, y = pval)) +
  geom_line(color = "#1b263b", linewidth = 1) +
  geom_point(color = "#1b263b", size = 3) +
  scale_y_log10() + 
  labs(
    x = "p",
    y = "p-value"
  ) +
  theme_minimal()
ggsave(filename="pval_12_26.png",plot=plot_pval_12_26,dpi=300,width=6,height=4,units="in")


results <- data.frame(n_pcs = p_list, pval = NA)
for(i in seq_along(p_list)){
  p <- p_list[i]
  pval <- test_hier_clusters_exact(data_hcl[,1:p], link="ward.D", K=len, k1=13, k=8, hcl=hcl)$pval
  results$pval[i] <-pval
  print(p)
}


plot_pval_13_8<-ggplot(results, aes(x = p_list, y = pval)) +
  geom_line(color = "#1b263b", linewidth = 1) +
  geom_point(color = "#1b263b", size = 3) +
  scale_y_log10() + 
  labs(
    x = "p",
    y = "p-value"
  ) +
  theme_minimal()
ggsave(filename="pval_13_8.png",plot=plot_pval_13_8,dpi=300,width=6,height=4,units="in")

