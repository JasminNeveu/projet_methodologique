library(dplyr)
library(fastcluster)
require(clusterpval)
library(ggplot2)
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
  geom_vline(xintercept = 5,linetype="dashed",linewidth=0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.5) +
  annotate(
    "text",
    x = Inf,
    y = 0.05,
    label = "0.05",
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
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
  pval <- test_hier_clusters_exact(data_hcl[,1:p], link="ward.D", K=len, k1=8, k=13, hcl=hcl)$pval
  results$pval[i] <-pval
  print(p)
}


plot_pval_8_13<-ggplot(results, aes(x = p_list, y = pval)) +
  geom_line(color = "#1b263b", linewidth = 1) +
  geom_point(color = "#1b263b", size = 3) +
  geom_vline(xintercept = 5,linetype="dashed",linewidth=0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.5) +
  annotate(
    "text",
    x = Inf,
    y = 0.05,
    label = "0.05",
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  scale_y_log10() + 
  labs(
    x = "p",
    y = "p-value"
  ) +
  theme_minimal()
ggsave(filename="pval_8_13.png",plot=plot_pval_8_13,dpi=300,width=6,height=4,units="in")




results <- data.frame(n_pcs = p_list, pval = NA)
for(i in seq_along(p_list)){
  p <- p_list[i]
  pval <- test_hier_clusters_exact(data_hcl[,1:p], link="ward.D", K=len, k1=5, k=7, hcl=hcl)$pval
  results$pval[i] <-pval
  print(p)
}


plot_pval_5_7<-ggplot(results, aes(x = p_list, y = pval)) +
  geom_line(color = "#1b263b", linewidth = 1) +
  geom_point(color = "#1b263b", size = 3) +
  geom_vline(xintercept = 5,linetype="dashed",linewidth=0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.5) +
  annotate(
    "text",
    x = Inf,
    y = 0.05,
    label = "0.05",
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  scale_y_log10() + 
  labs(
    x = "p",
    y = "p-value"
  ) +
  theme_minimal()
plot_pval_5_7
ggsave(filename="pval_5_7.png",plot=plot_pval_5_7,dpi=300,width=6,height=4,units="in")





results <- data.frame(n_pcs = p_list, pval = NA)
for(i in seq_along(p_list)){
  p <- p_list[i]
  pval <- test_hier_clusters_exact(data_hcl[,1:p], link="ward.D", K=len, k1=9, k=16, hcl=hcl)$pval
  results$pval[i] <-pval
  print(p)
}

plot_pval_9_16<-ggplot(results, aes(x = p_list, y = pval)) +
  geom_line(color = "#1b263b", linewidth = 1) +
  geom_point(color = "#1b263b", size = 3) +
  geom_vline(xintercept = 5,linetype="dashed",linewidth=0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.5) +
  annotate(
    "text",
    x = Inf,
    y = 0.05,
    label = "0.05",
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  scale_y_log10() + 
  labs(
    x = "p",
    y = "p-value"
  ) +
  theme_minimal()
plot_pval_9_16
ggsave(filename="pval_9_16.png",plot=plot_pval_9_16,dpi=300,width=6,height=4,units="in")

