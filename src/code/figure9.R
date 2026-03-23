library(ggplot2)
library(dplyr)
library(fastcluster)
require(clusterpval)
library(reshape2)
devtools::install_github("lucylgao/clusterpval")


data <- readRDS("1KGP_100PC.Rda")

# HAC
data_hcl <- data %>% select(-pop,-ID)
data_hcl <- as.matrix(data_hcl)
len <- length(unique(data$pop))
hcl <- hclust(dist(data_hcl, method="euclidean"), method="ward.D") 
pvals <- matrix(NA_real_, nrow = len, ncol = len)

for(i in 1:(len-1)){
  for(j in (i+1):len){
    pvals[i,j] <- test_hier_clusters_exact(data_hcl[,1:10], link="ward.D", K=len, k1=i, k2=j, hcl=hcl)$pval
    print(paste(i, "-", j))
  }
}

res <- pvals
res[lower.tri(res, diag = TRUE)] <- NA

df <- melt(res, varnames = c("k1", "k2"), value.name = "pval")
df <- na.omit(df)
df$k1 <- factor(df$k1, levels = 1:len)
df$k2 <- factor(df$k2, levels = 1:len)

heatmap_pval <- ggplot(df, aes(x = k1, y = k2, fill = pval)) +
  geom_tile(color = "black", linewidth = 0.3) +
  scale_fill_gradient(low = "#ffffff", high = "#1b263b", name = "p-value") +
  coord_equal() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 270, hjust = 1)
  )
heatmap_pval
ggsave(filename="heatmap_pval.png",plot=heatmap_pval,dpi=300,width=6,height=4,units="in")

