library(ggplot2)
library(MASS)
require(clusterpval)
require(fastcluster)
library(KmeansInference)

set.seed(123)

n <- 100      
B <- 2000       
p <- 5
pvals <- numeric(B)
Sigma <- diag(p)
i <- 1:p
j <- 1:p
R <- 0.5^(abs(outer(i, j, "-")))

# HAC

pvals <- numeric(B)

for (b in 1:B) {
  X <- MASS::mvrnorm(n = n,rep(0,p),R)
  hcl <- hclust(dist(X, method="euclidean")^2, method="average") 
  clust_pair <- sample(c(1,2,3), 2)
  pvals[b] <-  test_hier_clusters_exact(X, link="average", K=3, k1=clust_pair[1], k2=clust_pair[2], hcl=hcl,sig=1,iso=TRUE)$pval
}

pval_df <- data.frame(pval = pvals)  
hac_pval_not_indep<-ggplot(pval_df, aes(x = pval)) +
  stat_ecdf(geom = "step", color = "blue",size=0.8,pad = FALSE) +
  geom_abline(intercept = 0, slope = 1,size=0.1,linetype="dashed") +
  labs(x = "p-value", y = "ECDF") +
  theme_minimal() +
  ggtitle("HAC - average linkage") +
  coord_cartesian(xlim = c(0, 1),ylim=c(0,1)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave(filename="hac_pvalues_not_indep.png",plot=hac_pval_not_indep,dpi=300,width=6,height=4,units="in")




#kmeans
pvals <- numeric(B)

for (b in 1:B) {
  X <- MASS::mvrnorm(n = n,rep(0,p),R)
  clust_pair <- sample(c(1,2,3), 2)
  cl_1_2_inference_demo <- kmeans_inference(X, k=3, clust_pair[1], clust_pair[2],
                                            sig=1, iter.max = 20,seed=b)
  pvals[b] <- cl_1_2_inference_demo$pval
}

pval_df <- data.frame(pval = pvals)
kmean_pval_not_indep <- ggplot(pval_df, aes(x = pval,group=1)) +
  stat_ecdf(geom = "step", color = "blue", size = 0.8, pad = FALSE) +
  geom_abline(intercept = 0, slope = 1, size = 0.1, linetype = "dashed") +
  labs(x = "p-value", y = "ECDF", title = "k-means") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave(filename="kmeans_pvalues_not_indep.png",plot=kmean_pval_not_indep,dpi=300,width=6,height=4,units="in")

