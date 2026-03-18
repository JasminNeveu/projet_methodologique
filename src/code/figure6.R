library(ggplot2)
library(MASS)
require(clusterpval)
require(fastcluster)
library(KmeansInference)

set.seed(123)

n <- 100      
B <- 2000    
p <- 5

rho <- seq(0.1, 1, by = 0.2)
i <- 1:p
j <- 1:p

# HAC
pvals <- matrix(NA_real_, nrow = B, ncol = length(rho))

for (k in seq_along(rho)) {
  
  r <- rho[k]
  R <- r^(abs(outer(i, j, "-")))
  
  for (b in 1:B) {
    X <- MASS::mvrnorm(n = n, rep(0, p), R)
    hcl <- hclust(dist(X, method="euclidean")^2,
                  method="average")
    clust_pair <- sample(c(1,2,3), 2)
    pvals[b, k] <-
      test_hier_clusters_exact(
        X,
        link="average",
        K=3,
        k1=clust_pair[1],
        k2=clust_pair[2],
        hcl=hcl,
        sig=1,
        iso=TRUE
      )$pval
  }
}


pval_df <- data.frame(
  pval = as.vector(pvals),
  rho  = factor(rep(rho, each = B))
)

colors <- c("#afb9c4", "#778da9", "#415a77", "#1b263b", "#0d1b2a")

hac_pval_grad_rho <-
  ggplot(pval_df, aes(x = pval, color = rho, group = rho)) + 
  stat_ecdf(geom = "step", size = 0.8, pad = FALSE) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed",
              size = 0.3) +
  labs(x = "p-value",
       y = "ECDF",
       color = expression(rho)) +
  scale_color_manual(values = colors) +   
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_minimal() +
  ggtitle("HAC - average linkage") +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
hac_pval_grad_rho
ggsave(filename="hac_pvalues_grad_rho.png",plot=hac_pval_grad_rho,dpi=300,width=6,height=4,units="in")


# kmeans

pvals <- matrix(NA_real_, nrow = B, ncol = length(rho))

for (k in seq_along(rho)) {
  r <- rho[k]
  R <- r^(abs(outer(i, j, "-")))
  for (b in 1:B) {
    X <- MASS::mvrnorm(n = n, rep(0, p), R)
    hcl <- hclust(dist(X, method="euclidean")^2,
                  method="average")
    clust_pair <- sample(c(1,2,3), 2)
    cl_1_2_inference_demo <- kmeans_inference(X, k=3, clust_pair[1], clust_pair[2],
                                              sig=1, iter.max = 20,seed=b)
    pvals[b,k] <- cl_1_2_inference_demo$pval
    
  }
}


pval_df <- data.frame(
  pval = as.vector(pvals),
  rho  = factor(rep(rho, each = B))
)

colors <- c("#afb9c4", "#778da9", "#415a77", "#1b263b", "#0d1b2a")

kmeans_pval_grad_rho <-
  ggplot(pval_df, aes(x = pval, color = rho, group = rho)) + 
  stat_ecdf(geom = "step", size = 0.8, pad = FALSE) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed",
              size = 0.3) +
  labs(x = "p-value",
       y = "ECDF",
       color = expression(rho)) +
  scale_color_manual(values = colors) +   
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_minimal() +
  ggtitle("k-means") +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
kmeans_pval_grad_rho
ggsave(filename="kmeans_pvalues_grad_rho.png",plot=kmeans_pval_grad_rho,dpi=300,width=6,height=4,units="in")


