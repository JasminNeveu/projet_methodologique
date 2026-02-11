library(ggplot2)
library(MASS)

set.seed(123)

n <- 50       
B <- 2000       
p <- 5
pvals <- numeric(B)
Sigma <- diag(p)

#kmeans
for (b in 1:B) {
  x <- MASS::mvrnorm(n = n,rep(0,p),Sigma)
  km <- kmeans(x, centers = 3, nstart = 50)
  clusters <-km$cluster
  clusters
  
  clus_pair <- sample(c(1,2,3), 2)
  x1 <- x[clusters == clus_pair[1],]
  x2 <- x[clusters == clus_pair[2],]
  
  pvals[b] <- t.test(x1, x2)$p.value
}

pval_df <- data.frame(pval = pvals)  
kmean_pval <- ggplot(pval_df, aes(x = pval)) +
  stat_ecdf(geom = "step", color = "blue",size=0.8,pad = FALSE) +
  geom_abline(intercept = 0, slope = 1,size=0.1,linetype="dashed") +
  labs(x = "p-value", y = "ECDF") +
  theme_minimal() +
  ggtitle("k-means") +
  coord_cartesian(xlim = c(0, 1),ylim=c(0,1))

ggsave(filename="kmeans_pvalues.png",plot=kmean_pval,dpi=300,width=6,height=4,units="in")


# HAC

pvals <- numeric(B)

for (b in 1:B) {
  x <- MASS::mvrnorm(n = n,rep(0,p),Sigma)
  d <- dist(x, method = "euclidean")  
  hc <- hclust(d, method = "average")
  clusters <- cutree(hc, k = 3)
  clusters
  
  clus_pair <- sample(c(1,2,3), 2)
  x1 <- x[clusters == clus_pair[1],]
  x2 <- x[clusters == clus_pair[2],]
  
  pvals[b] <- t.test(x1, x2)$p.value
}

pval_df <- data.frame(pval = pvals)  
hac_pval<-ggplot(pval_df, aes(x = pval)) +
  stat_ecdf(geom = "step", color = "blue",size=0.8,pad = FALSE) +
  geom_abline(intercept = 0, slope = 1,size=0.1,linetype="dashed") +
  labs(x = "p-value", y = "ECDF") +
  theme_minimal() +
  ggtitle("HAC - average linkage") +
  coord_cartesian(xlim = c(0, 1),ylim=c(0,1))
ggsave(filename="hac_pvalues.png",plot=hac_pval,dpi=300,width=6,height=4,units="in")

