library(ggplot2)
library(MASS)
library(glmnet)

set.seed(123)


n <- 100
B <- 2000   

pvals <- numeric(B)
beta <- c(3,1.5,0,0,2,0,0,0)
beta_0 <- which(beta == 0)

p <- length(beta)
b <- 0
while(b <= B){
  Sigma <- diag(p)
  X <- MASS::mvrnorm(n = n,rep(0,p),Sigma)
  epsilon <- rnorm(n)
  y = X%*%beta + epsilon
  
  lasso_model <- cv.glmnet(X, y, alpha = 1)
  beta_hat <- as.numeric(coef(lasso_model))[-1]

  lasso_selected_coefs = which(beta_hat != 0)
  coefs_possible <- intersect(beta_0,lasso_selected_coefs)
  coefs_possible
  if(length(coefs_possible) >= 1){
    j <- sample(coefs_possible, 1)
    resid <- y - X %*% beta_hat
    sigma_hat <- sqrt(mean(resid^2))
    
    se_j <- sigma_hat / sqrt(sum(X[, j]^2))
    z_j <- beta_hat[j] / se_j
    
    pvals[b] <- 2 * (1 - pnorm(abs(z_j)))
    b <- b+ 1
    print(b)

  }
  
}

pval_df <- data.frame(pval = pvals)  

lasso_pval<-ggplot(pval_df, aes(x = pval)) +
  stat_ecdf(geom = "step", color = "blue",size = 0.8,pad = FALSE) +
  geom_abline(intercept = 0, slope = 1,size=0.1,linetype="dashed") +
  labs(x = "p-value", y = "ECDF") +
  theme_minimal() +
  ggtitle("") +
  coord_cartesian(xlim = c(0, 1),ylim=c(0,1))+
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
ggsave(filename="lasso_pvalues.png",plot=lasso_pval,dpi=300,width=6,height=4,units="in")

# regularization path

Sigma <- diag(p)
X <- MASS::mvrnorm(n = n,rep(0,p),Sigma)
epsilon <- rnorm(n)
y = X%*%beta + epsilon

lasso_model <- cv.glmnet(X, y, alpha = 1)
path_model <- lasso_model$glmnet.fit

png("lasso_path.png", width = 6, height = 4, units = "in", res = 300, pointsize = 12)
plot(path_model, xvar = "lambda")
dev.off()
