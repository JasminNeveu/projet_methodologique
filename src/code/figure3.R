library(MASS)
library(ggplot2)
library(latex2exp)

n  <- 10
p  <- 2
Sigma <- diag(p)/2
phi <- 8

mu1 <- c(0,  2.5)
mu2 <- c( 4,  0)
mu3 <- c( 0,  -2.5)

X1 <- mvrnorm(n, mu1, Sigma)
X2 <- mvrnorm(n, mu2, Sigma)
X3 <- mvrnorm(n, mu3, Sigma)
X <- rbind(X1, X2, X3)   

nu <- c(rep(1/n, n), rep(-1/n, n), rep(0, n))
nu_norm2_sq <- sum(nu^2)
proj   <- t(nu) %*% X
proj_norm <- sqrt(sum(proj^2))                  
dir_proj  <- proj / proj_norm 
proj_norm
X_pert <- X + outer(nu, as.vector(dir_proj)) / nu_norm2_sq * (phi - proj_norm)

cluster <- factor(rep(1:3, each = n))
df <- data.frame(
  x       = X_pert[, 1],
  y       = X_pert[, 2],
  cluster = cluster
)

perturbated_data_phi <- ggplot(df, aes(x, y, color = cluster)) +
  geom_point(size = 2) +
  labs(
    title = TeX(r"(a) Original data)"),
    x = "Feature 1",
    y = "Feature 2"
  ) +
  coord_cartesian(xlim = c(-2.5, 7.5), ylim = c(-4.5, 4.5)) +
  theme_classic() + 
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
  )
perturbated_data_phi
ggsave(filename="original_data.png",plot=perturbated_data_phi,dpi=300,width=6,height=4,units="in")

