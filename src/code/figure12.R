library(dplyr)
library(ggplot2)
df <- readRDS("data/df.rds")
cluster_map <- data.frame(
  cluster = c(
    clusters_eur,
    clusters_eas,
    clusters_amr,
    clusters_sas,
    clusters_afr
  ),
  superpop = c(
    rep("EUR", length(clusters_eur)),
    rep("EAS", length(clusters_eas)),
    rep("AMR", length(clusters_amr)),
    rep("SAS", length(clusters_sas)),
    rep("AFR", length(clusters_afr))
  )
)

df2 <- df

df2$k1 <- as.numeric(as.character(df2$k1))
df2$k2 <- as.numeric(as.character(df2$k2))


df2 <- df2 %>% 
  left_join(
    cluster_map %>%  rename(super1 = superpop),
    by = c("k1" = "cluster")
  )

df2 <- df2 %>% 
  left_join(
    cluster_map %>%  rename(super2 = superpop),
    by = c("k2" = "cluster")
  )
df_between <- df2 %>% 
  filter(super1 != super2)

df_inner <- df2 %>% 
  filter(super1 == super2)

ecdfplotdiff <- ggplot(df_between, aes(x = pval)) +
  stat_ecdf(linewidth = 1) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  labs(x = "p-value",
       y = "ECDF"
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )


ecdfplotinner <- ggplot(df_inner, aes(x = pval)) +
  stat_ecdf(linewidth = 1) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  labs(x = "p-value",
       y = "ECDF"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
ecdfplotinner
ecdfplotdiff
ggsave(filename="ecdf_diff_super.png",plot=ecdfplotdiff,dpi=300,width=6,height=4,units="in")
ggsave(filename="ecdf_inner_super.png",plot=ecdfplotinner,dpi=300,width=6,height=4,units="in")
