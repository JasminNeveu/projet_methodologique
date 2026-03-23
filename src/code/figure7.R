library(ggplot2)
library(dplyr)
library(ggrepel)
library(fastcluster)
source("data/1KGP_colours.R")

data <- readRDS("data/1KGP_100PC.Rda")

named_colors <- setNames(
  sapply(colour_list, `[`, 2),
  sapply(colour_list, `[`, 1)
)

centroids <- data %>%
  group_by(pop) %>%
  summarise(x = mean(PC1) , y = mean(PC2))

# color pop
ACPplotpop<-ggplot(data,aes(x=PC1,y=PC2,color=pop)) +
  geom_point(size=0.7) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_blank()) +
  scale_color_manual(values = named_colors,guide="none") +
  geom_text_repel(data = centroids, aes(x = x, y = y, label = pop, color = pop),
             inherit.aes = FALSE,fontface="bold", segment.color= NA,max.overlaps=Inf,point.padding=6,box.padding=0.3)
ACPplotpop
ggsave(filename="PCA.png",plot=ACPplotpop,dpi=300,width=6,height=4,units="in")

#color cluster
data_hcl <- data %>% select(-pop,-ID)
len <- length(unique(data$pop))
hcl <- hclust(dist(data_hcl, method="euclidean"), method="ward.D") 
clusters <- cutree(hcl,len)
data <- data %>% mutate(cluster = factor(clusters))

ACPplotcluster<-ggplot(data,aes(x=PC1,y=PC2,color=cluster)) +
  geom_point(size=0.7) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_blank())
ACPplotcluster

ggsave(filename="PCAcluster.png",plot=ACPplotcluster,dpi=300,width=6,height=4,units="in")


