library(ggplot2)
x <- seq(0, 1, length.out = 100)
df <- data.frame(
  x = x,
  y_super = x^2,   
  y_sub = sqrt(x)  
)
theme_clean_grid <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    axis.title.y = element_text(size = 28)  
  )
p1 <- ggplot(df, aes(x = x)) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey60") +
  geom_line(aes(y = y_super), color = "#0000b4", size = 1) +
  labs(y = "CDF") +
  theme_clean_grid +
  coord_fixed() 
p2 <- ggplot(df, aes(x = x)) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey60") +
  geom_line(aes(y = y_sub), color = "#0000b4", size = 1) +
  labs(y = "CDF") +
  theme_clean_grid +
  coord_fixed()
p1
p2
ggsave("super_uniform.png", plot = p1, width = 10, height = 5, dpi = 300)
ggsave("sous_uniform.png", plot = p2, width = 10, height = 5, dpi = 300)
