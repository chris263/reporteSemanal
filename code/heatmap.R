#heatmap

# The mtcars dataset:
data <- as.matrix(mtcars)

# Default Heatmap
heatmap(data)

heatmap(data, Colv = NA, Rowv = NA, scale="column")




#ggplot

library(ggplot2)
library(hrbrthemes)

# Dummy data
x <- LETTERS[1:20]
y <- paste0("var", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 9)

# Heatmap
ggplot(data, aes(X, Y, fill= Z)) +
  geom_tile()

# Color Brewer palette

library(viridis)
myDF2 <- myDF %>% filter(local %in% dfLocal[31])
myDF2 <- myDF %>% filter(local %in% "22WNUBYG501BU22")
data <- data.frame(X=myDF2$range,Y=myDF2$row,Nota=myDF2$Nota, genotipo=myDF2$genotipo)

ggplot(data, aes(X, Y, fill= Nota)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme_ipsum()



data %>%
  group_by(genotipo) %>%
  # mutate(rescale = scales::rescale(Z)) %>%
  ggplot(., aes(x = factor(Y), y = genotipo)) +
  geom_tile(aes(alpha = X, fill = Z), color = "white") +
  scale_alpha(range = c(0.1, 1)) +
  theme(legend.position = "none")


