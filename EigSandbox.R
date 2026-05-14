
theta <- 0:359/360 * 2 * pi
x <- cos(theta)
y <- sin(theta)

A <- matrix(c(2,2,3,3.5), nrow = 2, byrow = TRUE)
X <- rbind(x, y)
AX <- A %*% X

eig    <- eigen(A)
evals  <- Re(eig$values)
evecs  <- Re(eig$vectors)   # columns are unit eigenvectors

# scale eigenvectors by their eigenvalues to get A*v = lambda*v
evecs_scaled <- sweep(evecs, 2, evals, `*`)

library(ggplot2)
library(gganimate)

colour_map  <- c("X" = "red", "AX" = "blue",
                 "eigvec 1" = "darkgreen", "A·eigvec 1" = "darkgreen",
                 "eigvec 2" = "orange",    "A·eigvec 2" = "orange")
linetype_map <- c("X" = "solid", "AX" = "solid",
                  "eigvec 1" = "solid",  "A·eigvec 1" = "solid",
                  "eigvec 2" = "solid",  "A·eigvec 2" = "solid")
linewidth_map <- c("X" = 0.6, "AX" = 0.6,
                   "eigvec 1" = 1.2, "A·eigvec 1" = 0.6,
                   "eigvec 2" = 1.2, "A·eigvec 2" = 0.6)

df_curves <- rbind(
  data.frame(x = X[1,],  y = X[2,],  series = "X",  theta = theta),
  data.frame(x = AX[1,], y = AX[2,], series = "AX", theta = theta)
)

df_arrows <- data.frame(
  x    = 0,
  y    = 0,
  xend = c(evecs[1,], evecs_scaled[1,]),
  yend = c(evecs[2,], evecs_scaled[2,]),
  series = c("eigvec 1", "eigvec 2", "A·eigvec 1", "A·eigvec 2")
)

df_labels <- data.frame(
  x     = evecs_scaled[1,],
  y     = evecs_scaled[2,],
  label = sprintf("λ%d = %.2f", 1:2, evals),
  series = c("A·eigvec 1", "A·eigvec 2")
)

p <- ggplot() +
  geom_path(data = df_curves,
            aes(x = x, y = y, colour = series, linetype = series, linewidth = series)) +
  geom_segment(data = df_arrows,
               aes(x = x, y = y, xend = xend, yend = yend,
                   colour = series, linetype = series, linewidth = series),
               arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("label", x = Inf, y = 0, hjust = 1.05, vjust = -1.5,
           label = sprintf("λ1 = %.2f", evals[1]), colour = "darkgreen", size = 3) +
  annotate("label", x = Inf, y = 0, hjust = 1.05, vjust = -0.1,
           label = sprintf("λ2 = %.2f", evals[2]), colour = "orange",    size = 3) +
  scale_colour_manual(values = colour_map,     name = NULL,
                      breaks = c("X", "AX", "eigvec 1", "A·eigvec 1", "eigvec 2", "A·eigvec 2")) +
  scale_linetype_manual(values = linetype_map, name = NULL,
                        breaks = c("X", "AX", "eigvec 1", "A·eigvec 1", "eigvec 2", "A·eigvec 2")) +
  scale_linewidth_manual(values = linewidth_map, name = NULL,
                         breaks = c("X", "AX", "eigvec 1", "A·eigvec 1", "eigvec 2", "A·eigvec 2")) +
  coord_fixed(clip = "off") +
  labs(x = "x", y = "y", title = "X and AX") +
  guides(colour = guide_legend(ncol = 3)) +
  theme_classic() +
  theme(legend.position = "bottom") +
  transition_reveal(theta)

anim <- animate(p, nframes = 180, fps = 30, width = 700, height = 700, renderer = gifski_renderer())
anim_save("EigSandbox_anim.gif", anim)

