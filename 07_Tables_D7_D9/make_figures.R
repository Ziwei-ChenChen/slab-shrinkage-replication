# Figures D1 and D2 of Appendix D (cyber-sickness data). They do not depend
# on any estimator and are written to figures/ under the file names used in
# the paper.
#
#   Figure D1  initialclass_distribution  FMS scores (2.5 and 3.5 rounded down)
#              4class_distribution        the four severity classes
#   Figure D2  heatmap13features          correlations of the 13 features, M = 1
#              heatmap39of130features     correlations of the lags t-3, t-2, t-1
#
# Run from this folder:  Rscript make_figures.R

library(ggplot2)
library(corrplot)

d1 <- read.csv("data/data_for_regression1.csv", check.names = FALSE)
d10 <- read.csv("data/data_for_regression10.csv", check.names = FALSE)
fms <- d1[[ncol(d1)]]
dir.create("figures", showWarnings = FALSE)

save_both <- function(plot, name, width, height) {
  ggsave(file.path("figures", paste0(name, ".png")), plot = plot, dpi = 600, width = width, height = height, units = "in")
  setEPS()
  postscript(file.path("figures", paste0(name, ".eps")), width = width, height = height, family = "Times")
  print(plot)
  dev.off()
}

bar_plot <- function(counts, title, xlab, bar_width) {
  ggplot(counts, aes(x = Class, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue", width = bar_width) +
    labs(title = title, x = xlab, y = "Count") +
    theme_classic(base_size = 14, base_family = "Times") +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "plain", family = "Times"),
      axis.title = element_text(size = 18, face = "plain", family = "Times"),
      axis.text = element_text(size = 16, face = "plain", family = "Times"),
      strip.text = element_text(size = 16, face = "plain", family = "Times"),
      legend.position = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
}


# Figure D1 -------------------------------------------------------------------

counts <- data.frame(Class = factor(0:10), Count = tabulate(floor(fms) + 1, nbins = 11))
save_both(bar_plot(counts, "The Initial Distribution of the FMS", "FMS", 0.9), "initialclass_distribution", 8, 6)

# same grouping as the Poisson GLM in run_real_data.R
classes <- findInterval(fms, c(1, 4, 7))
counts <- data.frame(Class = factor(0:3), Count = tabulate(classes + 1, nbins = 4))
save_both(bar_plot(counts, "Distribution of the Classes", "Class", 0.35), "4class_distribution", 8, 6)


# Figure D2 -------------------------------------------------------------------

heat_colors <- colorRampPalette(c("#D73027", "#FC8D59", "white", "#91BFDB", "#4575B4"))(200)
features <- c("HR", "PC_HR", "HR_MIN", "HR_MAX", "HRV", "PC_HRV", "HRV_MIN",
              "HRV_MAX", "PC_BR", "GaSR", "PC_GaSR", "GaSR_MIN", "GaSR_MAX")

draw_13 <- function(r) {
  corrplot(r, method = "color", col = heat_colors, addCoef.col = "black", tl.col = "black",
           tl.cex = 0.8, cl.cex = 0.8, cl.pos = "r", number.cex = 0.8, addgrid.col = "grey")
}
draw_39 <- function(r) {
  corrplot(r, method = "color", col = heat_colors, tl.col = "black", tl.cex = 1.8, cl.cex = 1.8, cl.pos = "r")
}

r13 <- cor(d1[, 2:14])
dimnames(r13) <- list(features, features)
png("figures/heatmap13features.png", width = 8, height = 8, units = "in", res = 600)
draw_13(r13)
dev.off()
setEPS()
postscript("figures/heatmap13features.eps", width = 8, height = 8, horizontal = FALSE)
draw_13(r13)
dev.off()

# the last 39 feature columns are the lags t-3, t-2 and t-1
cols <- (ncol(d10) - 39):(ncol(d10) - 1)
lag_names <- paste0(rep(features, 3), "(t_", rep(3:1, each = 13), ")")
r39 <- cor(d10[, cols])
dimnames(r39) <- list(lag_names, lag_names)
png("figures/heatmap39of130features.png", width = 18, height = 18, units = "in", res = 600)
draw_39(r39)
dev.off()
setEPS()
postscript("figures/heatmap39of130features.eps", width = 18, height = 18)
draw_39(r39)
dev.off()
