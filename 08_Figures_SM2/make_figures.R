# Figures SM.2.1-SM.2.9 of the Supplementary Material from the output of
# run_simulation.R in results/. The figures go to figures/ under the file
# names used in the paper, and the median, mean and sd of each box in
# SM.2.5-SM.2.9 to summaries/.
#
# Run from this folder:  Rscript make_figures.R

library(ggplot2)

dir.create("figures", showWarnings = FALSE)
dir.create("summaries", showWarnings = FALSE)

save_both <- function(plot, name, width, height) {
  ggsave(file.path("figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 600)
  setEPS()
  postscript(file.path("figures", paste0(name, ".eps")), width = width, height = height, family = "Times")
  print(plot)
  dev.off()
}

times_theme <- theme_classic(base_size = 14, base_family = "Times") +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "plain", family = "Times"),
    axis.title = element_text(size = 18, face = "plain", family = "Times"),
    axis.text = element_text(size = 16, face = "plain", family = "Times"),
    strip.text = element_text(size = 16, face = "plain", family = "Times")
  )

ratio_summary <- function(x) {
  c(median = median(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE), std_dev = sd(x, na.rm = TRUE))
}


# SM.2.1 and SM.2.2: sample eigenvalues ---------------------------------------

toeplitz_cov <- function(rho, p) rho^abs(outer(1:p, 1:p, "-"))

# half of the eigenvalues equal 1 and the other half 0.5
half_cov <- function(p) diag(rep(c(1, 0.5), c(ceiling(p / 2), p - ceiling(p / 2))))

eigen_plot <- function(file, true_eigen, title, half_width) {
  d <- read.csv(file)
  cols <- paste0("E", seq_along(true_eigen))
  long <- data.frame(eigen_index = factor(rep(cols, each = nrow(d)), levels = cols), e_val = unlist(d[cols]))
  truth <- data.frame(eigen_index = factor(cols, levels = cols), true_val = sort(true_eigen, decreasing = TRUE))
  ggplot(long, aes(x = eigen_index, y = e_val)) +
    geom_boxplot(fill = "lightblue", color = "black", width = 0.75) +
    geom_segment(
      data = truth,
      aes(x = as.numeric(eigen_index) - half_width, xend = as.numeric(eigen_index) + half_width,
          y = true_val, yend = true_val),
      color = "red", linewidth = 0.8
    ) +
    labs(title = title, x = "Eigenvalue Index", y = "Eigenvalue") +
    times_theme
}

for (p in c(5, 10)) {
  for (rho in c(-0.5, 0, 0.5)) {
    plot <- eigen_plot(
      sprintf("results/SM2_1_2/toeplitz_p%d_rho%g.csv", p, rho),
      eigen(toeplitz_cov(rho, p))$values,
      bquote("Empirical Eigenvalues for " ~ p == .(p) ~ " and " ~ rho == .(rho)),
      0.375
    )
    save_both(plot, sprintf("plot_1_p=%d_rho=%g", p, rho), 8, 6)
  }
  plot <- eigen_plot(
    sprintf("results/SM2_1_2/independent_p%d.csv", p),
    diag(half_cov(p)),
    bquote("Empirical Eigenvalues for " ~ p == .(p) * ", Independent"),
    0.35
  )
  save_both(plot, sprintf("plot_1_p=%d_Independent", p), 8, 6)
  plot <- eigen_plot(
    sprintf("results/SM2_1_2/dependent_p%d.csv", p),
    diag(half_cov(p)),
    bquote("Empirical Eigenvalues for " ~ p == .(p) * ", Dependent"),
    0.35
  )
  save_both(plot, sprintf("plot_2_p=%d_Dependent", p), 8, 6)
}


# SM.2.3 and SM.2.4: n = 250 ---------------------------------------------------

n <- 250
p_n <- c(0.05, 0.1, 0.25, 0.5, 0.75)
p_vals <- ceiling(n * p_n)

# figure names: theta "alt" is beta = (1, -1, ...) and "all" is beta = (1, ..., 1)
for (theta in c("alt", "all")) {
  beta_design <- if (theta == "alt") "alt" else "ones"
  for (lambda_max in c(1, 10)) {
    d <- do.call(rbind, lapply(p_vals, function(p) {
      file <- sprintf("results/SM2_3_4/beta_%s_lambdamax%g_p%d.csv", beta_design, lambda_max, p)
      cbind(read.csv(file)[c("ratio_beta_OLS", "ratio_lambda")], p = p)
    }))
    d$ratio <- factor(d$p, levels = p_vals, labels = paste0(100 * p_n, "%"))

    beta_plot <- ggplot(d, aes(x = ratio, y = ratio_beta_OLS)) +
      geom_boxplot(fill = "lightblue", color = "black") +
      labs(
        title = expression(paste("OLS ", L[1], " Error Ratio for Estimated Parameters vs. ", p/n, ", n = ", 250)),
        x = expression(paste(p/n)),
        y = expression(L[1] * " Error Ratio")
      ) +
      times_theme
    lambda_plot <- ggplot(d, aes(x = ratio, y = ratio_lambda)) +
      geom_boxplot(fill = "lightcoral", color = "black") +
      labs(
        title = expression(paste(L[1], " Error Ratio for Empirical Eigenvalues vs. ", p/n, ", n = ", 250)),
        x = expression(paste(p/n)),
        y = expression(L[1] * " Error Ratio")
      ) +
      times_theme
    save_both(beta_plot, sprintf("Boxplot_Ratio_Beta_theta_%s_lambdaMax_%g", theta, lambda_max), 8, 6)
    save_both(lambda_plot, sprintf("Boxplot_Ratio_Lambda_theta_%s_lambdaMax_%g", theta, lambda_max), 8, 6)
  }
}


# SM.2.5 to SM.2.9: n = 1,000 to 10,000 ----------------------------------------

n_vals <- c(1000, 2500, 5000, 10000)
p_n <- c(0.25, 0.5, 0.75, 0.9, 0.95)
n_labels <- c("n == \"1,000\"", "n == \"2,500\"", "n == \"5,000\"", "n == \"10,000\"")

# one row per replication with columns n_val (facet label) and p_ratio
read_large <- function(file_of, cols) {
  out <- NULL
  for (n_val in n_vals) for (ratio in p_n) {
    d <- read.csv(file_of(n_val, ceiling(n_val * ratio)))[cols]
    d$n <- n_val
    d$p_ratio <- paste0(round(ratio * 100), "%")
    out <- rbind(out, d)
  }
  out$n_val <- factor(out$n, levels = n_vals, labels = n_labels)
  out
}

faceted <- function(d, y, fill, title, ylab) {
  ggplot(d, aes(x = factor(p_ratio), y = .data[[y]])) +
    geom_boxplot(fill = fill, color = "black") +
    labs(title = title, x = expression(paste(p/n)), y = ylab) +
    facet_wrap(~ n_val, nrow = 2, ncol = 2, labeller = label_parsed) +
    times_theme
}

# median, mean and sd of each box, as quoted in the text
box_summary <- function(d, y, figure, estimator) {
  s <- aggregate(d[[y]], list(n = d$n, p_ratio = d$p_ratio), ratio_summary)
  data.frame(figure = figure, estimator = estimator, s[1:2], s$x)
}

methods <- c("OLS", "St", "DSh", "GSR")
d <- read_large(function(n_val, p) sprintf("results/SM2_5_8/n%d_p%d.csv", n_val, p), paste0("ratio_beta_", methods))
summaries <- NULL
for (m in methods) {
  y <- paste0("ratio_beta_", m)
  title <- bquote(.(paste0(m, " ")) * L[1] * " Error Ratio for Estimated Parameters (Indep/Dep) vs. " * p/n)
  plot <- faceted(d, y, "lightblue", title, expression("L"[1] * " Error Ratio"))
  save_both(plot, sprintf("Faceted_%s_Boxplot_Ratio_Beta_n1000_n10000", m), 12, 6)
  summaries <- rbind(summaries, box_summary(d, y, "SM.2.5-SM.2.8", m))
}

d <- read_large(function(n_val, p) sprintf("results/SM2_9/n%d_p%d.csv", n_val, p),
                c("ratio_St_OLS", "ratio_DSh_OLS", "ratio_GSR_OLS"))
fills <- c(St = "lightblue", DSh = "lightgreen", GSR = "lightcoral")
for (m in names(fills)) {
  y <- paste0("ratio_", m, "_OLS")
  title <- bquote(L[1] * .(paste0(" Error Ratio for Estimated Parameters (", m, "/OLS) vs. ")) * p/n)
  plot <- faceted(d, y, fills[[m]], title, expression(L[1] * " Error Ratio"))
  save_both(plot, sprintf("Faceted_Boxplot_Ratio_%s_vs_OLS", m), 12, 6)
  summaries <- rbind(summaries, box_summary(d, y, "SM.2.9", m))
}

write.csv(summaries, "summaries/SM2_5_9_ratio_summaries.csv", row.names = FALSE)
