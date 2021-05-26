
do_normalization_short <- function(df, method = "EigenMS"){

  ### data wrangling
  label <- as.matrix(df[1,])
  grps = as.factor(t(label))

  df = df[-1,] %>%
    rownames_to_column(.,var = "rowname") %>%
    mutate(across(-rowname, as.numeric)) %>%
    column_to_rownames(., var = "rowname")

  choices <- function(method) {
    switch(method,
      "Cubic Spline" = cubic_norm(df),
      "EigenMS" = EigenMS_norm(df,grps),
      stop("Normalization method not included.")
    )
  }

  norm_df <- choices(method)

  ### add label back to normalized data frame
  norm_df <- rbind.data.frame(label, norm_df)

  return(norm_df)
}


# ------------------------ Visualization

#' ggplot version of \code{MASS::truehist()}
#'
#' @param data A numeric vector.
#' @param title Character. Plot title.
#'
#' @return A histogram.
#'
ggplot_truehist <- function(data, title) {
  data <- as.numeric(data)
  ggplot() +
    aes(data) +
    geom_histogram(aes(y = ..density..), bins = 50,
                   fill = "cornflowerblue", color = "gray30") +
    # geom_density(fill = NA, size = 1, color = "black") +
    labs(title = title) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      aspect.ratio = 1
    )
}


#' ggplot version of \code{car::qqPlot()}
#'
#' Code from https://www.tjmahr.com/quantile-quantile-plots-from-scratch/.
#'
#' @param data A numeric vector.
#'
#' @return A qqplot.
#'
ggplot_carqq <- function(data, title) {
  # the function se_z() applies the equation to a z-score (the theoretical
  # values on the x axis)
  se_z <- function(z, n) {
    sqrt(pnorm(z) * (1 - pnorm(z)) / n) / dnorm(z)
  }

  df <- tibble(
    x_sample = sort(as.numeric(data)),
    quantile = ppoints(length(x_sample)),
    z_theoretical = qnorm(quantile, 0, 1)
  )

  limit <- max(df$z_theoretical)

  conf_band <- tibble(
    z = seq(-(limit + 0.2), limit + 0.2, length.out = 300),
    n = length(df$x_sample),

    robust_sd = IQR(df$x_sample, na.rm = TRUE) / 1.349,
    robust_line = median(df$x_sample) + robust_sd * z,
    robust_se = robust_sd * se_z(z, n),
    robust_upper = robust_line + 2 * robust_se,
    robust_lower = robust_line - 2 * robust_se,
  )

  # ----- plot
 ggplot(df) +

    # plot data points
    geom_point(aes(x = z_theoretical, y = x_sample)) +

    # add reference line
    geom_abline(
      aes(intercept = mean, slope = sd),
      data = tibble(
        sd = IQR(df$x_sample, na.rm = TRUE) / 1.349,
        mean = median(df$x_sample)
      ),
      color = "blue",
      size = 1
    ) +

    # add confidence band
    geom_ribbon(
      aes(x = z, ymax = robust_upper, ymin = robust_lower),
      data = conf_band,
      fill = NA,
      color = "blue",
      linetype = 2,
      show.legend = FALSE,
      size = 1
    ) +

    labs(
      x = "theoretical quantiles",
      y = "sample quantiles",
      title = title
    ) +
    theme_bw() +
   scale_x_continuous(limits = c(-limit, limit)) +
   scale_y_continuous(limits = c(min(df$x_sample), max(df$x_sample)))
}


#' ggplot of \code{pcaMethods::pca()}
#'
#' @param data A data matrix.
#' @param labels A data matrix of group labels for samples.
#' @param title Character. Plot title
#'
#' @return A scatter plot with group ellipses.
#'
ggplot_pca <- function(data, labels, title) {
  data <- as.matrix(data)
  class(data) <- "numeric"
  pc1 <- pca(t(data), scale = "pareto")
  pc1merged <- merge(t(rbind.data.frame(t(labels), data)),
                     scores(pc1), by=0)
  ggplot(pc1merged, aes(PC1, PC2, colour=Label)) +
    geom_point() +
    stat_ellipse() +
    xlab(paste("PC1", round((pc1@R2[1] * 100), digits = 1), "% of the variance")) +
    ylab(paste("PC2", round((pc1@R2[2] * 100), digits = 1), "% of the variance")) +
    ggtitle(label = title)
}
