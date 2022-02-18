
#'------------------------------------------------------------------------------
#' Title
#'
#' @param filename
#' @param featurecol
#' @param samplecol
#'
#' @return
#' @export
#'
#' @examples
upload.file <- function(filename, featurecol) {
  fh <- read.csv(filename, header = TRUE)
  samplecol = 2:length(fh)
  fileparsed <- fh[,c(featurecol, samplecol)]
  fileparsed <- as.data.frame(t(fileparsed));
  colnames(fileparsed) <- fileparsed[1,]; fileparsed <- fileparsed[-1,]
  return(fileparsed)
}


#'------------------------------------------------------------------------------
#' Title
#'
#' @param input_file
#'
#' @return
#' @export
#'
#' @examples
presence.absence <- function(input_file){
  library(tidyr)
  library(dplyr)
  t1 <- input_file %>%
    na_if(0)
  pheno_groups <- unique(t1$Label); group_num <- length(pheno_groups)

  for (i in 1:group_num){
    t2 <- t1 %>%
      filter(Label == pheno_groups[i])

    t1[paste0('missing', pheno_groups[i]),] <- as.matrix(apply(t2, 2, function(x) (sum(is.na(x)))/dim(t2)[1]))
    t1 <- rbind.data.frame(slice(t1[paste0('missing', pheno_groups[i]),]),t1[-dim(t1)[1],])
  }
  t1 <- as.data.frame(t(t1[,-1]))
  keep_file <- data.frame()
  reject_file <- data.frame()
  #  write(reject_file, "absentFeatures.txt", append = TRUE)
  for (i in 1:dim(t1)[1]){
    if( mean(as.numeric(unlist(t1[i,1:group_num]))) < .50 ){
      keep_file <- rbind.data.frame(t1[i,], keep_file)
    }
    else {
      reject_file <- rbind.data.frame(t1[i,], reject_file)
    }
  }
  return(keep_file)
}


#'------------------------------------------------------------------------------
#' Track missing signal at different stages
#'
#' @param stage_name Character. The stage of workflow.
#' @param no_of_missing Numeric. The number of missing data.
#'
#' @return A data frame.
#'
track_missing <- function(stage_name, no_of_missing) {
  missing_df <- missing_df %>%
    dplyr::add_row(
      Stage = stage_name,
      No_of_missing = no_of_missing
    )
  return(missing_df)
}


#'------------------------------------------------------------------------------
#' do normalization on log transformed data matrix
#'
#' @param df
#' @param labels_d1
#' @param method
#'
#' @return A data frame
#'
do_normalization_short <- function(df, labels_d1, method = "EigenMS", ...){

  ### data wrangling
  grps = as.factor(t(labels_d1))

  df = df[-1,] %>%
    rownames_to_column(.,var = "rowname") %>%
    mutate(across(-rowname, as.numeric)) %>%
    column_to_rownames(., var = "rowname")

  choices <- function(method) {
    switch(method,
      "Cubic Spline" = cubic_norm(df, ...),
      "EigenMS" = EigenMS_norm(df, grps),
      "Invariant" = invariant_norm(df, ...),
      stop("Normalization method not included.")
    )
  }

  norm_df <- choices(method)

  ### add label back to normalized data frame
  norm_df <- rbind.data.frame("Label" = as.character(labels_d1),
                              norm_df)
  return(norm_df)
}


#'------------------------------------------------------------------------------
#' Title
#'
#' @param df
#' @param label
#' @param WORKING_DIR
#'
#' @return
#' @export
#'
#' @examples
pls_da <- function(df, label, WORKING_DIR) {
  my.plsda <- ropls::opls(t(df), label, fig.pdfC = "none",
                          info.txtC = file.path(WORKING_DIR, "PLS-DA_info.txt") )
  # PLS-DA_overview
  plsda_plot <- ropls::plot(my.plsda, parAsColFcVn = label, fig.pdfC = "interactive",parLabVc = rep('o', length(label)))
  my.vip <- sort(ropls::getVipVn(my.plsda), decreasing = T)
  return(list(my.plsda, plsda_plot, my.vip))
}


#'------------------------------------------------------------------------------
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


#'------------------------------------------------------------------------------
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


#'------------------------------------------------------------------------------
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

#'------------------------------------------------------------------------------
#' Draw a multi-way Venn diagram
#'
#' \code{make_multi_venn} is a convenient wrapper around \code{VennDiagram::venn.diagram()} to draw a multi-way Venn diagram.
#' @param v_data
#'
#' @return
#' @export
#'
#' @examples
make_multi_venn <- function(v_data) {
  my_col <- RColorBrewer::brewer.pal(8, "Pastel2")
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  multi_venn_set <- VennDiagram::venn.diagram(
    x = v_data,
    filename = NULL,
    for (a in 1:length(v_data)) {
       category.names = c(sprintf("%s\n %i", names(v_data[a]), length(v_data[[a]])))
    },
    margin = 0.2,
    resolution = 300,
    euler.d = FALSE,
    scaled = FALSE,
    fill = my_col[1:3],
    lwd = 1,
    alpha = rep(0.7, times = length(v_data)),
    col = rep("grey60", times = length(v_data)),
    cat.fontface = rep("bold", times = length(v_data)),
    # cat.fontfamily = c("Lucida Sans", "Lucida Sans"),
    cat.dist = rep(0.2, times = length(v_data)),
    ext.text = TRUE
    # main = sprintf("Overlapped features, FDR: %.2f", fdr),
    # main.fontface = "bold",
    # main.fontfamily = "Lucida Sans",
    # main.cex = 1.5,
    # main.pos = c(0.5, 0.8)
  )
  return(ggplotify::as.ggplot(grid::grobTree(multi_venn_set)))
}

#-------------------------------------------------------------------------------
# Define save plot function ----
# save a ggplot graph in both PDF and jpg format. The filename should be
# PATH+filename_without_extension(.jpg, .pdf). 'fig' is the
# variable returned by ggplot.
ggsave_both <- function(filename_no_ext, fig_width = 15, fig_height = 10, fig) {
        ggsave(paste0(filename_no_ext, '.pdf'), device = 'pdf',
               width = fig_width, height = fig_height, fig)

        ggsave(paste0(filename_no_ext, '.jpg'), device = 'jpg',
               width = fig_width, height = fig_height, fig)
}
