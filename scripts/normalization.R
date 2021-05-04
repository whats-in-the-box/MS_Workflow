
library(ProteoMM)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

##### Normalization

# declare functions
# put this in helpers/utils.R, internal

# ------------------------------------------------------------------------------
#' \code{cubic_norm} is a convenient wrapper around \code{affy::normalize.qspline()},
#'   to generate a normalized data matrix from a preferably log2 transformed data matrix,
#'   \code{rownames()} and \code{colnames()} of the output data matrix needs to
#'   be specified as the same as the input matrix.
#'
#' @param df
#' @return
cubic_norm <- function(df){
  df= as.matrix(df)

  norm_df <- affy::normalize.qspline(df, p.min = 0, p.max = 10, fit.iters = 50,
                                     min.offset = 5, samples = 75)
  norm_df <- as.data.frame(norm_df)
  rownames(norm_df) = rownames(df)
  colnames(norm_df) = colnames(df)
  return(norm_df)
}

# ------------------------------------------------------------------------------
#' \code{EigenMS_norm} normalize a data matrix, preferably a log2 transformed data matrix,
#'   with eigenMS method from ProteoMM package.
#'
#' @param df
#' @param grps
#' @return
EigenMS_norm <- function(df,grps){

  df <- as.data.frame(df) #feature x sample
  df <- cbind.data.frame(row.names(df), df)

  intsCols <- 2:dim(df)[2]
  m_logInts = make_intencities(df, intsCols)
  #m_logInts = convert_log2(type.convert(m_logInts))
  m_mets.info = as.matrix(make_meta(df, 1))#[-1])

  set.seed(123)
  # first portion of EigenMS: identify bias trends
  hs_m_ints_eig1 = eig_norm1(m=m_logInts,treatment = grps,
                             prot.info=m_mets.info)
  # second portion of EigenMS: EigenMS normalization
  hs_m_ints_norm_1bt = eig_norm2(rv=hs_m_ints_eig1)
  norm_df <- as.data.frame(hs_m_ints_norm_1bt[["norm_m"]])
  return(norm_df)
}

# ------------------------------------------------------------------------------
#' \code{make_hist} is a visualization function that produce side-by-side histogram
#'   from log2 transformed data matrix and normalized data matrix using \code{hist()} and
#'   \code{MASS::truehist()} functions, respectively.
#'
#' @param df
#' @param norm_df
#' @import
#' @return
 make_hist = function(df,norm_df){
  df = as.data.frame(df)
  norm_df = as.data.frame(norm_df)

  par(mfrow = c(1,2))
  hist(unlist(t(df)),main = "Before normalization",nbins=20,probability = T)
  MASS::truehist(unlist(t(norm_df)), main = "After normalization", nbins = 50)
  par(mfrow=c(1,1))
}

# ------------------------------------------------------------------------------
#' \code{make_qqplot} is a visualization function that produce side_by_side qqplot
#'   from log2 transformed data matrix and normalized data matrix using \code{car::qqplot()}.
#'
#' @param df
#' @param norm_df
#' @return
make_qqplot = function(df,norm_df){
  df = as.data.frame(df)
  norm_df = as.data.frame(norm_df)

  par(mfrow = c(1,2))
  car::qqPlot(unlist(t(df)),main = "Before normalization", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
  car::qqPlot(unlist(t(norm_df)),main = "After normalization", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
  par(mfrow=c(1,1))
}

# ------------------------------------------------------------------------------
#' \code{make_pca} is a visualization function that wraps \code{pcaMethods::pca()}
#'   and \code{ggplot()} together. The function takes a normalized data matrix and
#'   a \code{labels} vector specifying membership of the samples.
#'
#' @param norm_df
#' @param labels
#' @return
 make_pca = function(norm_df,labels){
  # run PCA
  pc = pcaMethods::pca(t(norm_df), nPcs = 3,method = "svd", scale = "pareto",center = TRUE)
  norm_df_2 = rbind.data.frame(labels,norm_df)%>%
    t(.)%>%
    as.data.frame(.)%>%
    rownames_to_column(.)%>%
    as_tibble(.)%>%
    select(rowname,Label)%>%
    column_to_rownames(., var = "rowname")%>%
    merge(.,pcaMethods::scores(pc), by=0)
  # draw PCA
  p_pca = ggplot(norm_df_2, aes(PC1, PC2,colour=Label)) + geom_point() + stat_ellipse() +
    xlab(paste("PC1", round(pc@R2[1] * 100, digits = 1), "% of the variance")) +
    ylab(paste("PC2", round(pc@R2[2] * 100, digits = 1), "% of the variance")) #+ ggtitle(label = "After normalization")
  return(p_pca)
}

# ------------------------------------------------------------------------------
#' \code{normalization} is a wrapper function to achieve normalization of the data matrix,
#'   feature x sample the column \code{Label} specifying membership of the samples.
#'
#' @param df
#' @param method
#' @return
normalization <- function(df, method = "EigenMS"){

  ### data wrangling
  df = as.data.frame(df)
  labels <- as.matrix(df[1,])
  grps = as.factor(t(labels))

  method = method
  choices = switch(method, "Cubic Spline" = "Cubic Spline", "EigenMS" = "EigenMS")

  df = df[-1,]%>%
    rownames_to_column(.,var = "rowname")%>%
    mutate(across(-rowname, as.numeric)) %>%
    column_to_rownames(., var = "rowname")

  ### make normalized data matrix
  if (choices == "Cubic Spline"){
    norm_df = cubic_norm(df)
  } else if (choices == "EigenMS"){
    norm_df = EigenMS_norm(df,grps)
  }

  ### Histogram
  pdf('hist.pdf')
  make_hist(df,norm_df)
  dev.off()

  ### QQ plot
  pdf('qqplot.pdf')
  make_qqplot(df,norm_df)
  dev.off()

  ### PCA
  make_pca(norm_df,labels)
  ggsave("pca.pdf", width = 15, height = 10, p_pca)
  dev.off()

  ### add label back to normalized data frame
  norm_df = rbind.data.frame(labels, norm_df)
  return(norm_df)
}
