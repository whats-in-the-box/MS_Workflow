
library(ProteoMM)
library(tidyverse)
library(ggplot2)

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
#' @param samples for about 200 data points use samples = 0.33, for about 2000 data points use samples = 0.05, for about 10000 data points use samples = 0.02
#'
#' @return a cubic-spline normalized data matrix
cubic_norm <- function(df, ...){
  df = as.matrix(df)
  if (min(df)<0){
    print(paste("Data contains negative values, min(df) = ", min(df),sep = "") )
    norm_df <- affy::normalize.qspline(df + abs(min(df)), ...)
    norm_df = as.data.frame(norm_df)
    rownames(norm_df) = rownames(df)
    colnames(norm_df) = colnames(df)

    return(norm_df)
    print(paste("Cubic Spline normalization is completed after shifting dataframe by ", abs(min(df)), sep = "") )
    print("Attention: Fold change is not accurate after shifting dataframe.")
  } else {
    norm_df <- affy::normalize.qspline(df, ...)
    norm_df = as.data.frame(norm_df)
    rownames(norm_df) = rownames(df)
    colnames(norm_df) = colnames(df)

    return(norm_df)
  }
}

# ------------------------------------------------------------------------------
#' \code{EigenMS_norm} normalize a data matrix, preferably a log2 transformed data matrix,
#'   with eigenMS method from ProteoMM package.
#'
#' @param df
#' @param grps
#' @return a EigenMS normalized data matrix
EigenMS_norm <- function(df,grps){

  df <- as.data.frame(df) #feature x sample
  df <- cbind.data.frame(row.names(df), df)

  intsCols <- 2:dim(df)[2]
  m_logInts = make_intencities(df, intsCols)
  #m_logInts = convert_log2(type.convert(m_logInts))
  m_mets.info = as.matrix(make_meta(df, 1))#[-1])

  # set.seed(123)
  # first portion of EigenMS: identify bias trends
  hs_m_ints_eig1 = eig_norm1(m=m_logInts,treatment = grps,
                             prot.info=m_mets.info)
  # second portion of EigenMS: EigenMS normalization
  hs_m_ints_norm_1bt = eig_norm2(rv=hs_m_ints_eig1)
  norm_df <- as.data.frame(hs_m_ints_norm_1bt[["norm_m"]])
  return(norm_df)
}

# ------------------------------------------------------------------------------
#' \code{invariant_norm} is a convenient wrapper around \code{lumi::rankinvariant()},
#'   to generate a normalized data matrix from a preferably log2 transformed data matrix,
#'   \code{rownames()} and \code{colnames()} of the output data matrix needs to
#'   be specified as the same as the input matrix.
#'
#' @param df
#' @return an invariant normalized data matrix
invariant_norm <- function(df, ...){
  df= as.matrix(df)

  norm_df <- lumi::rankinvariant(as.matrix(df), ...)
  norm_df <- as.data.frame(norm_df)
#  rownames(norm_df) = rownames(df)
#  colnames(norm_df) = colnames(df)
  return(norm_df)
}

