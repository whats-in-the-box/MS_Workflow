
library(tidyverse)

##### Normality check

# declare functions

# ------------------------------------------------------------------------------
#' \code{normality_check} checks if a data matrix, feature x sample with a row \code{Label}
#'   specifying membership of the samples, is normally distributed using histogram and QQ plot.
#'
#' @param df
#' @return a graph exported in working directory
normality_check <- function(df){

  ### data wrangling
  df = as.data.frame(df)
  labels <- as.matrix(df[1,])
  grps = as.factor(t(labels))

  df = df[-1,]%>%
    rownames_to_column(.,var = "rowname")%>%
    mutate(across(-rowname, as.numeric)) %>%
    column_to_rownames(., var = "rowname")

  pdf("normality_check.pdf",width = 15, height = 10)
  par(mfrow = c(1,2))
  hist(unlist(t(df)),main = "Histogram",nbins=20,probability = T)
  lines(density(unlist(df)))
  car::qqPlot(unlist(t(df)),main = "QQ plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
  dev.off()
}
