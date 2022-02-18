library(tidyverse)
library(ropls)

##### Multivariate Statistics

# declare functions

# ------------------------------------------------------------------------------
#' \code{do_multivariate} takes in a data matrix, feature by sample with a row
#' specifying group information, shows an overview of the data matrix, and do
#' diagnostic PCA, PLS-DA and OPLS-DA. The function returns overview plots and model details
#' of PCA, PLS-DA and OPLS-DA, and VIP scores of variables from PLS-DA.
#'
#' @param df data matrix, feature by sample with row \code{Label} specifying membership of the samples.
#' @param remove_outlier logical if or not remove inputted outlier names.
#' @param outliers a character vector of outlier names.
#'
#' @return my.vip dataframe containing VIP ranking of variables after PLS-DA
#'
#' @export "input_overview.pdf" overview of input data matrix
#' @export "outlier diagnostics.pdf" overview of diagnostic PCA
#' @export "pca_info.txt" model info of PCA
#' @export "PLS-DA_overview.pdf" overview of PLS-DA
#' @export "PLS-DA_info.txt" model info of PLS-DA
#' @export "PLS-DA scores.pdf" score plot of PLS-DA
#' @export "OPLS-DA_overview.pdf" overview of OPLS-DA
#' @export "OPLS-DA_info.txt" model info of OPLS-DA
#' @export "OPLS-DA scores.pdf" score plot of OPLS-DA
do_multivariate = function(df, remove_outlier = FALSE, outliers = NULL){

  # remove outliers based on Mahanabois distance
  if (isTRUE(remove_outlier)) {
    if (!is.null(outliers) & is.vector(outliers) == TRUE){

      # data wrangling
      outliers = outliers
      label_new = as.character(df[1, !(colnames(df) %in% outliers) ])
      df_new = df[-1,]%>%
        rownames_to_column(.,var = "rowname")%>%
        mutate(across(-rowname, as.numeric)) %>%
        column_to_rownames(., var = "rowname")%>%
        select(!(all_of(outliers)))

      # input overview
      pdf("input_overview.pdf",width = 15,height = 10)
      attach(df_new)
      ropls::view(df_new)
      dev.off()

      # new diagnostic PCA after removing outliers
      my.pca_new <- ropls::opls(as.data.frame(t(df_new)), parAsColFcVn = label_new, fig.pdfC = "none",info.txtC = "pca_info.txt")
      pdf("outlier_diagnostics.pdf",width = 15,height = 10)
      plot(my.pca_new, typeVc = "outlier", parAsColFcVn = label_new,fig.pdfC = "interactive")
      dev.off()

      # PLS-DA
      my.plsda_new <- ropls::opls(t(df_new),label_new,fig.pdfC = "none",info.txtC = "PLS-DA_info.txt")
      pdf("PLS-DA_overview.pdf", width = 15, height = 10)
      plot(my.plsda_new, parAsColFcVn = label_new, fig.pdfC = "interactive")
      dev.off()

      my.vip_new <- as.data.frame(sort(ropls::getVipVn(my.plsda_new), decreasing = T))
      colnames(my.vip_new)="VIP"
      return(my.vip_new)

      pdf("PLS-DA scores.pdf", width = 15, height = 10)
      plot(my.plsda_new, typeVc = "x-score", parAsColFcVn = label_new, parLabVc = rep('o', length(label_new)), fig.pdfC = "interactive")
      dev.off()

      # OPLS-DA
      my.oplsda_new <- ropls::opls(t(df_new), label_new, predI = 1, orthoI = NA, permI=10,fig.pdfC = "none",info.txtC = "OPLS-DA_info.txt")
      pdf("OPLS-DA_overview.pdf", width = 15, height = 10)
      plot(my.plsda_new, parAsColFcVn = label, fig.pdfC = "interactive")
      dev.off()

      pdf("OPLS-DA scores.pdf", width = 15, height = 10)
      plot(my.oplsda_new, typeVc = "x-score", parAsColFcVn = label_new, parLabVc = rep('o', length(label_new)), fig.pdfC = "interactive")
      dev.off()

    } else if (is.null(outliers)){
      print("Error: please supply outliers.")
    } else {
      print("Error: input is not valid, please provide a vector of outlier sample names.")
    }
  } else {
    # data wrangling
    label = as.character(df[1, ])
    df = df[-1,]%>%
      rownames_to_column(.,var = "rowname")%>%
      mutate(across(-rowname, as.numeric)) %>%
      column_to_rownames(., var = "rowname")

    # input overview
    pdf("input_overview.pdf",width = 15,height = 10)
    attach(df)
    ropls::view(df)
    dev.off()

    # Outlier diagnostics
    my.pca <- ropls::opls(as.data.frame(t(df)), parAsColFcVn = label, fig.pdfC = "none",info.txtC = "pca_info.txt")
    pdf("outlier_diagnostics.pdf",width = 15,height = 10)
    plot(my.pca, typeVc = "outlier", parAsColFcVn = label,fig.pdfC = "interactive")
    dev.off()

    # PLS-DA
    my.plsda <- ropls::opls(t(df),label,fig.pdfC = "none",info.txtC = "PLS-DA_info.txt")
    pdf("PLS-DA_overview.pdf", width = 15, height = 10)
    plot(my.plsda, parAsColFcVn = label, fig.pdfC = "interactive")
    dev.off()

    my.vip <- as.data.frame(sort(ropls::getVipVn(my.plsda), decreasing = T))
    colnames(my.vip)="VIP"
    return(my.vip)

    pdf("PLS-DA scores.pdf", width = 15, height = 10)
    plot(my.plsda, typeVc = "x-score", parAsColFcVn = label, parLabVc = rep('o', length(label)), fig.pdfC = "interactive")
    dev.off()

    # OPLS-DA
    my.oplsda <- ropls::opls(t(df), label, predI = 1, orthoI = NA, permI=10,fig.pdfC = "none",info.txtC = "OPLS-DA_info.txt")
    pdf("OPLS-DA_overview.pdf", width = 15, height = 10)
    plot(my.plsda, parAsColFcVn = label, fig.pdfC = "interactive")
    dev.off()

    pdf("OPLS-DA scores.pdf", width = 15, height = 10)
    plot(my.oplsda, typeVc = "x-score", parAsColFcVn = label, parLabVc = rep('o', length(label)), fig.pdfC = "interactive")
    dev.off()

  }

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
pls_da <- function(df, label, WORKING_DIR, ...) {
        my.plsda <- ropls::opls(t(df), label, fig.pdfC = "none",
                                info.txtC = file.path(WORKING_DIR, "PLS-DA_info.txt") )
        # PLS-DA_overview
        plsda_plot <- ropls::plot(my.plsda, parAsColFcVn = label, fig.pdfC = "interactive",parLabVc = rep('o', length(label)))
        my.vip <- sort(ropls::getVipVn(my.plsda), decreasing = T)
        return(list(my.plsda, plsda_plot, my.vip))
}
