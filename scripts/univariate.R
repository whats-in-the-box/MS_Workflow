library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(extrafont)
library(RColorBrewer)
# also need to import VennDiagram

##### Univariate Statistics and Fold Change

# declare functions
# put this in helpers/utils.R, internal
#
#' Get fold changes from a data matrix
#'
#' \code{get_fc} calculates fold changes from a data matrix with the
#'     column \code{Label} specifying membership of the samples.
#'
#' @param d3_mod A data frame or tibble. Normalized data matrix with columns of
#'     features and rows of samples.
#' @param linear Logical if the resulting fold changes should be linear or log2.
#'     Default to linear.
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @examples
#' iris %>%
#'     rename(Label = Species) %>%
#'     get_fc()
get_fc <- function(d3_mod, linear = TRUE) {
  if (isTRUE(linear)) {
    fc_col_name <- "FC(lin)"
    temp_d <- d3_mod %>%
      mutate(across(-Label, gtools::logratio2foldchange)) %>%
      group_by(Label) %>%
      summarise(across(everything(), mean))
    (temp_d[2, -1] / temp_d[1, -1]) %>%
      pivot_longer(everything(), names_to = "variable", values_to = fc_col_name)
  } else {
    fc_col_name <- "FC(log2)"
    temp_d <- d3_mod %>%
      group_by(Label) %>%
      summarise(across(everything(), mean))
    (temp_d[2, -1] - temp_d[1, -1]) %>%
      pivot_longer(everything(), names_to = "variable", values_to = fc_col_name)
  }
}

# -----------------------------------------------------------------------------
#' Compute univariate stats for a data matrix
#'
#' \code{do_univariate} computes univariate stats between two groups from a
#'     data matrix with the column \code{Label} specifying membership of the
#'     samples.
#'
#' @param d3_mod A data frame or tibble. Normalized data matrix with columns of
#'     features and rows of samples.
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @export
#'
#' @examples
#' iris %>%
#'     rename(Label = Species) %>%
#'     filter(Label %in% sample(unique(Label), 2)) %>%
#'     do_univariate()
do_univariate <- function(d3_mod) {
  d3_mod %>%
    pivot_longer(cols = -Label, names_to = "variable", values_to = "value") %>%
    group_nest(variable) %>%
    mutate(
      pT = purrr::map(data, ~ t.test(value ~ Label, data = .x)$p.value),
      pW = purrr::map(data, ~ wilcox.test(value ~ Label, data = .x)$p.value)
    ) %>%
    tidyr::unnest(cols = c("pT", "pW")) %>%
    mutate(
      BHT = p.adjust(pT, method = "BH"),
      BHW = p.adjust(pW, method = "BH")
    ) %>%
    select(-data) %>%
    {
      purrr::reduce(list(., get_fc(d3_mod), get_fc(d3_mod, linear = FALSE)), left_join)
    } %>%
    select("variable", "pT", "BHT", "pW", "BHW", "FC(lin)", "FC(log2)")
}

# -----------------------------------------------------------------------------
# from visual_functions.R in generate_html_reports repo
#
#' Make a volcano plot
#'
#' Another function to make a volcano plot. It's customized to the way CS
#' likes it.
#'
#' @param res A data frame or tibble. DE results table.
#'
#' @param df_filt A data frame or tibble. A smaller DE results table with only
#'     DE features.
#' @param fdr Numeric. The FDR cutoff to color DE features.
#' @param log2fc Numeric. The log2 fold change cutoff to color DE features.
#' @param feature_col Character. The column name for features.Used to retrieve
#'     labels for top 10 DE features.
#' @param padj_col Character, Optional. The column name for adjusted P values.
#'     Used for plotting. Default to "padj".
#' @param log2fc_col Character, Optional. The column name for log2 fole changes.
#'     Used for plotting. Default to "log2FoldChange".
#'
#' @import dplyr
#' @import ggplot2
#' @export
#'
#' @examples
#' todo
deseq2_volcano <- function(res, df_filt, fdr = fdr, log2fc = log2fc, feature_col,
                           padj_col = NULL, log2fc_col = NULL) {
  if (!is.null(padj_col)) padj_col <- padj_col else padj_col <- "padj"
  if (!is.null(log2fc_col)) log2fc_col <- log2fc_col else log2fc_col <- "log2FoldChange"
  if (!"-log10padj" %in% colnames(res)) res <- res %>% mutate(`-log10padj` = -log10(padj_col))
  # data wrangling
  # # todo
  # # need to change the below line to a df with feature name as a column to
  # # better generalize across different data
  # dge_df <- as.data.frame(res) %>% rownames_to_column(var = "gene")
  dge_df <- res
  df_filt_up <- df_filt %>% filter(status == "Up")
  df_filt_down <- df_filt %>% filter(status == "Down")
  up_no <- nrow(filter(df_filt, status == "Up"))
  down_no <- nrow(filter(df_filt, status == "Down"))
  # truncate feature names if too long
  df_filt <- df_filt %>%
    mutate(across(variable, ~ stringr::str_trunc(., width = 15, ellipsis = "")))
  up_top10 <- df_filt_up %>%
    arrange(padj_col) %>%
    dplyr::slice(1:10)
  down_top10 <- df_filt_down %>%
    arrange(padj_col) %>%
    dplyr::slice(1:10)
  ###
  volcano_xlim <- max(na.omit(dge_df[log2fc_col]))
  volcano_ylim <- max(-log10(na.omit(dge_df[padj_col])))
  x_anno <- 0.7 * -volcano_xlim
  y_anno <- 0.9 * volcano_ylim
  ###
  feature_col <- sym(feature_col)
  log2fc_col <- sym(log2fc_col)
  # draw plot
  ggplot(dge_df, aes(!!log2fc_col, `-log10padj`)) +
    geom_point(alpha = 0.4, size = 1.5, colour = "grey50", na.rm = TRUE) +
    scale_x_continuous(limits = c(-volcano_xlim, volcano_xlim)) +
    scale_y_continuous(limits = c(0, volcano_ylim)) +
    geom_point(
      data = df_filt_down, shape = 21, alpha = 0.6,
      size = 1.5, fill = "blue", colour = "blue", na.rm = TRUE
    ) +
    geom_point(
      data = df_filt_up, shape = 21, alpha = 0.6,
      size = 1.5, fill = "red", colour = "red", na.rm = TRUE
    ) +
    ###
    geom_point(
      data = up_top10, shape = 21, fill = "red",
      colour = "black", size = 2, na.rm = TRUE
    ) +
    ggrepel::geom_text_repel(
      data = up_top10, aes(label = !!feature_col), na.rm = TRUE,
      size = 4, family = "Lucida Sans", max.overlaps = 20
    ) +
    geom_point(
      data = down_top10, shape = 21, fill = "blue",
      colour = "black", size = 2, na.rm = TRUE
    ) +
    ggrepel::geom_text_repel(
      data = down_top10, aes(label = !!feature_col), na.rm = TRUE,
      size = 4, family = "Lucida Sans", max.overlaps = 20
    ) +
    ###
    theme_bw(base_size = 14) +
    labs(
      x = "log2-fold change",
      y = "-log 10 (padj)",
      title = sprintf("FDR: %.2f, log2FC: %.1f", fdr, log2fc)
    ) +
    theme(
      text = element_text(family = "Lucida Sans", face = "plain"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.title = element_text(size = 12.0),
      plot.title = element_text(size = 15.0, hjust = 0.5),
      aspect.ratio = 1
    ) +
    annotate("text",
      x = x_anno, y = y_anno, size = 5,
      label = sprintf("Up %i\nDown %i", up_no, down_no)
    )
}

# -----------------------------------------------------------------------------
# from visual_functions.R in generate_html_reports repo
#
#' Make a heatmap
#'
#' A convenient wrapper around \code{pheatmap} to make a Z-score pretty heatmap
#'     that takes matrix size into consideration.
#'
#' @param transformed_count A data frame or tibble. Normalized data matrix
#'     with columns of features and rows of samples.
#' @param df_filt A data frame or tibble. A smaller DE results table with only
#'     DE features.
#' @param feature_col Character. The column name for features.Used to retrieve
#'     labels for top 10 DE features.
#' @param anno A data frame. Named data frame of one column of grouping info.
#' @param top_n Numeric, Optional. The number of top features to plot
#'     (sorted by adjusted P values).Default to plot all DE features.
#' @param col_order A vector of column names, Optional. If provided, columns are
#'     not clustered and in the order provided. Default to clustered columns.
#' @param save Logical, Optional. Default to save.
#' @param padj_col Character, Optional. The column name for adjusted P values.
#'     Used for plotting. Default to "padj".
#'
#' @import dplyr
#' @export
#'
#' @examples
#' todo
deseq2_hm <- function(transformed_count, df_filt, feature_col, anno,
                      top_n = NULL, col_order = NULL,
                      save = TRUE, padj_col = NULL) {

  if (!is.null(padj_col)) padj_col <- padj_col else padj_col <- "padj"
  cluster_cols <- TRUE

  ##### colors for hm and anno
  inc_col <- colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"))(200)
  anno_col <- list()
  temp_a <- sample(RColorBrewer::brewer.pal(8, "Dark2"), n_distinct(anno[1]))
  names(temp_a) <- levels(anno[[1]])
  anno_col[[names(anno)]] <- temp_a

  ##### Plot sub-hm if top_n is provided
  if (!is.null(top_n)) {
    top_n_df_filt <- df_filt %>%
      arrange(padj_col) %>%
      dplyr::slice(1:top_n)
    sub_mat <- top_n_df_filt %>%
      select(feature_col) %>%
      left_join(transformed_count) %>%
      column_to_rownames("variable")
    # transformed_count[top_n_df_filt[[feature_col]], ]
  } else {
    sub_mat <- df_filt %>%
      select(feature_col) %>%
      left_join(transformed_count) %>%
      column_to_rownames("variable")
    # transformed_count[df_filt[[feature_col]], ]
  }

  ##### Specify parameters depending on hm size
  if (nrow(sub_mat) > 200) rowname_switch <- FALSE else rowname_switch <- TRUE
  if (ncol(sub_mat) > 60) {
    cellwidth <- 10
    fontsize <- 6
  } else {
    cellwidth <- 30
    fontsize <- 8
  }

  # if column order is specified
  if (!is.null(col_order)) {
    col_order <- col_order
    cluster_cols <- FALSE

    p3 <- pheatmap::pheatmap(
      sub_mat,
      annotation_col = anno,
      annotation_colors = anno_col,
      color = inc_col,
      cluster_cols = cluster_cols,
      labels_col = col_order,
      cellwidth = cellwidth,
      fontsize_row = fontsize,
      fontsize_col = fontsize,
      show_rownames = rowname_switch,
      treeheight_col = 25,
      treeheight_row = 25,
      angle_col = 45,
      fontsize_number = 6,
      scale = "row",
      annotation_names_col = FALSE,
      main = paste(
        "Clustering of normalized count (Z-score): ",
        dim(sub_mat)[1], "features"
      )
    )

  } else {

    p3 <- pheatmap::pheatmap(
      sub_mat,
      annotation_col = anno,
      annotation_colors = anno_col,
      color = inc_col,
      cellwidth = cellwidth,
      fontsize_row = fontsize,
      fontsize_col = fontsize,
      show_rownames = rowname_switch,
      treeheight_col = 25,
      treeheight_row = 25,
      angle_col = 45,
      fontsize_number = 6,
      scale = "row",
      annotation_names_col = FALSE,
      main = paste(
        "Clustering of normalized count (Z-score): ",
        dim(sub_mat)[1], "features"
      )
    )

  }
  if (isTRUE(save)) {
    ggsave("DE_features_heatmap_all_replicates.png",
      device = "png",
      plot = p3, dpi = 300, width = 17, height = 10
    )
    ggsave("DE_features_heatmap_all_replicates.pdf",
           device = "pdf",
           plot = p3, dpi = 300, width = 17, height = 10
    )
  }
}

# -----------------------------------------------------------------------------
#' Draw a pairwise venn diagram
#'
#' \code{make_venn} is a convenient wrapper around \code{VennDiagram::venn.diagram()} to draw a pairwise venn diagram.
#'
#' @param v_data A named list of vectors (e.g., integers, chars), with each component
#' corresponding to a separate circle in the Venn diagram.
#'
#' @return a \code{ggplot} of the venn diagram
#' @export
#'
#' @examples
#' make_venn(list("A" = 1:100, "B" = 91:150))
make_venn <- function(v_data) {
  my_col <- RColorBrewer::brewer.pal(8, "Pastel2")
  pairwise_venn_set12 <- VennDiagram::venn.diagram(
    x = v_data,
    filename = NULL,
    category.names = c(
      sprintf("%s\n %i", names(v_data[1]), length(v_data[[1]])),
      sprintf("%s\n %i", names(v_data[2]), length(v_data[[2]]))
    ),
    margin = 0.2,
    resolution = 300,
    euler.d = FALSE,
    scaled = FALSE,
    fill = my_col[1:2],
    lwd = 1,
    alpha = c(0.7, 0.7),
    col = c("grey60", "grey60"),
    cat.fontface = c("bold", "bold"),
    # cat.fontfamily = c("Lucida Sans", "Lucida Sans"),
    cat.dist = c(0.2, 0.2),
    ext.text = TRUE
    # main = sprintf("Overlapped features, FDR: %.2f", fdr),
    # main.fontface = "bold",
    # main.fontfamily = "Lucida Sans",
    # main.cex = 1.5,
    # main.pos = c(0.5, 0.8)
  )
  return(ggplotify::as.ggplot(grid::grobTree(pairwise_venn_set12)))
}
################

# load("normalization.rda")
#
# # first do data formatting
# d3_mod <- t(norm_d1) %>%
#   as_tibble() %>%
#   mutate(across(-Label, as.numeric)) %>%
#   rename_with(str_trim)
#
# # then do univariate analysis
# uni_res <- do_univariate(d3_mod)
#
# readr::write_csv(uni_res[,1:7], file = "univariate_results.csv")
#
#
# #-------------------------------- visualizations
# fdr <- 0.05
# log2fc <- 0
#
# # format univariate results tibble for plotting
# uni_res <- uni_res %>%
#   rowwise() %>%
#   mutate(padj = min(c(BHT, BHW))) %>%
#   # get the lowest padj
#   ungroup() %>%
#   mutate(`-log10padj` = -log(padj))
#
# # get DE features only tibble
# uni_res_filt <- uni_res %>%
#   # filter(BHT < fdr | BHW < fdr) %>%
#   filter(BHT < fdr & BHW < fdr) %>%
#   mutate(status = if_else("FC(log2)" < 0, "Down", "Up"))
#
#
# # 1) volcano plot
# deseq2_volcano(uni_res, uni_res_filt,
#   fdr = fdr, log2fc = log2fc,
#   "variable", padj_col = "padj", log2fc_col = "FC(log2)"
# )
#
# # -----------------------------------------------------------------------------
# # prep the annotation for hm
# anno <- data.frame(Label = as.factor(t(norm_d1)[, "Label"]))
#
# # prep the transformed matrix for hm
# d1_mod <- norm_d1 %>%
#   rownames_to_column("variable") %>%
#   mutate(across(-variable, as.numeric))
#
# # sanity check the variable names matched
# all(uni_res_filt$variable %in% d1_mod$variable)
#
# # 2) heatmap/multiple heatmaps?
# deseq2_hm(d1_mod, uni_res_filt, "variable", anno,
#   top_n = NULL, col_order = NULL,
#   save = FALSE, padj_col = NULL
# )
#
# # deseq2_hm(d1_mod, uni_res_filt, "variable", anno,
# #           top_n = NULL, col_order = NULL,
# #           save = TRUE, padj_col = NULL
# # )
#
# # -----------------------------------------------------------------------------
# # 3) venn diagram between t-test and wilcox test
# plot_set1 <- uni_res %>%
#   filter(BHT < fdr) %>%
#   pull(variable)
# plot_set2 <- uni_res %>%
#   filter(BHW < fdr) %>%
#   pull(variable)
# # v_data <- list("T Test" = plot_set1, "Wilcoxon Test" = plot_set2)
# v_data <- list("Wilcoxon Test" = plot_set2, "T Test" = plot_set1)
#
# # loadfonts(device="postscript")
# v1 <- make_venn(v_data)
# ggsave("univariate_venn.png", v1, device = "png", width = 6, height = 6)
# ggsave("univariate_venn.pdf", v1, device = "pdf", width = 6, height = 6)
#
#
# # pdf(file = "univariate_venn.pdf", family = "Lucida Sans Unicode", width = 6, height = 6)
# # plot(v1)
# # dev.off()
