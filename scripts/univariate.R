library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(extrafont)
library(RColorBrewer)

##### Univariate Statistics and Fold Change

save(d3, d3_mod, wt, file = "for_univariate.Rda")

# declare function
#
# only works for pairwise foldchange
# put this in helpers/utils.R, internal
#' @importFrom magrittr %>%
#' @import dplyr
get_fc <- function(d3_mod, linear = TRUE) {
  if (isTRUE(linear)) {
    fc_col_name = "FC(lin)"
    temp_d <- d3_mod %>%
      mutate(across(-Label, gtools::logratio2foldchange)) %>%
      group_by(Label) %>%
      summarise(across(everything(), mean))
    (temp_d[2, -1] / temp_d[1, -1]) %>%
      pivot_longer(everything(), names_to = "variable", values_to = fc_col_name)
  } else {
    fc_col_name = "FC(log2)"
    temp_d <- d3_mod %>%
      group_by(Label) %>%
      summarise(across(everything(), mean))
    (temp_d[2, -1] - temp_d[1, -1]) %>%
      pivot_longer(everything(), names_to = "variable", values_to = fc_col_name)
  }
}

# export this function
#' @importFrom magrittr %>%
#' @import dplyr
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
    { purrr::reduce(list(., get_fc(d3_mod), get_fc(d3_mod, linear = FALSE)), left_join) } %>%
    select("variable", "pT", "BHT", "pW", "BHW", "FC(lin)", "FC(log2)")
}

# from visual_functions.R in generate_html_reports repo
# modified to incorporate the difference of MS data
deseq2_volcano <- function(res, df_filt, fdr = fdr, log2fc = log2fc, feature_col,
                           padj_col = NULL, log2fc_col = NULL) {
  if (!is.null(padj_col)) padj_col = padj_col else padj_col = "padj"
  if (!is.null(log2fc_col)) log2fc_col = log2fc_col else log2fc_col = "log2FoldChange"
  if (!"-log10padj" %in% colnames(res)) res <- res %>% mutate(`-log10padj` = -log10(padj_col))
  # data wrangling
  # # todo
  # # need to change the below line to a df with feature name as a column to
  # # better generalize across different data
  # dge_df <- as.data.frame(res) %>% rownames_to_column(var = "gene")
  dge_df <- res
  df_filt_up <- df_filt %>% filter(status == "Up")
  df_filt_down <- df_filt %>% filter(status == "Down")
  up_no = nrow(filter(df_filt, status == 'Up'))
  down_no = nrow(filter(df_filt, status == 'Down'))
  # truncate feature names if too long
  df_filt <- df_filt %>%
    mutate(across(variable, ~ str_trunc(., width = 15, ellipsis = "")))
  up_top10 <- df_filt_up %>% arrange(padj_col) %>% dplyr::slice(1:10)
  down_top10 <- df_filt_down %>% arrange(padj_col) %>% dplyr::slice(1:10)
  ###
  volcano_xlim <- max(na.omit(dge_df[log2fc_col]))
  volcano_ylim <- max(-log10(na.omit(dge_df[padj_col])))
  x_anno <- 0.7 * -volcano_xlim
  y_anno <- 0.9 * volcano_ylim
  ###
  feature_col = sym(feature_col)
  log2fc_col = sym(log2fc_col)
  # draw plot
  ggplot(dge_df, aes(!!log2fc_col, `-log10padj`)) +
    geom_point(alpha = 0.4, size = 1.5, colour = "grey50", na.rm = TRUE) +
    scale_x_continuous(limits = c(-volcano_xlim, volcano_xlim)) +
    scale_y_continuous(limits = c(0, volcano_ylim)) +
    geom_point(data = df_filt_down, shape = 21, alpha = 0.6,
               size = 1.5, fill = "blue", colour = "blue", na.rm = TRUE) +
    geom_point(data = df_filt_up, shape = 21, alpha = 0.6,
               size = 1.5, fill = "red", colour = "red", na.rm = TRUE) +
    ###
    geom_point(data = up_top10, shape = 21, fill = "red",
               colour = "black", size = 2, na.rm = TRUE) +
    geom_text_repel(data = up_top10, aes(label = !!feature_col), na.rm = TRUE,
                    size = 4, family = "Lucida Sans", max.overlaps = 20) +
    geom_point(data = down_top10, shape = 21, fill = "blue",
               colour = "black", size = 2, na.rm = TRUE) +
    geom_text_repel(data = down_top10, aes(label = !!feature_col), na.rm = TRUE,
                    size = 4, family = "Lucida Sans", max.overlaps = 20) +
    ###
    theme_bw(base_size = 14) +
    labs(x = "log2-fold change",
         y = "-log 10 (padj)",
         title = sprintf("FDR: %.2f, log2FC: %.1f", fdr, log2fc)) +
    theme(
      text = element_text(family = "Lucida Sans", face = "plain"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.title = element_text(size = 12.0),
      plot.title = element_text(size = 15.0, hjust = 0.5),
      aspect.ratio = 1
    ) +
    annotate("text", x = x_anno, y = y_anno, size = 5,
             label=sprintf("Up %i\nDown %i", up_no, down_no))
}

# from visual_functions.R in generate_html_reports repo
deseq2_hm <- function(transformed_count, df_filt, feature_col,
                      anno, top_n = NULL, col_order = NULL,
                      save = TRUE, padj_col = NULL) {
  if (!is.null(padj_col)) padj_col = padj_col else padj_col = "padj"
  cluster_cols <- TRUE
  # colors for hm and anno
  inc_col <- colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(200)
  temp_a <- sample(brewer.pal(8, "Dark2"), n_distinct(anno[1]))
  names(temp_a) <- levels(anno[[1]])
  anno_col[[names(anno)]] <- temp_a
  if(!is.null(top_n)) {
    top_n_df_filt <- df_filt %>% arrange(padj_col) %>% dplyr::slice(1:top_n)
    sub_mat <- transformed_count[top_n_df_filt[[feature_col]], ]
  } else {
    sub_mat <- df_filt %>%
      select(feature_col) %>%
      left_join(transformed_count) %>%
      column_to_rownames("variable")
      # transformed_count[df_filt[[feature_col]], ]
  }
  if (nrow(sub_mat) > 200) rowname_switch <- FALSE else rowname_switch <- TRUE
  if (ncol(sub_mat) > 60) {
    cellwidth = 10; fontsize = 6
  } else {
    cellwidth = 30; fontsize = 8}

  # draw hm
  if (!is.null(col_order)) { col_order = col_order; cluster_cols = FALSE;

    p3 <- pheatmap(
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
      main = paste("Clustering of normalized count (Z-score): ",
                   dim(sub_mat)[1], "features"))
  } else {

    p3 <- pheatmap(
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
      main = paste("Clustering of normalized count (Z-score): ",
                   dim(sub_mat)[1], "features"))
  }
  if (isTRUE(save)) ggsave("DE_features_heatmap_all_replicates.png", device = "png",
                           plot = p3, dpi = 300, width = 10, height = 10)
}

################


# first do data formatting
d3_mod <- t(d3) %>%
  as_tibble() %>%
  mutate(across(-Label, as.numeric)) %>%
  rename_with(str_trim)

# then do univariate analysis
uni_res <- do_univariate(d3_mod)


########### check with Bridget
# BHW are kinda different between my code and Bridget's
head(wt)
uni_res["P68871-M.VHLTPEEKSAVTAL.W" == uni_res$variable, ]
uni_res["P68871-L.WGKVNVDEVGGEALGRLL.V" == uni_res$variable, ]

uni_res["AAHLPAEFTPAV" == uni_res$variable, ]
uni_res["AAHLPAEFTPAVHAS" == uni_res$variable, ]

# > head(wt)
# pT         BHT           pW
# P68871-M.VHLTPEEKSAVTAL.W     0.0005579333 0.001556703 0.0003792417
# P68871-L.WGKVNVDEVGGEALGRLL.V 0.2217473194 0.306513114 0.2410024539
# P68871-L.WGKVNVDEVGGEAL.G     0.0003873439 0.001101128 0.0009007453
# P68871-V.AGVANALAHKYH         0.0094188851 0.019529065 0.0107771150
# P68871-D.GLAHLDNLKGTF.A       0.0032427103 0.008097818 0.0001870327
# P68871-V.HLTPEEK.S            0.4059796106 0.494616668 0.6151583083
# BHW   FC(lin)
# P68871-M.VHLTPEEKSAVTAL.W     0.004343393  6.187240
# P68871-L.WGKVNVDEVGGEALGRLL.V 0.419922624  2.391418
# P68871-L.WGKVNVDEVGGEAL.G     0.003130252  3.638086
# P68871-V.AGVANALAHKYH         0.040491455  1.542126
# P68871-D.GLAHLDNLKGTF.A       0.020222175  2.103264
# P68871-V.HLTPEEK.S            0.602605752 -1.169379
# > uni_res["P68871-M.VHLTPEEKSAVTAL.W" == uni_res$variable, ]
# # A tibble: 1 x 6
# variable                        pT     BHT       pW     BHW `FC(lin)`
# <chr>                        <dbl>   <dbl>    <dbl>   <dbl>     <dbl>
#   1 P68871-M.VHLTPEEKSAVTAL.W 0.000558 0.00156 0.000379 0.00103      6.19
# > uni_res["P68871-L.WGKVNVDEVGGEALGRLL.V" == uni_res$variable, ]
# # A tibble: 1 x 6
# variable                         pT   BHT    pW   BHW `FC(lin)`
# <chr>                         <dbl> <dbl> <dbl> <dbl>     <dbl>
#   1 P68871-L.WGKVNVDEVGGEALGRLL.V 0.222 0.307 0.241 0.325      2.39


# visualizations
fdr = 0.05
log2fc = 0

# format univariate results tibble for plotting
uni_res <- uni_res %>%
  rowwise() %>%
  mutate(padj = min(c(BHT, BHW))) %>% # get the lowest padj
  ungroup() %>%
  mutate(`-log10padj` = -log(padj))
  # # temp for now for testing, remove the outlier peptide
  # filter(`-log10padj` < 130)

# get DE features only tibble
uni_res_filt <- uni_res %>%
  filter(BHT < fdr | BHW < fdr) %>%
  mutate(status = if_else("FC(log2)" < 0, "Down", "Up"))


# 1) volcano plot
deseq2_volcano(uni_res, uni_res_filt, fdr = fdr, log2fc = log2fc,
               "variable", padj_col = "padj", log2fc_col = "FC(log2)")

# prep the annotation for hm
anno <- data.frame(Label = as.factor(t(d3)[, "Label"]))

# prep the transformed matrix for hm
d1 <- d1 %>%
  rownames_to_column("variable") %>%
  mutate(across(-variable, as.numeric))

# sanity check the variable names matched
all(uni_res_filt$variable %in% d1$variable)

# 2) heatmap/multiple heatmaps?
# build user_supplied_col_order into the function, default to colnames()
deseq2_hm(d1, uni_res_filt, "variable",
          anno, top_n = NULL, col_order = NULL,
          save = FALSE, padj_col = NULL)

# 3) venn diagram between t-test and wilcox test
