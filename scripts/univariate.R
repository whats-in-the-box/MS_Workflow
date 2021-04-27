library(tidyverse)

##### Univariate Statistics and Fold Change

save(d3, d3_mod, wt, file = "for_univariate.Rda")

# declare function
# only works for pairwise foldchange
get_fc_linear <- function(d3) {
  temp_d <- d3_mod %>%
    mutate(across(-Label, gtools::logratio2foldchange)) %>%
    group_by(Label) %>%
    summarise(across(everything(), mean))
    # # works for now
    # # it's really clunky from here onward, don't know how to fix it
    # pivot_longer(-Label) %>%
    # group_by(name) %>%
    # arrange(name) %>%
    # summarise('FC(lin)' = value/lag(value)) %>%
    # drop_na('FC(lin)')
  # tried a base R approach instead
  (temp_d[2, -1]/temp_d[1, -1]) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "FC(lin)")
}

do_univariate <- function(d3_mod) {
  d3_mod %>%
    pivot_longer(cols = -Label, names_to = "variable", values_to = "value") %>%
    group_nest(variable) %>%
    mutate(pT = map(data, ~t.test(value ~ Label, data = .x)$p.value),
           pW = map(data, ~wilcox.test(value ~ Label, data = .x)$p.value)) %>%
    unnest(cols = c("pT", "pW")) %>%
    mutate(BHT = p.adjust(pT, method = "BH"),
           BHW = p.adjust(pW, method = "BH")) %>%
    select(-data) %>%
    left_join(., get_fc_linear(d3_mod)) %>%
    select("variable", "pT", "BHT", "pW", "BHW", "FC(lin)")
}
################

# first do data formatting
d3_mod <- t(d3) %>%
  as_tibble() %>%
  mutate(across(-Label, as.numeric))  %>%
  rename_with(str_trim)

# then do univariate analysis
uni_res <- do_univariate(d3_mod)



uni_res["P68871-M.VHLTPEEKSAVTAL.W" == uni_res$variable, ]
uni_res["P68871-L.WGKVNVDEVGGEALGRLL.V" == uni_res$variable, ]
# BHW are kinda different between my method and Bridget's

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
