Differential expression analysis proteomics
================
Kaspar Bresser

- [Import and tidy data](#import-and-tidy-data)
- [number of proteins](#number-of-proteins)
- [Abundance LFQ](#abundance-lfq)
- [PCA](#pca)
- [Test differences](#test-differences)
- [Plotting](#plotting)

Pre-processing of the MS data from the activation time-course.

First load packages

``` r
library(tidyverse)
library(broom) 
library(limma)
library(ggrepel)
library(readxl)
library(circlize)
library(ComplexHeatmap)
library(cluster)
library(lemon)
#library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(gghighlight)
library(ggpubr)
library(ggrastr)
library(ggridges)
library(corrplot)
```

### Import and tidy data

Selecting the proteins that have a “proteotypic” count of more than 1.

``` r
read_tsv("Data/20240625_Timecourse_report.pr_matrix.tsv") %>% 
  dplyr::select(Protein.Group, Genes, Proteotypic, Precursor.Id, matches('\\.d$')) %>%
    gather('Sample', 'Intensity', matches('\\.d$')) %>%
    filter(!is.na(Intensity)) %>%
    filter(Proteotypic == T) %>% #
    group_by(Protein.Group) %>%
    distinct(Precursor.Id, .keep_all = T) %>%
    summarise(count = n()) %>%
    filter(count > 1)  %>%
    dplyr::pull(Protein.Group) -> select.prots

str(select.prots)
```

    ##  chr [1:6777] "A0A075B6T6" "A0A0B4J245" "A0A0U1RRE5" "A0AVK6" "A0AVT1" ...

import proteomic quantification data

``` r
dat <- read_tsv("Data/20240625_Timecourse_report.pg_matrix.tsv")

dat
```

    ## # A tibble: 8,037 × 19
    ##    Protein.Group     Protein.Ids      Protein.Names Genes First.Protein.Descri…¹
    ##    <chr>             <chr>            <chr>         <chr> <chr>                 
    ##  1 A0A024RBG1;Q9NZJ9 Q9NZJ9;Q8NFP7;Q… NUD4B_HUMAN;… NUDT… Diphosphoinositol pol…
    ##  2 A0A075B6T6        A0A075B6T6       TVAL2_HUMAN   TRAV… T cell receptor alpha…
    ##  3 A0A0B4J241        A0A0B4J241       TVAM1_HUMAN   TRAV… T cell receptor alpha…
    ##  4 A0A0B4J245        A0A0B4J245       TVAL1_HUMAN   TRAV… T cell receptor alpha…
    ##  5 A0A0B4J279        A0A0B4J279       TVA21_HUMAN   TRAV… T cell receptor alpha…
    ##  6 A0A0B4J2D5;P0DPI2 P0DPI2;A0A0B4J2… GAL3A_HUMAN;… GATD… Glutamine amidotransf…
    ##  7 A0A0B4J2E0;P01733 P01733;A0A0B4J2… TVBL3_HUMAN;… TRBV… T cell receptor beta …
    ##  8 A0A0B4J2F0        A0A0B4J2F0       PIOS1_HUMAN   PIGB… Protein PIGBOS1       
    ##  9 A0A0B4J2F2;P57059 A0A0B4J2F2;P570… SIK1B_HUMAN;… SIK1… Probable serine/threo…
    ## 10 A0A0U1RRE5        A0A0U1RRE5       NBDY_HUMAN    NBDY  Negative regulator of…
    ## # ℹ 8,027 more rows
    ## # ℹ abbreviated name: ¹​First.Protein.Description
    ## # ℹ 14 more variables:
    ## #   `C:\\Users\\massspecuser\\Desktop\\20240622_Carmen\\RES00360_Anouk_Jurgens\\20240620_AnoukJurgens_Timecourse\\1_0h_Donor1_2_S2-A4_1_9368.d` <dbl>,
    ## #   `C:\\Users\\massspecuser\\Desktop\\20240622_Carmen\\RES00360_Anouk_Jurgens\\20240620_AnoukJurgens_Timecourse\\2_2h_Donor1_S2-B5_1_9369.d` <dbl>,
    ## #   `C:\\Users\\massspecuser\\Desktop\\20240622_Carmen\\RES00360_Anouk_Jurgens\\20240620_AnoukJurgens_Timecourse\\3_4h_Donor1_S2-C5_1_9370.d` <dbl>,
    ## #   `C:\\Users\\massspecuser\\Desktop\\20240622_Carmen\\RES00360_Anouk_Jurgens\\20240620_AnoukJurgens_Timecourse\\4_6h_Donor1_S2-D5_1_9371.d` <dbl>, …

Tidy up a bit, and convert abundances to log2

``` r
dat %>% 
  filter(Protein.Group %in% select.prots) %>% 
  pivot_longer(contains("Users"), names_to = "sample", values_to = "LFQ") %>% 
  mutate(sample = str_extract(sample, "\\d+_\\d+h_.+onor\\d+")) %>% 
  separate(sample, into = c("sample.nr", "timepoint",  "donor")) %>% 
  mutate(donor = str_to_lower(donor)) %>% 
  mutate(LFQ = log2(LFQ)) %>% 
  mutate(timepoint = paste0("T", timepoint)) %>%   
  mutate(timepoint = case_when(timepoint == "T16h" ~ "T24h", TRUE ~ timepoint)) %>% 
  mutate(timepoint = factor(timepoint, levels = c("T0h", "T2h", "T4h", "T6h", "T24h"))) %>% 
  dplyr::rename(gene.name = Genes) %>% 
  dplyr::select(-c(Protein.Group, Protein.Names, First.Protein.Description)) -> dat

dat
```

    ## # A tibble: 94,878 × 6
    ##    Protein.Ids gene.name sample.nr timepoint donor    LFQ
    ##    <chr>       <chr>     <chr>     <fct>     <chr>  <dbl>
    ##  1 A0A075B6T6  TRAV12-2  1         T0h       donor1  11.7
    ##  2 A0A075B6T6  TRAV12-2  2         T2h       donor1  11.3
    ##  3 A0A075B6T6  TRAV12-2  3         T4h       donor1  11.2
    ##  4 A0A075B6T6  TRAV12-2  4         T6h       donor1  10.7
    ##  5 A0A075B6T6  TRAV12-2  5         T24h      donor1  11.1
    ##  6 A0A075B6T6  TRAV12-2  6         T0h       donor2  11.2
    ##  7 A0A075B6T6  TRAV12-2  7         T2h       donor2  11.5
    ##  8 A0A075B6T6  TRAV12-2  8         T4h       donor2  11.4
    ##  9 A0A075B6T6  TRAV12-2  9         T6h       donor2  11.4
    ## 10 A0A075B6T6  TRAV12-2  10        T24h      donor2  11.3
    ## # ℹ 94,868 more rows

### number of proteins

``` r
library(RColorBrewer)
library(ggh4x)

dat %>% 
  na.omit() %>% 
  count(timepoint, donor) %>% 
  ggplot(aes(x = fct_cross(timepoint, donor), y = n, fill = timepoint))+
    geom_bar(stat = "identity", color = "black")+
    scale_fill_manual(values = brewer.pal(5, "Purples"))+
    theme_classic()+
    theme(panel.grid.major.y = element_line())+
    facet_grid2(~timepoint, scales = "free_x", space = "free_x")
```

<img src="7-Supplementary_protein_files/figure-gfm/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />

``` r
ggsave("Figs/Supplementary1/protein_detected.pdf", width = 5.4, height = 2.3)
```

### Abundance LFQ

``` r
dat %>% 
  ggplot(aes(x = fct_cross(timepoint, donor), y = LFQ, fill = timepoint))+
    geom_boxplot( color = "black")+
    scale_fill_manual(values = brewer.pal(5, "Purples"))+
    theme_classic()+
    theme(panel.grid.major.y = element_line())+
    facet_grid2(~timepoint, scales = "free_x", space = "free_x")
```

    ## Warning: Removed 569 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

<img src="7-Supplementary_protein_files/figure-gfm/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

``` r
ggsave("Figs/Supplementary1/protein_abundance.pdf", width = 5.4, height = 2.3)
```

    ## Warning: Removed 569 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

### PCA

``` r
dat %>% 
  group_by(gene.name, timepoint, donor) %>% 
  summarise(LFQ = median(LFQ)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = c(donor, timepoint), values_from = LFQ) %>% 
  na.omit() %>% 
  column_to_rownames("gene.name") -> dat.pca
```

``` r
dat.pca %>% 
  t() %>% 
  as.data.frame() %>% 
  prcomp(center = T, scale = F) -> pca
```

``` r
library(RColorBrewer)
as.data.frame(pca$x) %>% 
  rownames_to_column("sample") %>% 
  separate(sample, into = c("donor", "timepoint")) %>% 
  mutate(timepoint = factor(timepoint, levels = c("T0h", "T2h", "T4h", "T6h", "T24h"))) %>% 
ggplot(aes(x = PC2, y = PC1, fill = timepoint, shape = donor))+
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24))+
  scale_fill_manual(values = brewer.pal(5, "Purples"))+
  ggtitle("pca of LC-MS data")+
  theme_classic()+
  theme(panel.grid.major = element_line())+
  guides(fill = guide_legend("timepoint", override.aes = list(shape = 21)))
```

<img src="7-Supplementary_protein_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

``` r
ggsave("Figs/Supplementary1/protein_pca.pdf", width = 3.5, height = 2.5)
```

### Test differences

will test for differences using limma

first set up a design matrix, we want to compare between time-points

``` r
dat %>% 
#  mutate(sample.nr = paste0("sample.", sample.nr)) %>% 
  dplyr::select(sample.nr, timepoint, donor) %>% 
  distinct() -> pheno

design <- model.matrix(~0+ pheno$timepoint)
colnames(design) <- unique(pheno$timepoint)

design
```

    ##    T0h T2h T4h T6h T24h
    ## 1    1   0   0   0    0
    ## 2    0   1   0   0    0
    ## 3    0   0   1   0    0
    ## 4    0   0   0   1    0
    ## 5    0   0   0   0    1
    ## 6    1   0   0   0    0
    ## 7    0   1   0   0    0
    ## 8    0   0   1   0    0
    ## 9    0   0   0   1    0
    ## 10   0   0   0   0    1
    ## 11   0   1   0   0    0
    ## 12   0   0   1   0    0
    ## 13   0   0   0   1    0
    ## 14   0   0   0   0    1
    ## attr(,"assign")
    ## [1] 1 1 1 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$`pheno$timepoint`
    ## [1] "contr.treatment"

Get the LFQ values for all samples and fit linear model.

lmFit needs samples to be columns, and genes to be row names.

``` r
dat %>% 
  pivot_wider(names_from = c(sample.nr, timepoint, donor), values_from = LFQ) %>% 
  na.omit() %>% 
  distinct(gene.name, .keep_all = T) %>% 
  column_to_rownames("gene.name") %>% 
  dplyr::select(contains("donor")) %>% 
  lmFit(., design) -> lm.fit

design
```

    ##    T0h T2h T4h T6h T24h
    ## 1    1   0   0   0    0
    ## 2    0   1   0   0    0
    ## 3    0   0   1   0    0
    ## 4    0   0   0   1    0
    ## 5    0   0   0   0    1
    ## 6    1   0   0   0    0
    ## 7    0   1   0   0    0
    ## 8    0   0   1   0    0
    ## 9    0   0   0   1    0
    ## 10   0   0   0   0    1
    ## 11   0   1   0   0    0
    ## 12   0   0   1   0    0
    ## 13   0   0   0   1    0
    ## 14   0   0   0   0    1
    ## attr(,"assign")
    ## [1] 1 1 1 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$`pheno$timepoint`
    ## [1] "contr.treatment"

Make a contrast matrix for the comparisons of interest, will compare
everything against 0h and everything against 24h.

``` r
contrast.matrix <- makeContrasts(
  
  T2h - T0h,
  T4h - T0h,
  T6h - T0h,
  T24h - T0h,
  T2h - T24h,
  T4h - T24h,
  T6h - T24h,

  
  levels = colnames(design))

contrast.matrix
```

    ##       Contrasts
    ## Levels T2h - T0h T4h - T0h T6h - T0h T24h - T0h T2h - T24h T4h - T24h
    ##   T0h         -1        -1        -1         -1          0          0
    ##   T2h          1         0         0          0          1          0
    ##   T4h          0         1         0          0          0          1
    ##   T6h          0         0         1          0          0          0
    ##   T24h         0         0         0          1         -1         -1
    ##       Contrasts
    ## Levels T6h - T24h
    ##   T0h           0
    ##   T2h           0
    ##   T4h           0
    ##   T6h           1
    ##   T24h         -1

Compute differential expression

``` r
# fit contrasts
lm.fit2 <- contrasts.fit(lm.fit, contrast.matrix)
# ebayes
lm.fit2 <- eBayes(lm.fit2)
```

Check number of significant DE proteins

``` r
results <- decideTests(lm.fit2, p.value = 0.05, adjust.method = 'BH', lfc = 0.5)
summary(results)
```

    ##        T2h - T0h T4h - T0h T6h - T0h T24h - T0h T2h - T24h T4h - T24h
    ## Down          96        82       134        666        522        227
    ## NotSig      6321      6325      6255       5406       5452       5828
    ## Up            93       103       121        438        536        455
    ##        T6h - T24h
    ## Down          158
    ## NotSig       5982
    ## Up            370

Set function to extract all results

``` r
get_DE <- function(fit, comp){
  
  topTable(fit, number = Inf, coef = comp, sort.by = 'none') %>% 
    mutate(comparison = comp)
}
```

Extract all results

``` r
comparisons <- colnames(lm.fit2$coefficients)

comparisons %>% 
  map( ~get_DE(lm.fit2, .)) %>% 
  map(rownames_to_column, "gene.name") %>% 
  list_rbind() %>% 
  as_tibble() %>% 
  dplyr::rename(logFC.protein = logFC, adj.P.Val.protein = adj.P.Val) -> all.DE.protein

all.DE.protein
```

    ## # A tibble: 45,570 × 8
    ##    gene.name logFC.protein AveExpr       t P.Value adj.P.Val.protein     B
    ##    <chr>             <dbl>   <dbl>   <dbl>   <dbl>             <dbl> <dbl>
    ##  1 TRAV12-2        -0.148     11.2 -0.691  0.503               0.771 -6.59
    ##  2 NBDY            -0.753     10.7 -3.08   0.00960             0.122 -3.12
    ##  3 UBA6            -0.125     14.8 -1.82   0.0934              0.385 -5.27
    ##  4 ESYT2            0.0606    15.5  0.738  0.475               0.755 -6.56
    ##  5 MED19            0.203     12.1  1.49   0.161               0.478 -5.75
    ##  6 UHRF1BP1L        0.136     11.6  0.829  0.423               0.719 -6.48
    ##  7 TMEM223          0.0115    12.2  0.0941 0.927               0.973 -6.84
    ##  8 ARHGAP10         0.378     10.2  1.54   0.150               0.467 -5.69
    ##  9 ILVBL            0.195     15.0  2.79   0.0163              0.166 -3.64
    ## 10 ZC3H12D         -0.651     10.3 -2.42   0.0325              0.236 -4.30
    ## # ℹ 45,560 more rows
    ## # ℹ 1 more variable: comparison <chr>

### Plotting

Plot example proteins

rasterise(geom_point(), dpi = 300, scale = 0.5)

Function to plot volcano’s with counts

``` r
plt_volcanos <- function(dat, title){
  
  dat %>% 
    group_by(col) %>%
    count(adj.P.Val.protein < 0.05) %>% 
    pull(n) %>% 
    as.character() -> cnts
  
  ggplot(dat, aes( x = logFC.protein, y = log10(adj.P.Val.protein)))+ 
  rasterise(geom_point(aes(fill = col), shape = 21, size = 2), dpi = 300, scale = 1.0)+
  scale_fill_manual(values = c( "#4796C5","#C23A3B"))+
  scale_y_reverse()+
  labs(title = "")+
  gghighlight(adj.P.Val.protein < 0.05 & abs(logFC.protein) > 0.5,  
              unhighlighted_params = list(shape = 21, fill = "grey"),  use_group_by = F)+
  theme_classic()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line())+
    xlim(-4,4)+
    ylim(0, -8)+
  geom_vline(xintercept = c(0.5, -0.5), linetype ="dotted")+
  geom_hline(yintercept = log10(0.05), linetype ="dotted")+
  annotate(geom = "text", label = cnts[2], x=-3.5, y=-8, color = "#4796C5", fontface =2)+
  annotate(geom = "text", label =cnts[4], x=3.5, y=-8, color = "#C23A3B", fontface =2)+
  ggtitle(title)

}
```

Plot volcano’s with counts

``` r
all.DE.protein %>% 
  filter(str_detect(comparison, "0h$")) %>% 
  mutate(comparison = fct_relevel(comparison, "T2h - T0h", "T4h - T0h", "T6h - T0h"),
         comparison = fct_recode(comparison, T2h = "T2h - T0h", T4h = "T4h - T0h", T6h = "T6h - T0h", T24h = "T24h - T0h"),
         col = case_when(logFC.protein > 0 ~ "up", TRUE ~ "down")) -> for.plot


for.plot %>% 
  nest(data = -comparison) %>% 
  mutate(plots = map2(data, comparison, plt_volcanos)) %>% 
  pull(plots) %>% 
  ggpubr::ggarrange(plotlist = . , ncol = 4)
```

<img src="7-Supplementary_protein_files/figure-gfm/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

``` r
ggsave("Figs/Supplementary1/protein_volcanos_DE.pdf", width = 12, height = 4)
```

``` r
all.DE.protein %>% 
  filter(adj.P.Val.protein < 0.05) %>% 
  pull(gene.name) -> genes


dat %>% 
  group_by(gene.name, timepoint, donor) %>% 
  summarise(LFQ = median(LFQ)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = c(donor, timepoint), values_from = LFQ) %>% 
  na.omit() %>% 
  filter(gene.name %in% genes) %>% 
  column_to_rownames("gene.name") %>% 
  t() %>%
  scale() %>%
  t() -> dat.for.cluster
```

``` r
### heatmap
Heatmap(dat.for.cluster,
    cluster_columns = F,
    cluster_column_slices = F,
  #  clustering_distance_columns = "euclidean",
  use_raster = T,
  raster_quality = 10, 
  row_km = 4,
    show_row_names = F,
    show_row_dend = T,
    show_column_dend = T,
    show_column_names = T,
    # col = colorRamp2(c(-18, 0, 8),c('blue', 'white', 'red')),
 #   column_split = pheno$ct,
    row_names_gp = grid::gpar(fontsize = 8))#,cluster_columns = T, km = 3)
```

<img src="7-Supplementary_protein_files/figure-gfm/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

``` r
###350x350
```
