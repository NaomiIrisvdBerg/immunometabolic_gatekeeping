##Single cell RNAseq analysis publicly available CELLxGENE data across all tissues for human, normal
##In this script, we first explore the inter Sex sources of variance, in order to determine 
#whether the default download of CELLxGENE (one aggregated Expression value for each unique combination of Cell Type, Gene Symbol and Tissue)
#can reliably be used to investigate cross-tissue expression patterns for a given subset of genes and cell types.
##We then analyse cross-tissue patterns including solid tissues other than those listed in the manuscript, to explore robustness of the patterns
##From line ~708 onwards, we explore trends exclusively in the 14 focus tissues
##run script chronologically due to reuse of object names


#load libraries:
library(readr) ##we ran this script with 2.1.5
library(dplyr) ##v 1.1.4
library(tidyr) ##v 1.3.1
library(ggplot2) ##v 4.0.0
library(stringr) ##v 1.5.2
library(data.table) ##v 1.17.8
library(patchwork) ##v 1.3.2
library(pheatmap) ##v 1.0.13
library(randomcoloR) ##v 1.1.0.1
library(RColorBrewer) ##v 1.1.3
library(lme4) ##v 1.1.37


CELLxGENE_gene_expression_020626_by_publication <- read_csv("Desktop/CELLxGENE_Feb26/CELLxGENE_gene_expression_020626.csv", skip = 9)
CELLxGENE_gene_expression_020626_by_sex <- read_csv("Desktop/CELLxGENE_Feb26/CELLxGENE_gene_expression_020626 (1).csv", skip = 9)
CELLxGENE_gene_expression_020626 <- read_csv("Desktop/CELLxGENE_Feb26/CELLxGENE_gene_expression_020626 (2).csv", skip = 9)


##publication info is collapsed into cell type, (we lost hierarchical structuring when downloading as csv), so need to recover that first:

clean_cellxgene <- function(df) {
  pub_regex <- "\\([12][0-9]{3}\\)"  ##e.g. "(2020)"
  
  df %>%
    dplyr::mutate(
      is_pub_row      = stringr::str_detect(`Cell Type`, pub_regex) | `Cell Type` %in% c("No Publication"),
      is_placeholder  = `Cell Type` %in% c("aggregated", "cell"),
      cell_type_anchor = if_else(is_pub_row | is_placeholder, NA_character_, `Cell Type`)
    ) %>%
    dplyr::group_by(Tissue) %>%
    fill(cell_type_anchor, .direction = "down") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Publication = if_else(is_pub_row, `Cell Type`, Publication),
      `Cell Type` = cell_type_anchor
    ) %>%
    dplyr::filter(!is.na(`Cell Type`))
}



datasets <- list(
  CELLxGENE_gene_expression_020626_by_publication
)

names(datasets) <- c(
  "Feb_2026_by_pub"
)

cleaned_datasets <- lapply(datasets, clean_cellxgene)
list2env(cleaned_datasets, envir = .GlobalEnv)



Feb_2026_by_pub




##We exclude: non-solid tissues (circulating or fluid-derived), embryonic or developmental tissues, or tissue subcompartments rather than whole solid organs
Feb_2026_by_pub_solid_tissues <- subset(Feb_2026_by_pub, Tissue != "blood" 
                                                                        & Tissue != "bone marrow" & Tissue != "lymph node" 
                                                                        & Tissue != "milk" & Tissue != "pleural fluid" & Tissue != "saliva" & Tissue != "yolk sac" & Tissue != "embryo" 
                                                                        & Tissue != "placenta" & Tissue != "vasculature" & Tissue != "mucosa" & Tissue != "lamina propria" & Tissue != "omentum")


Feb_2026_by_pub_solid_tissues_a <- subset(Feb_2026_by_pub_solid_tissues, Publication != "aggregated")
Feb_2026_by_pub_solid_tissues_a <- subset(Feb_2026_by_pub_solid_tissues_a, `Number of Cells Expressing Genes` >=10)
min(Feb_2026_by_pub_solid_tissues_a$`Number of Cells Expressing Genes`)



##investigate number of publications per tissue:
check_publications <- subset(Feb_2026_by_pub_solid_tissues_a, `Cell Type` == "endothelial cell" | `Cell Type` == "fibroblast" | `Cell Type` == "T cell") %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(n_publications = n_distinct(Publication), .groups = "drop")

##by cell type:
check_publications1 <- Feb_2026_by_pub_solid_tissues_a %>%
  dplyr::filter(`Cell Type` %in% c("endothelial cell",
                                   "fibroblast",
                                   "T cell")) %>%
  dplyr::group_by(Tissue, `Cell Type`) %>%
  dplyr::summarise(
    n_publications = dplyr::n_distinct(Publication),
    .groups = "drop"
  )



#need to an exception for urinary bladder as we combine this with "bladder organ" as a representative for "bladder", to arrive at a total of 3 datasets:
at_least_2_pubs <- subset(
  check_publications,
  n_publications > 1 | Tissue == "urinary bladder"
)

Feb_2026_by_pub_solid_tissues_b <- Feb_2026_by_pub_solid_tissues_a %>%
  dplyr::filter(Tissue %in% unique(at_least_2_pubs$Tissue))


unique(Feb_2026_by_pub_solid_tissues_b$Tissue)
unique(Feb_2026_by_pub_solid_tissues_b$`Cell Type`)






###Explore sex-driven variance (i.e., is aggregation by sex as done by default otherwise alright?):
clean_cellxgene_sex <- function(df) {
  df %>%
    dplyr::mutate(
      ##rows where Cell Type is actually a sex label
      is_sex_row = `Cell Type` %in% c("male", "female", "unknown"),
      ##placeholder rows that just say "aggregated" or "cell"
      is_placeholder = `Cell Type` %in% c("aggregated", "cell"),
      ##last real cell type to carry downward
      cell_type_anchor = if_else(is_sex_row | is_placeholder,
                                 NA_character_, `Cell Type`)
    ) %>%
    dplyr::group_by(Tissue) %>%
    fill(cell_type_anchor, .direction = "down") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      `Cell Type` = cell_type_anchor
    ) %>%
    filter(!is.na(`Cell Type`))
}

datasets <- list(
  CELLxGENE_gene_expression_020626_by_sex
)

names(datasets) <- c(
  "Feb_2026_by_sex"
)

cleaned_datasets_sex <- lapply(datasets, clean_cellxgene_sex)
list2env(cleaned_datasets_sex, envir = .GlobalEnv)



head(Feb_2026_by_sex)


##take out non-solid/developmental/non-specific tissues like before:
Feb_2026_by_sex_solid_tissues <- subset(Feb_2026_by_sex, Tissue != "blood"
                                                                & Tissue != "bone marrow" & Tissue != "lymph node"
                                                                & Tissue != "milk" & Tissue != "pleural fluid" & Tissue != "saliva" & Tissue != "yolk sac" & Tissue != "embryo"
                                                                & Tissue != "placenta" & Tissue != "vasculature" & Tissue != "mucosa" & Tissue != "lamina propria" & Tissue != "omentum")


Feb_2026_by_sex_solid_tissues_a <- subset(Feb_2026_by_sex_solid_tissues, Sex != "aggregated")
Feb_2026_by_sex_solid_tissues_a <- subset(Feb_2026_by_sex_solid_tissues_a, `Number of Cells Expressing Genes` >= 10)
min(Feb_2026_by_sex_solid_tissues_a$`Number of Cells Expressing Genes`)

#View(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_a)


Feb_2026_by_sex_solid_tissues_b <- Feb_2026_by_sex_solid_tissues_a %>%
  dplyr::filter(Tissue %in% unique(at_least_2_pubs$Tissue))

Feb_2026_by_sex_solid_tissues_b1 <- 
  Feb_2026_by_sex_solid_tissues_b %>%
  dplyr::mutate(Tissue = case_when(
    Tissue %in% c("urinary bladder", "bladder organ") ~ "bladder",
    TRUE ~ Tissue
  ))

unique(Feb_2026_by_sex_solid_tissues_b1$Tissue)
unique(Feb_2026_by_sex_solid_tissues_b1$`Cell Type`)


Feb_2026_by_sex_solid_tissues_b_endothelial_cell <- subset(Feb_2026_by_sex_solid_tissues_b1, `Cell Type` == "endothelial cell")

Feb_2026_by_sex_solid_tissues_b_fibroblast <- subset(Feb_2026_by_sex_solid_tissues_b1, `Cell Type` == "fibroblast")


##inspect sex-based differences
ggplot(
  subset(Feb_2026_by_sex_solid_tissues_b_fibroblast, Sex != "unknown"),
  aes(x = Expression, colour = Sex, fill = Sex)
) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ Tissue, scales = "free") + ##one panel per tissue
  xlim(0.5, 3.5) +
  labs(
    x = "Expression",
    y = "Density",
    title = "Distribution of Expression per Sex by Tissue"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9)
  )




####investigate role of sex variance for stromal cells across all solid tissues, excluding T cell exhaustion genes
Feb_2026_by_sex_solid_tissues_b1_stromal <- subset(Feb_2026_by_sex_solid_tissues_b1, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")
Feb_2026_by_sex_solid_tissues_b1_stromal <- subset(Feb_2026_by_sex_solid_tissues_b1_stromal, `Gene Symbol` != "TOX" 
                                                                           & `Gene Symbol` != "LAG3" & `Gene Symbol` != "PDCD1" & `Gene Symbol` != "CTLA4" & `Gene Symbol` != "HAVCR2")


##take out sex-specific tissues:
Feb_2026_by_sex_solid_tissues_b1_stromal <- subset(Feb_2026_by_sex_solid_tissues_b1_stromal, Tissue != "prostate gland" & Tissue != "breast" & Tissue != "testis"
                                                                           & Tissue != "fallopian tube" & Tissue != "uterus" & Tissue != "ovary")

##sex differences are real but smaller than biological differences across genes or tissues:
fit_stromal_cells_all_solid_tissues_by_sex <- lmer(Expression ~ 1 + (1 | Sex) + (1 | `Gene Symbol`) + (1 | Tissue),
                                                   data = Feb_2026_by_sex_solid_tissues_b1_stromal)

VarCorr(fit_stromal_cells_all_solid_tissues_by_sex)



####investigate role of sex-mediated variance for T cells across all solid tissues, excluding T cell exhaustion genes
Feb_2026_by_sex_solid_tissues_b1_T <- subset(Feb_2026_by_sex_solid_tissues_b1, `Cell Type` == "T cell")

##take out sex-specific tissues:
Feb_2026_by_sex_solid_tissues_b1_T <- subset(Feb_2026_by_sex_solid_tissues_b1_T, Tissue != "prostate gland" & Tissue != "breast" & Tissue != "testis"
                                                                     & Tissue != "fallopian tube" & Tissue != "uterus" & Tissue != "ovary")


##Variance by sex is real but much smaller than biological differences across genes and/or tissues.
fit_T_cells_all_solid_tissues_by_sex <-  lmer(Expression ~ 1 + (1 | Sex) + (1 | `Gene Symbol`) + (1 | Tissue),
                                              data = Feb_2026_by_sex_solid_tissues_b1_T)

VarCorr(fit_T_cells_all_solid_tissues_by_sex)
##hence, it is asssumed safe to aggregate each combination of Tissue, Cell Type and Gene Symbol by sex.
#more specifically: under the conditions tested here (adult non–sex-specific solid tissues, stromal cells and T cells), sex-driven variance is substantially smaller than gene- or tissue-driven variance, supporting the use of sex-aggregated CELLxGENE expression values for cross-tissue analyses
###this is the default by which CELLxGENE data is downloaded when selecting genes, tissues and cell types of interest.
##note that variance by publication is difficult to assess accordingly since that is indistinguishable from exploring variance by Tissue, 
#since most Publications included contain data for a single (or limited set of) Tissue. So since publication effects cannot be disentangled 
#from tissue effects in the default aggregated download, they are addressed explicitly in the second script in the Github Repository (260212_cellxgene_publication_aware_aggregation.R) using publication-aware aggregation






###inspect relative expressions in the aggregated data (DEFAULT DOWNLOAD FROM CELLXGENE):
##Note that in the script below we do not just calculate the relative (weighted & normalised) expressions as shown in Table 1 of the manuscript,
##we also explore cross-tissue correlations for solid tissues beyond those listed in Table 1.
##for simple relative ordering of Expressions as shown in Table 1, the block at line ~670 contains a helpful visualisation



CELLxGENE_gene_expression_020626

CELLxGENE_gene_expression_020626_solid_tissues <- subset(CELLxGENE_gene_expression_020626, Tissue != "blood" 
                                                         & Tissue != "bone marrow" & Tissue != "lymph node" 
                                                         & Tissue != "milk" & Tissue != "pleural fluid" & Tissue != "saliva" & Tissue != "yolk sac" & Tissue != "embryo" 
                                                         & Tissue != "placenta" & Tissue != "vasculature" & Tissue != "mucosa" & Tissue != "lamina propria" & Tissue != "omentum")

CELLxGENE_gene_expression_020626_solid_tissues_a <- subset(CELLxGENE_gene_expression_020626_solid_tissues, `Cell Type` != "aggregated")

##QC 1: enough publications
CELLxGENE_gene_expression_020626_solid_tissues_b <- CELLxGENE_gene_expression_020626_solid_tissues_a %>%
  dplyr::filter(Tissue %in% unique(at_least_2_pubs$Tissue))



unique(CELLxGENE_gene_expression_020626_solid_tissues_b$Tissue)



##QC 2: at least 10 Cells Expressing Genes:
CELLxGENE_gene_expression_020626_solid_tissues_b1 <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b, `Number of Cells Expressing Genes` >= 10)
min(CELLxGENE_gene_expression_020626_solid_tissues_b1$`Number of Cells Expressing Genes`)

t_cell_data <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b1, `Cell Type` == "T cell")

##group urinary bladder and bladder organ together:
t_cell_data <- t_cell_data %>%
  dplyr::mutate(
    Tissue = case_when(
      Tissue %in% c("bladder organ", "urinary bladder") ~ "bladder",
      TRUE                                              ~ Tissue
    )
  )

weighted_expression <- t_cell_data %>%
  dplyr::filter(!is.na(Expression), !is.na(`Cell Count`)) %>%
  dplyr::group_by(`Tissue`, `Gene Symbol`) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * `Cell Count`) / sum(`Cell Count`),
    Total_Cell_Count = sum(`Cell Count`),
    .groups = 'drop'
  )


##as explained in the main manuscript, TOX is selected as a representative T cell exhaustion marker because it is a  master transcriptional regulator 
##that enforces and stabilises the exhausted epigenetic program, effectively “locking” T cells into an exhaustion state.
##we keep all other metabolic genes in there to assess TOX co-expression with metabolic genes across tissues
filtered_data_T_exhaustion <- weighted_expression %>%
  dplyr::filter(`Gene Symbol` == "TOX" | `Gene Symbol` == "LAG3" | `Gene Symbol` == "CTLA4" | `Gene Symbol` == "HAVCR2" | `Gene Symbol` == "PDCD1") %>% 
  dplyr::group_by(`Gene Symbol`) %>%
  dplyr::mutate(
    normalised_Expression =
      Weighted_Expression / max(Weighted_Expression, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()


filtered_data_T_exhaustion_heatmap <- as.data.table(filtered_data_T_exhaustion)

heatmap_matrix_T_exhaustion <- filtered_data_T_exhaustion_heatmap %>%
  dcast(Tissue ~ `Gene Symbol`, value.var = "normalised_Expression") %>%
  tibble::column_to_rownames(var = "Tissue") %>%
  as.matrix()

pheatmap(
  heatmap_matrix_T_exhaustion,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "blue"))(100),
  main = "Gene expression by 'T cell' in normal tissue"
)



##summarise focus exhaustion gene (not truly necessary step since we have one cell type and one gene, but for consistency of formatting):
mean_expression_T_exhaustion <- filtered_data_T_exhaustion %>%
  dplyr::filter(`Gene Symbol` == "TOX") %>% 
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)   #should be same per 
  )




##if including complete set of genes:
filtered_data_T_exhaustion_all <- weighted_expression %>%
  dplyr::group_by(`Gene Symbol`) %>%
  dplyr::mutate(
    normalised_Expression =
      Weighted_Expression / max(Weighted_Expression, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()


filtered_data_T_exhaustion_heatmap <- as.data.table(filtered_data_T_exhaustion_all)

heatmap_matrix_T_exhaustion <- filtered_data_T_exhaustion_heatmap %>%
  dcast(Tissue ~ `Gene Symbol`, value.var = "normalised_Expression") %>%
  tibble::column_to_rownames(var = "Tissue") %>%
  as.matrix()

pheatmap(
  heatmap_matrix_T_exhaustion,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "blue"))(100),
  main = "Gene expression by 'T cell' in normal tissue"
)



##summarise focus exhaustion gene (not truly necessary step since we have one cell type and one gene, but for consistency of formatting):
mean_expression_T_exhaustion <- filtered_data_T_exhaustion %>%
  dplyr::filter(`Gene Symbol` == "TOX") %>% 
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)   
  )

mean_expression_T_exhaustion_all <- filtered_data_T_exhaustion_all %>%
  dplyr::filter(`Gene Symbol` == "TOX" | `Gene Symbol` == "LAG3" | `Gene Symbol` == "CTLA4" | `Gene Symbol` == "HAVCR2" | `Gene Symbol` == "PDCD1") %>% 
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)   
  )


##Next, assess metabolic intensities across tissue intrinsic/stromal cells:
tissue_intrinsic_data <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b1, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")


#uncomment if wanting to group urinary bladder and bladder organ together:
tissue_intrinsic_data <- tissue_intrinsic_data %>%
  dplyr::mutate(
    Tissue = case_when(
      Tissue %in% c("bladder organ", "urinary bladder") ~ "bladder",
      TRUE                                              ~ Tissue
    )
  )


weighted_expression_intrinic <- tissue_intrinsic_data %>%
  dplyr::filter(!is.na(Expression), !is.na(`Cell Count`)) %>%
  dplyr::group_by(`Tissue`, `Gene Symbol`) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * `Cell Count`) / sum(`Cell Count`),
    Total_Cell_Count = sum(`Cell Count`),
    .groups = 'drop'
  )


##take out exhaustion genes
weighted_expression_intrinic_1 <- subset(weighted_expression_intrinic, `Gene Symbol` != "TOX" & `Gene Symbol` != "LAG3" & `Gene Symbol` != "CTLA4" & `Gene Symbol` != "HAVCR2" & `Gene Symbol` != "PDCD1")

filtered_data_intrinic_1 <- weighted_expression_intrinic_1 %>%
  dplyr::group_by(`Gene Symbol`) %>%
  dplyr::mutate(
    normalised_Expression =
      Weighted_Expression / max(Weighted_Expression, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

filtered_data_intrinic_1_heatmap <- as.data.table(filtered_data_intrinic_1)


heatmap_matrix_intrinic_1 <- filtered_data_intrinic_1_heatmap %>%
  dcast(Tissue ~ `Gene Symbol`, value.var = "normalised_Expression") %>%
  tibble::column_to_rownames(var = "Tissue") %>%
  as.matrix()

pheatmap(
  heatmap_matrix_intrinic_1,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "blue"))(100),
  main = "Gene expression by fibroblasts + endothelial cells in normal tissue"
)

filtered_data_intrinic_1$Class <- ifelse(filtered_data_intrinic_1$`Gene Symbol` != "SLC2A1", "pH_homeostasis", "SLC2A1")

##summarise in across genes per gene panel:
mean_expression_intrinsic_1 <- filtered_data_intrinic_1 %>%
  dplyr::group_by(Tissue, Class) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)    
  )


mean_expression_intrinsic_2 <- subset(filtered_data_intrinic_1, `Gene Symbol` == "SLC2A1") %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean SLC2A1",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)   
  )








#########
CELLxGENE_gene_expression_020626_solid_tissues_b <- CELLxGENE_gene_expression_020626_solid_tissues_a %>%
  dplyr::filter(Tissue %in% unique(at_least_2_pubs$Tissue))

CELLxGENE_gene_expression_020626_solid_tissues_b1 <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b, `Number of Cells Expressing Genes` >= 10)
min(CELLxGENE_gene_expression_020626_solid_tissues_b1$`Number of Cells Expressing Genes`)

CELLxGENE_gene_expression_020626_solid_tissues_b1

tissue_intrinsic_data <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b1, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")

tissue_intrinsic_data <- tissue_intrinsic_data %>%
  dplyr::mutate(
    Tissue = case_when(
      Tissue %in% c("bladder organ", "urinary bladder") ~ "bladder",
      TRUE                                              ~ Tissue
    )
  )

weighted_expression_intrinic <- tissue_intrinsic_data %>%
  dplyr::filter(!is.na(Expression), !is.na(`Cell Count`)) %>%
  dplyr::group_by(`Tissue`, `Gene Symbol`) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * `Cell Count`) / sum(`Cell Count`),
    Total_Cell_Count = sum(`Cell Count`),
    .groups = 'drop'
  )


weighted_expression_intrinic_1 <- subset(weighted_expression_intrinic, `Gene Symbol` != "TOX" & `Gene Symbol` != "LAG3" & `Gene Symbol` != "CTLA4" & `Gene Symbol` != "HAVCR2" & `Gene Symbol` != "PDCD1")

filtered_data_intrinic_1 <- weighted_expression_intrinic_1 %>%
  dplyr::group_by(`Gene Symbol`) %>%
  dplyr::mutate(
    normalised_Expression =
      Weighted_Expression / max(Weighted_Expression, na.rm = TRUE)
  ) 


filtered_data_intrinic_1$Class <- ifelse(filtered_data_intrinic_1$`Gene Symbol` != "SLC2A1", "pH_homeostasis", "SLC2A1")
mean_expression_intrinsic_1 <- filtered_data_intrinic_1 %>%
  dplyr::group_by(Tissue, Class) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)    
  )

##Next, we prep the correlation matrix:
##combine intrinsic panels into wide form
intrinsic_wide <- subset(mean_expression_intrinsic_1, Tissue != "adrenal gland") %>% ##outlier excluded
  dplyr::select(Tissue, Class, Weighted_Expression) %>%
  tidyr::pivot_wider(names_from = Class, values_from = Weighted_Expression)


###add the T-cell exhaustion panel (so beyond just TOX) as a column
all_panels <- intrinsic_wide %>%
  dplyr::inner_join(
    mean_expression_T_exhaustion_all %>%
      dplyr::select(Tissue, T_cell_exhaustion = Weighted_Expression),
    by = "Tissue"
  )

##remove Tissue column & compute a full Spearman correlation matrix:
expr_mat <- all_panels %>%
  dplyr::ungroup() %>%           
  dplyr::select(-Tissue) %>%
  as.data.frame()         

cor_mat <- cor(expr_mat, method = "spearman", use = "pairwise.complete.obs")


get_p <- function(x) {
  n <- nrow(expr_mat)
  t_val <- x * sqrt((n - 2) / (1 - x^2))
  2 * pt(-abs(t_val), df = n - 2)
}
p_mat <- matrix(get_p(cor_mat),
                nrow = ncol(expr_mat),
                dimnames = list(colnames(expr_mat), colnames(expr_mat)))

sig_mat <- ifelse(p_mat < 0.001, "***",
                  ifelse(p_mat < 0.01,  "**",
                         ifelse(p_mat < 0.05,  "*", "")))

diag(cor_mat) <- NA 
diag(sig_mat) <- ""


pheatmap(
  cor_mat,
  color = colorRampPalette(c("white","royalblue3"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = sig_mat,
  number_color = "black",
  main = "Spearman correlations between mean weighted expressions, 
all solid tissues (N = 29) "
)


my_colours <- distinctColorPalette(28) 

tissues <- c(
  "adipose tissue","adrenal gland","bladder","brain","breast",
  "colon","endocrine gland","esophagus","exocrine gland","eye",
  "fallopian tube","heart","intestine","kidney","large intestine",
  "liver","lung","musculature","ovary","pancreas","prostate gland",
  "skin of body","small intestine","spleen","stomach","ureter",
  "uterus"
)

my_colours <- setNames(distinctColorPalette(length(tissues)), tissues)

base_cols <- brewer.pal(12, "Paired")           
my_colours <- setNames(colorRampPalette(base_cols)(length(tissues)), tissues)



##investigate correlations between individual panels:
ggplot(all_panels,
       aes(x = pH_homeostasis,
           y = T_cell_exhaustion,
           color = Tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = my_colours) +  
  labs(
    x = "pH_homeostasis",
    y = "T_cell_exhaustion",
    title = "Per-tissue correlation",
    color = "Tissue"
  ) +
  theme_classic(base_size = 14) +            
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.background = element_blank()
  )



ggplot(all_panels,
       aes(x = SLC2A1,
           y = T_cell_exhaustion,
           color = Tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = my_colours) +  
  labs(
    x = "SLC2A1",
    y = "T_cell_exhaustion",
    title = "Per-tissue correlation",
    color = "Tissue"
  ) +
  theme_classic(base_size = 14) +            
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.background = element_blank()
  )




####SLC2A1 as proxy for metabolism!
##focus genes only:
merged <- mean_expression_intrinsic_2 %>%
  rename(
    Weighted_Expression_intrinsic = Weighted_Expression,
    Total_Cell_Count_intrinsic    = Total_Cell_Count
  ) %>%
  dplyr::inner_join(
    mean_expression_T_exhaustion_all %>% ##_all for complete exhaustion panel, 
      rename(
        Weighted_Expression_T_exhaustion = Weighted_Expression,
        Total_Cell_Count_T_exhaustion    = Total_Cell_Count
      ),
    by = "Tissue"
  )

merged <- merged %>%
  rename(
    Gene_Symbol_intrinsic  = `Gene Symbol.x`,
    Gene_Symbol_T_exhaustion = `Gene Symbol.y`
  )



distinct_colours <- c(
  "#1f77b4", "#ff7f0e", "green2", "#d62728", "#9467bd", "#8c564b", "blue", "grey40",
  "#bcbd82", "#17becf", "darkgreen", "coral3", "black", "#66c2a5", "#fc8d62", "#8da0cb",
  "yellow", "#a6d854", "#ffc32f", "#e5c494", "grey", "#ff69b4", "cyan", "royalblue",
  "purple", "#dda0dd", "orange", "firebrick", "seagreen", "gold", "turquoise", "magenta",
  "darkorchid", "brown", "khaki4", "deeppink", "lightskyblue", "chartreuse3", "orchid", "dodgerblue",
  "sienna", "steelblue", "violetred", "cadetblue", "olivedrab", "plum3", "slateblue", "tomato",
  "darkcyan", "chocolate3", "maroon3", "palegreen3", "skyblue3", "navy", "darkorange2", "darkred",
  "lightpink"
)

##take out outlier adrenal gland: 
merged1 <- subset(merged, Tissue != "adrenal gland")

ggplot(merged1,
       aes(x = Weighted_Expression_intrinsic,
           y = Weighted_Expression_T_exhaustion,
           color = Tissue)) +
  geom_smooth(
    method  = "lm",
    color   = "grey40",
    linetype = "dashed",
    se      = FALSE
  ) +
  geom_point(aes(size = log2(Total_Cell_Count_T_exhaustion / Total_Cell_Count_intrinsic)), alpha = 0.85) +
  scale_size_continuous(name = "Cell-count ratio") +
  scale_color_manual(values = distinct_colours) +  
  labs(
    x = "SLC2A1 expression in fibroblasts + endothelial cells",
    y = "T cell exhaustion panel expression in T cells",
    title = " ",
    color = "Solid tissue"
  ) +
  theme_classic(base_size = 14) +            
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.background = element_blank()
  )



summary(lm(merged1$Weighted_Expression_T_exhaustion ~ merged1$Weighted_Expression_intrinsic))
cor.test(merged1$Weighted_Expression_intrinsic, merged1$Weighted_Expression_T_exhaustion, method = "spearman")




###14 focus tissues only:
CELLxGENE_gene_expression_020626_solid_tissues_b1 <- CELLxGENE_gene_expression_020626_solid_tissues_b1 %>%
  dplyr::mutate(
    Tissue = case_when(
      Tissue %in% c("bladder organ", "urinary bladder") ~ "bladder",
      TRUE                                              ~ Tissue
    )
  )


CELLxGENE_gene_expression_020626_solid_tissues_b1 <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b1, Tissue == "colon" | Tissue == "brain" | Tissue == "breast"
                                                            | Tissue == "esophagus" | Tissue == "eye" | Tissue == "fallopian tube" | Tissue == "kidney" | Tissue == "liver"
                                                            | Tissue == "lung" | Tissue == "pancreas" | Tissue == "prostate gland" | Tissue == "skin of body" | Tissue == "stomach" | Tissue == "bladder")


t_cell_data <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b1, `Cell Type` == "T cell")


weighted_expression <- t_cell_data %>%
  dplyr::filter(!is.na(Expression), !is.na(`Cell Count`)) %>%
  dplyr::group_by(`Tissue`, `Gene Symbol`) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * `Cell Count`) / sum(`Cell Count`),
    Total_Cell_Count = sum(`Cell Count`),
    .groups = 'drop'
  )

filtered_data_T_exhaustion <- weighted_expression %>%
  dplyr::filter(`Gene Symbol` %in% c("TOX")) %>%
  dplyr::group_by(`Gene Symbol`) %>%
  dplyr::mutate(
    normalised_Expression =
      Weighted_Expression / max(Weighted_Expression, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

filtered_data_T_exhaustion_all <- weighted_expression %>%
  dplyr::filter(`Gene Symbol` %in% c("CTLA4","LAG3","PDCD1","HAVCR2","TOX")) %>%
  dplyr::group_by(`Gene Symbol`) %>%
  dplyr::mutate(
    normalised_Expression =
      Weighted_Expression / max(Weighted_Expression, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()


labs4 <- c("Low", "Intermediate",
           "High", "Very high")


##summarise across genes in the panel:
mean_expression_T_exhaustion_all <- filtered_data_T_exhaustion_all %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)   
  )

mean_expression_T_exhaustion_all <- mean_expression_T_exhaustion_all %>%
  dplyr::mutate(
    Tissue_bin = dplyr::ntile(Weighted_Expression, 4),
    Tissue_bin = factor(labs4[Tissue_bin], levels = labs4)
  )


##TOX only
mean_expression_T_exhaustion <- filtered_data_T_exhaustion %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "TOX",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)   
  )

mean_expression_T_exhaustion <- mean_expression_T_exhaustion %>%
  dplyr::mutate(
    Tissue_bin = dplyr::ntile(Weighted_Expression, 4),
    Tissue_bin = factor(labs4[Tissue_bin], levels = labs4)
  )



##Next, assess metabolic intensities across tissue intrinsic/stromal cells:
tissue_intrinsic_data <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b1, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")

weighted_expression_intrinic <- tissue_intrinsic_data %>%
  dplyr::filter(!is.na(Expression), !is.na(`Cell Count`)) %>%
  dplyr::group_by(`Tissue`, `Gene Symbol`) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * `Cell Count`) / sum(`Cell Count`),
    Total_Cell_Count = sum(`Cell Count`),
    .groups = 'drop'
  )


weighted_expression_intrinic_1 <- subset(weighted_expression_intrinic, `Gene Symbol` != "TOX" & `Gene Symbol` != "LAG3" & `Gene Symbol` != "CTLA4" & `Gene Symbol` != "HAVCR2" & `Gene Symbol` != "PDCD1")

filtered_data_intrinic_1 <- weighted_expression_intrinic_1 %>%
  dplyr::group_by(`Gene Symbol`) %>%
  dplyr::mutate(
    normalised_Expression =
      Weighted_Expression / max(Weighted_Expression, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

filtered_data_intrinic_1_heatmap <- as.data.table(filtered_data_intrinic_1)

filtered_data_intrinic_1_heatmap <- subset(filtered_data_intrinic_1_heatmap, `Gene Symbol` != "SDHA")

heatmap_matrix_intrinic_1 <- filtered_data_intrinic_1_heatmap %>%
  dcast(Tissue ~ `Gene Symbol`, value.var = "normalised_Expression") %>%
  tibble::column_to_rownames(var = "Tissue") %>%
  as.matrix()

##note split here in dendogram follows the overall classifications listed in Table 1
pheatmap(
  heatmap_matrix_intrinic_1,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "blue"))(100),
  main = "Gene expression by fibroblasts + endothelial cells in normal tissue"
)


filtered_data_intrinic_1$Class <- ifelse(filtered_data_intrinic_1$`Gene Symbol` != "SLC2A1", "pH_homeostasis", "SLC2A1")

##summarise in across genes per gene panel (Class):
mean_expression_intrinsic_1 <- filtered_data_intrinic_1 %>%
  dplyr::group_by(Tissue, Class) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)   
  )



str(filtered_data_intrinic_1)


mean_expression_intrinsic_metabolism <- filtered_data_intrinic_1 %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)   
  )


labs4 <- c("Low", "Intermediate",
           "High", "Very high")

mean_expression_intrinsic_2a <- subset(filtered_data_intrinic_1, `Gene Symbol` == "SLC2A1") %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "SLC2A1",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)    
  )


mean_expression_intrinsic_2a <- mean_expression_intrinsic_2a %>%
  dplyr::mutate(
    Tissue_bin = dplyr::ntile(Weighted_Expression, 4),
    Tissue_bin = factor(labs4[Tissue_bin], levels = labs4)
  )




mean_expression_intrinsic_3 <- subset(filtered_data_intrinic_1, Class == "pH_homeostasis") %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean pH homeostasis genes",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)  
  )


mean_expression_intrinsic_3 <- mean_expression_intrinsic_3 %>%
  dplyr::mutate(
    Tissue_bin = dplyr::ntile(Weighted_Expression, 4),
    Tissue_bin = factor(labs4[Tissue_bin], levels = labs4)
  )



##prep the correlation matrix:
#combine intrinsic panels into wide form
intrinsic_wide <- mean_expression_intrinsic_1 %>%
  dplyr::select(Tissue, Class, Weighted_Expression) %>%
  tidyr::pivot_wider(names_from = Class, values_from = Weighted_Expression)

###add the T-cell exhaustion panel as a column
all_panels <- intrinsic_wide %>%
  dplyr::inner_join(
    mean_expression_T_exhaustion_all %>% ##all for complete panel
      dplyr::select(Tissue, T_cell_exhaustion = Weighted_Expression),
    by = "Tissue"
  )

##remove Tissue column and compute a full Spearman correlation matrix
expr_mat <- all_panels %>%
  dplyr::ungroup() %>%           
  dplyr::select(-Tissue) %>%
  as.data.frame()         

cor_mat <- cor(expr_mat, method = "spearman", use = "pairwise.complete.obs")


get_p <- function(x) {
  n <- nrow(expr_mat)
  t_val <- x * sqrt((n - 2) / (1 - x^2))
  2 * pt(-abs(t_val), df = n - 2)
}
p_mat <- matrix(get_p(cor_mat),
                nrow = ncol(expr_mat),
                dimnames = list(colnames(expr_mat), colnames(expr_mat)))

sig_mat <- ifelse(p_mat < 0.001, "***",
                  ifelse(p_mat < 0.01,  "**",
                         ifelse(p_mat < 0.05,  "*", "")))

diag(cor_mat) <- NA 
diag(sig_mat) <- ""


pheatmap(
  cor_mat,
  color = colorRampPalette(c("white","royalblue3"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = sig_mat,
  number_color = "black",
  main = "Spearman correlations between mean weighted expressions, 
all focus (Table 1) tissues (N = 14)"
)


my_colours <- distinctColorPalette(28) 

tissues <- c(
  "brain","breast",
  "colon","esophagus","eye",
  "fallopian tube", "kidney",
  "liver","lung","pancreas","prostate gland",
  "skin of body","stomach",
  "bladder"
)

my_colours <- setNames(distinctColorPalette(length(tissues)), tissues)

base_cols <- brewer.pal(12, "Paired")           
my_colours <- setNames(colorRampPalette(base_cols)(length(tissues)), tissues)


ggplot(all_panels,
       aes(x = pH_homeostasis,
           y = T_cell_exhaustion,
           color = Tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = my_colours) +  
  labs(
    x = "pH_homeostasis",
    y = "T_cell_exhaustion",
    title = "Cross-tissue correlation",
    color = "Tissue"
  ) +
  theme_classic(base_size = 14) +            
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.background = element_blank()
  )



####SLC2A1 as proxy for metabolism!
##focus genes only:
merged <- mean_expression_intrinsic_2a %>%
  rename(
    Weighted_Expression_intrinsic = Weighted_Expression,
    Total_Cell_Count_intrinsic    = Total_Cell_Count
  ) %>%
  dplyr::inner_join(
    mean_expression_T_exhaustion_all %>%
      rename(
        Weighted_Expression_T_exhaustion = Weighted_Expression,
        Total_Cell_Count_T_exhaustion    = Total_Cell_Count
      ),
    by = "Tissue"
  )

merged <- merged %>%
  rename(
    Gene_Symbol_intrinsic  = `Gene Symbol.x`,
    Gene_Symbol_T_exhaustion = `Gene Symbol.y`
  )


ggplot(merged,
       aes(x = Weighted_Expression_intrinsic,
           y = Weighted_Expression_T_exhaustion,
           color = Tissue)) +
  geom_point(aes(size = log2(Total_Cell_Count_T_exhaustion / Total_Cell_Count_intrinsic)), alpha = 0.95) +
  geom_smooth(
    method  = "lm",
    color   = "grey20",
    linetype = "dashed",
    se      = FALSE
  ) +
  scale_size_continuous(name = "Cell-count ratio") +
  scale_color_manual(values = my_colours) +  
  labs(
    x = "SLC2A1 expression by fibroblasts + endothelial cells",
    y = "T cell exhaustion panel expression by T cells",
    title = "Per-tissue correlation",
    color = "Solid tissue"
  ) +
  theme_classic(base_size = 14) +            
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.background = element_blank()
  )

summary(lm(merged$Weighted_Expression_intrinsic ~ merged$Weighted_Expression_T_exhaustion))
cor.test(merged$Weighted_Expression_T_exhaustion, merged$Weighted_Expression_intrinsic, method = "spearman")



##SLC2A1 vs TOX only
merged <- mean_expression_intrinsic_2a %>%
  rename(
    Weighted_Expression_intrinsic = Weighted_Expression,
    Total_Cell_Count_intrinsic    = Total_Cell_Count
  ) %>%
  dplyr::inner_join(
    mean_expression_T_exhaustion %>%
      rename(
        Weighted_Expression_T_exhaustion = Weighted_Expression,
        Total_Cell_Count_T_exhaustion    = Total_Cell_Count
      ),
    by = "Tissue"
  )

merged <- merged %>%
  rename(
    Gene_Symbol_intrinsic  = `Gene Symbol.x`,
    Gene_Symbol_T_exhaustion = `Gene Symbol.y`
  )


ggplot(merged,
       aes(x = Weighted_Expression_intrinsic,
           y = Weighted_Expression_T_exhaustion,
           color = Tissue)) +
  geom_point(aes(size = log2(Total_Cell_Count_T_exhaustion / Total_Cell_Count_intrinsic)), alpha = 0.8) +
  geom_smooth(
    method  = "lm",
    color   = "grey20",
    linetype = "dashed",
    se      = FALSE
  ) +
  scale_size_continuous(name = "Cell-count ratio") +
  scale_color_manual(values = my_colours) +  
  labs(
    x = "SLC2A1 expression by fibroblasts + endothelial cells",
    y = "TOX expression by T cells",
    title = "Per-tissue correlation",
    color = "Solid tissue"
  ) +
  theme_classic(base_size = 14) +            
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.background = element_blank()
  )

summary(lm(merged$Weighted_Expression_intrinsic ~ merged$Weighted_Expression_T_exhaustion))
cor.test(merged$Weighted_Expression_T_exhaustion, merged$Weighted_Expression_intrinsic, method = "spearman")







####pH homeostasis genes as proxies for acidosis management/ response
##focus genes only:
merged_2 <- mean_expression_intrinsic_3 %>%
  rename(
    Weighted_Expression_intrinsic = Weighted_Expression,
    Total_Cell_Count_intrinsic    = Total_Cell_Count
  ) %>%
  dplyr::inner_join(
    mean_expression_T_exhaustion_all %>%
      rename(
        Weighted_Expression_T_exhaustion = Weighted_Expression,
        Total_Cell_Count_T_exhaustion    = Total_Cell_Count
      ),
    by = "Tissue"
  )

merged_2 <- merged_2 %>%
  rename(
    Gene_Symbol_intrinsic  = `Gene Symbol.x`,
    Gene_Symbol_T_exhaustion = `Gene Symbol.y`
  )

##Note that the ratio is dependent on the panels and cell types elected; in the paper (Table 1), we show an approximation based on ratios deduced when 
##taking the count of stromal cells (endothelial + fibroblasts) expressing the panel of pH homeostasis genes versus the count of T cells expressing the larger panel of T cell exhaustion genes (to ensure the ratio 'captures' T cells that may not express TOX)
merged_2$ratio <- merged_2$Total_Cell_Count_intrinsic / merged_2$Total_Cell_Count_T_exhaustion

ggplot(merged_2,
       aes(x = Weighted_Expression_intrinsic,
           y = Weighted_Expression_T_exhaustion,
           color = Tissue)) +
  geom_point(aes(size = log2(Total_Cell_Count_T_exhaustion / Total_Cell_Count_intrinsic)), alpha = 0.9) +
  geom_smooth(
    method  = "lm",
    color   = "grey20",
    linetype = "dashed",
    se      = FALSE
  ) +
  scale_size_continuous(name = "Cell-count ratio") +
  scale_color_manual(values = my_colours) +  
  labs(
    x = "pH homeostasis panel expression by fibroblasts + endothelial cells",
    y = "T cell exhaustion panel expression by T cells",
    title = "Per-tissue correlation",
    color = "Solid tissue"
  ) +
  theme_classic(base_size = 14) +            
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.background = element_blank()
  )

summary(lm(merged_2$Weighted_Expression_T_exhaustion~merged_2$Weighted_Expression_intrinsic))
cor.test(merged_2$Weighted_Expression_T_exhaustion, merged_2$Weighted_Expression_intrinsic, method = "spearman")




