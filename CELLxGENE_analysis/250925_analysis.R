##Single cell RNAseq analysis publicly available CELLxGENE data across all tissues for human, normal
##In this script, we first explore the inter Publication and inter Sex sources of variance, in order to determine 
#whether the default download of CELLxGENE (one aggregated Expression value for each unique combination of Cell Type, Gene Symbol and Tissue)
#can reliably be used to investigate coarse-grained cross-tissue expression patterns for a given subset of genes and cell types.
##We then analyse cross-tissue patterns including solid tissues other than those listed in the manuscript, to explore robustness of the patterns
##Note further that we also include marker panels not discussed in the manuscript to further explore robustness of patterns.
##run script chronologically due to reuse of object names
##From line ~860 onwards, we explore trends exclusively in the 14 focus tissues (and for focus genes listed in Table 1).


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


##T cell exhaustion:
CELLxGENE_gene_expression_092525_T_cell_exhaustion_by_publication <- read_csv("Documents/CELLxGENE_Sept25/T_cell_exhaustion/CELLxGENE_gene_expression_092525.csv", skip = 9)
CELLxGENE_gene_expression_092525_T_cell_exhaustion_by_sex <- read_csv("Documents/CELLxGENE_Sept25/T_cell_exhaustion/CELLxGENE_gene_expression_092525 (1).csv", skip = 9)
CELLxGENE_gene_expression_092525_T_cell_exhaustion <- read_csv("Documents/CELLxGENE_Sept25/T_cell_exhaustion/CELLxGENE_gene_expression_092525 (2).csv", skip = 9)

##Metabolism ones:
CELLxGENE_gene_expression_092525_beta_oxidation_by_publication <- read_csv("Documents/CELLxGENE_Sept25/beta_oxidation/CELLxGENE_gene_expression_092525.csv", skip = 9)
CELLxGENE_gene_expression_092525_beta_oxidation_by_sex <- read_csv("Documents/CELLxGENE_Sept25/beta_oxidation/CELLxGENE_gene_expression_092525 (1).csv", skip = 9)
CELLxGENE_gene_expression_092525_beta_oxidation <- read_csv("Documents/CELLxGENE_Sept25/beta_oxidation/CELLxGENE_gene_expression_092525 (2).csv", skip = 9)

CELLxGENE_gene_expression_092525_glycolysis_and_glucose_uptake_by_publication <- read_csv("Documents/CELLxGENE_Sept25/glycolysis_and_glucose_uptake/CELLxGENE_gene_expression_092525.csv", skip = 9)
CELLxGENE_gene_expression_092525_glycolysis_and_glucose_uptake_by_sex <- read_csv("Documents/CELLxGENE_Sept25/glycolysis_and_glucose_uptake/CELLxGENE_gene_expression_092525 (1).csv", skip = 9)
CELLxGENE_gene_expression_092525_glycolysis_and_glucose_uptake <- read_csv("Documents/CELLxGENE_Sept25/glycolysis_and_glucose_uptake/CELLxGENE_gene_expression_092525 (2).csv", skip = 9)

CELLxGENE_gene_expression_092525_lactate_oxidation_by_publication <- read_csv("Documents/CELLxGENE_Sept25/lactate_oxidation/CELLxGENE_gene_expression_092525.csv", skip = 9)
CELLxGENE_gene_expression_092525_lactate_oxidation_by_sex <- read_csv("Documents/CELLxGENE_Sept25/lactate_oxidation/CELLxGENE_gene_expression_092525 (1).csv", skip = 9)
CELLxGENE_gene_expression_092525_lactate_oxidation <- read_csv("Documents/CELLxGENE_Sept25/lactate_oxidation/CELLxGENE_gene_expression_092525 (2).csv", skip = 9)

CELLxGENE_gene_expression_092525_OXPHOS_by_publication <- read_csv("Documents/CELLxGENE_Sept25/OXPHOS/CELLxGENE_gene_expression_092525.csv", skip = 9)
CELLxGENE_gene_expression_092525_OXPHOS_by_sex <- read_csv("Documents/CELLxGENE_Sept25/OXPHOS/CELLxGENE_gene_expression_092525 (1).csv", skip = 9)
CELLxGENE_gene_expression_092525_OXPHOS <- read_csv("Documents/CELLxGENE_Sept25/OXPHOS/CELLxGENE_gene_expression_092525 (2).csv", skip = 9)

CELLxGENE_gene_expression_092525_TCA_cycle_by_publication <- read_csv("Documents/CELLxGENE_Sept25/TCA_cycle/CELLxGENE_gene_expression_092525.csv", skip = 9)
CELLxGENE_gene_expression_092525_TCA_cycle_by_sex <- read_csv("Documents/CELLxGENE_Sept25/TCA_cycle/CELLxGENE_gene_expression_092525 (1).csv", skip = 9)
CELLxGENE_gene_expression_092525_TCA_cycle <- read_csv("Documents/CELLxGENE_Sept25/TCA_cycle/CELLxGENE_gene_expression_092525 (2).csv", skip = 9)

##pH homeostasis
CELLxGENE_gene_expression_092525_pH_homeostasis_by_publication <- read_csv("Documents/CELLxGENE_Sept25/pH_homeostasis/CELLxGENE_gene_expression_092525.csv", skip = 9)
CELLxGENE_gene_expression_092525_pH_homeostasis_by_sex <- read_csv("Documents/CELLxGENE_Sept25/pH_homeostasis/CELLxGENE_gene_expression_092525 (1).csv", skip = 9)
CELLxGENE_gene_expression_092525_pH_homeostasis <- read_csv("Documents/CELLxGENE_Sept25/pH_homeostasis/CELLxGENE_gene_expression_092525 (2).csv", skip = 9)

##ECM density
CELLxGENE_gene_expression_092525_ECM_density_by_publication <- read_csv("Documents/CELLxGENE_Sept25/ECM_density/CELLxGENE_gene_expression_092525.csv", skip = 9)
CELLxGENE_gene_expression_092525_ECM_density_by_sex <- read_csv("Documents/CELLxGENE_Sept25/ECM_density/CELLxGENE_gene_expression_092525 (1).csv", skip = 9)
CELLxGENE_gene_expression_092525_ECM_density <- read_csv("Documents/CELLxGENE_Sept25/ECM_density/CELLxGENE_gene_expression_092525 (2).csv", skip = 9)



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
  CELLxGENE_gene_expression_092525_T_cell_exhaustion_by_publication,
  CELLxGENE_gene_expression_092525_beta_oxidation_by_publication,
  CELLxGENE_gene_expression_092525_glycolysis_and_glucose_uptake_by_publication,
  CELLxGENE_gene_expression_092525_lactate_oxidation_by_publication,
  CELLxGENE_gene_expression_092525_OXPHOS_by_publication,
  CELLxGENE_gene_expression_092525_TCA_cycle_by_publication,
  CELLxGENE_gene_expression_092525_pH_homeostasis_by_publication,
  CELLxGENE_gene_expression_092525_ECM_density_by_publication
)

names(datasets) <- c(
  "T_exhaustion_by_pub",
  "beta_by_pub",
  "glycolysis_by_pub",
  "lactate_by_pub",
  "OXPHOS_by_pub",
  "TCA_by_pub",
  "pH_by_pub",
  "ECM_by_pub"
)

cleaned_datasets <- lapply(datasets, clean_cellxgene)
list2env(cleaned_datasets, envir = .GlobalEnv)




##investigate number of publications per tissue
T_exhaustion_by_pub$Class <- "T_cell_exhaustion"
beta_by_pub$Class <- "beta_oxidation"
glycolysis_by_pub$Class <- "glycolysis_and_glucose_uptake"
lactate_by_pub$Class <- "lactate_oxidation"
OXPHOS_by_pub$Class <- "OXPHOS"
TCA_by_pub$Class <- "TCA_cycle"
pH_by_pub$Class <- "pH_homeostasis"
ECM_by_pub$Class <- "ECM_density"


CELLxGENE_gene_expression_092525_by_publication <- rbind(T_exhaustion_by_pub, beta_by_pub,
                                                         glycolysis_by_pub, lactate_by_pub,
                                                         OXPHOS_by_pub, TCA_by_pub,
                                                         pH_by_pub, ECM_by_pub)



##We exclude: non-solid tissues (circulating or fluid-derived), embryonic or developmental tissues, or tissue subcompartments rather than whole solid organs
CELLxGENE_gene_expression_092525_by_publication_solid_tissues <- subset(CELLxGENE_gene_expression_092525_by_publication, Tissue != "blood" 
                                                             & Tissue != "bone marrow" & Tissue != "lymph node" 
                                                             & Tissue != "milk" & Tissue != "pleural fluid" & Tissue != "saliva" & Tissue != "yolk sac" & Tissue != "embryo" 
                                                             & Tissue != "placenta" & Tissue != "vasculature" & Tissue != "mucosa" & Tissue != "lamina propria" & Tissue != "omentum")


CELLxGENE_gene_expression_092525_by_publication_solid_tissues_a <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues, Publication != "aggregated")
CELLxGENE_gene_expression_092525_by_publication_solid_tissues_a <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_a, `Number of Cells Expressing Genes` >=10)
min(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_a$`Number of Cells Expressing Genes`)

#View(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_a)

##if not filtering by cell type:
check_publications <- CELLxGENE_gene_expression_092525_by_publication_solid_tissues_a %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(n_publications = n_distinct(Publication), .groups = "drop")

##only keeping publications that contain any of the cell types of interest:
check_publications <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_a, `Cell Type` == "fibroblast" | `Cell Type` == "T cell" | `Cell Type` == "endothelial cell") %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(n_publications = n_distinct(Publication), .groups = "drop")


#need to an exception for urinary bladder as we combine this with "bladder organ" as a representative for "bladder", to arrive at a total of 3 datasets:
at_least_2_pubs <- subset(
  check_publications,
  n_publications > 1 | Tissue == "urinary bladder"
)

CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b <- CELLxGENE_gene_expression_092525_by_publication_solid_tissues_a %>%
  dplyr::filter(Tissue %in% unique(at_least_2_pubs$Tissue))


unique(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b$Tissue)
unique(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b$`Cell Type`)



###Prepare to investigate inter-publication variability in Expression levels:
##first, to arrive at 3 publications representing the tissue of origin for bladder cancer, we here regard "urinary bladder" and "bladder organ" as "bladder":
CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b1 <- 
  CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b %>%
  dplyr::mutate(Tissue = case_when(
    Tissue %in% c("urinary bladder", "bladder organ") ~ "bladder",
    TRUE ~ Tissue
  ))


##Publication-level variation explored per tissue*gene*cell type combo:
CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised <- 
  CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b1 %>%
  dplyr::group_by(`Gene Symbol`, Tissue, `Cell Type`) %>%
  dplyr::summarise(
    median_expr = median(Expression, na.rm = TRUE),
    mad_expr = mad(Expression, na.rm = TRUE),
    n_studies = sum(!is.na(Expression))
  )


CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised_fibroblast <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised, `Cell Type` == "fibroblast")
CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised_endothelial_cell <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised, `Cell Type` == "endothelial cell")
CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised_T_cell <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised, `Cell Type` == "T cell")

CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised_stromal <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")
unique(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised_stromal$`Cell Type`)


###if only keeping focus tissues: 
CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised_stromal_focus <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised_stromal, Tissue == "colon" | Tissue == "brain" | Tissue == "breast"
                                                                      | Tissue == "brain" | Tissue == "esophagus" | Tissue == "eye" | Tissue == "fallopian tube" | Tissue == "kidney" | Tissue == "liver"
                                                                      | Tissue == "lung" | Tissue == "pancreas" | Tissue == "prostate gland" | Tissue == "skin of body" | Tissue == "stomach" | Tissue == "bladder")


CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised_stromal_focus_no_T_exhaustion_genes <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised_stromal_focus, `Gene Symbol` != "TOX" 
                                                                                                                            & `Gene Symbol` != "LAG3" & `Gene Symbol` != "PDCD1" & `Gene Symbol` != "TOX" & `Gene Symbol` != "HAVCR2")

max_median_expr <- max(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised_stromal_focus_no_T_exhaustion_genes$median_expr, na.rm = TRUE)
max_mad_expr <- max(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b_summarised_stromal_focus_no_T_exhaustion_genes$mad_expr, na.rm = TRUE)


####investigate role of inter-publication variance for stromal cells across all solid tissues, excluding T cell exhaustion genes
CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b1_stromal <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b1, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")
CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b1_stromal <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b1_stromal, `Gene Symbol` != "TOX" 
                                                                                      & `Gene Symbol` != "LAG3" & `Gene Symbol` != "PDCD1" & `Gene Symbol` != "CTLA4" & `Gene Symbol` != "HAVCR2")


##Publication-to-publication differences (batch effects, methods, biases) are real but smaller than biological differences across genes:
##Because the dominant source of signal is still the Gene Symbol itself rather than Publication, we can go ahead with reliably interpreting and comparing expressions across genes
##Note that Tissue is omitted from the model because it is strongly confounded with Publication (many publications sample only one tissue), so a separate Tissue random effect is not reliably estimable and its variance would mostly be absorbed into Publication
##Publication and Tissue are thus partially confounded, so the Publication variance likely overestimates purely technical/batch effect noise
fit_stromal_cells_all_solid_tissues <- lmer(Expression ~ 1 + (1 | Publication) + (1 | `Gene Symbol`),
                                          data = CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b1_stromal)
                                          
VarCorr(fit_stromal_cells_all_solid_tissues)




##now the same but for T cells:
####investigate role of inter-publication variance for T cells across all solid tissues, excluding T cell exhaustion genes
CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b1_T <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b1, `Cell Type` == "T cell")

#take out genes not reflective of T cell biology (i.e., stromal specific genes)
CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b1_T <- subset(CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b1_T, Class != "ECM_density")

fit_T_cells_all_solid_tissues <- lmer(Expression ~ 1 + (1 | Publication) + (1 | `Gene Symbol`),
                                            data = CELLxGENE_gene_expression_092525_by_publication_solid_tissues_b1_T)


##so: For stromal cells, biological variation is dominated by gene identity, Publication noise is in similar order of magnitude but still smaller: so aggregating across publications is assumed safe
##Re. T cells, Publication absorbs both technical and tissue-level biological variation; The signal is a mix of real tissue biology and study-specific noise; but we assume this does not invalidate tissue comparisons
##especially since our focus is on comparing stromal metabolic profiles between tissues.
VarCorr(fit_T_cells_all_solid_tissues)






###Same approach, now by sex instead of publication:
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
  CELLxGENE_gene_expression_092525_T_cell_exhaustion_by_sex,
  CELLxGENE_gene_expression_092525_beta_oxidation_by_sex,
  CELLxGENE_gene_expression_092525_glycolysis_and_glucose_uptake_by_sex,
  CELLxGENE_gene_expression_092525_lactate_oxidation_by_sex,
  CELLxGENE_gene_expression_092525_OXPHOS_by_sex,
  CELLxGENE_gene_expression_092525_TCA_cycle_by_sex,
  CELLxGENE_gene_expression_092525_pH_homeostasis_by_sex,
  CELLxGENE_gene_expression_092525_ECM_density_by_sex
)

names(datasets) <- c(
  "T_exhaustion_by_sex",
  "beta_by_sex",
  "glycolysis_by_sex",
  "lactate_by_sex",
  "OXPHOS_by_sex",
  "TCA_by_sex",
  "pH_by_sex",
  "ECM_by_sex"
)

cleaned_datasets_sex <- lapply(datasets, clean_cellxgene_sex)
list2env(cleaned_datasets_sex, envir = .GlobalEnv)


T_exhaustion_by_sex$Class <- "T_cell_exhaustion"
beta_by_sex$Class <- "beta_oxidation"
glycolysis_by_sex$Class <- "glycolysis_and_glucose_uptake"
lactate_by_sex$Class <- "lactate_oxidation"
OXPHOS_by_sex$Class <- "OXPHOS"
TCA_by_sex$Class <- "TCA_cycle"
pH_by_sex$Class <- "pH_homeostasis"
ECM_by_sex$Class <- "ECM_density"


CELLxGENE_gene_expression_092525_by_sex <- rbind(T_exhaustion_by_sex, beta_by_sex,
                                                         glycolysis_by_sex, lactate_by_sex,
                                                         OXPHOS_by_sex, TCA_by_sex,
                                                         pH_by_sex, ECM_by_sex)


##take out non-solid/developmental/non-specific tissues like before:
CELLxGENE_gene_expression_092525_by_sex_solid_tissues <- subset(CELLxGENE_gene_expression_092525_by_sex, Tissue != "blood"
                                                                        & Tissue != "bone marrow" & Tissue != "lymph node"
                                                                        & Tissue != "milk" & Tissue != "pleural fluid" & Tissue != "saliva" & Tissue != "yolk sac" & Tissue != "embryo"
                                                                        & Tissue != "placenta" & Tissue != "vasculature" & Tissue != "mucosa" & Tissue != "lamina propria" & Tissue != "omentum")


CELLxGENE_gene_expression_092525_by_sex_solid_tissues_a <- subset(CELLxGENE_gene_expression_092525_by_sex_solid_tissues, Sex != "aggregated")
CELLxGENE_gene_expression_092525_by_sex_solid_tissues_a <- subset(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_a, `Number of Cells Expressing Genes` >= 10)
min(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_a$`Number of Cells Expressing Genes`)

#View(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_a)


CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b <- CELLxGENE_gene_expression_092525_by_sex_solid_tissues_a %>%
  dplyr::filter(Tissue %in% unique(at_least_2_pubs$Tissue))

CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1 <- 
  CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b %>%
  dplyr::mutate(Tissue = case_when(
    Tissue %in% c("urinary bladder", "bladder organ") ~ "bladder",
    TRUE ~ Tissue
  ))

unique(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1$Tissue)
unique(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1$`Cell Type`)


CELLxGENE_gene_expression_092525_by_sex_solid_tissues_c_endothelial_cell <- subset(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1, 
                                                                                   (Class %in% c("glycolysis_and_glucose_uptake",
                                                                                                 "beta_oxidation",
                                                                                                 "OXPHOS",
                                                                                                 "TCA_cycle",
                                                                                                 "lactate_oxidation")) &
                                                                                     `Cell Type` == "endothelial cell")

CELLxGENE_gene_expression_092525_by_sex_solid_tissues_c_fibroblast <- subset(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1, 
                                                                             (Class %in% c("glycolysis_and_glucose_uptake",
                                                                                           "beta_oxidation",
                                                                                           "OXPHOS",
                                                                                           "TCA_cycle",
                                                                                           "lactate_oxidation")) &
                                                                               `Cell Type` == "fibroblast")



####investigate role of sex variance for stromal cells across all solid tissues, excluding T cell exhaustion genes
CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1_stromal <- subset(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")
CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1_stromal <- subset(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1_stromal, `Gene Symbol` != "TOX" 
                                                                                   & `Gene Symbol` != "LAG3" & `Gene Symbol` != "PDCD1" & `Gene Symbol` != "CTLA4" & `Gene Symbol` != "HAVCR2")


##take out sex-specific tissues:
CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1_stromal <- subset(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1_stromal, Tissue != "prostate gland" & Tissue != "breast" & Tissue != "testis"
                                                                           & Tissue != "fallopian tube" & Tissue != "uterus" & Tissue != "ovary")

##sex differences are real but smaller than biological differences across genes or tissues:
fit_stromal_cells_all_solid_tissues_by_sex <- lmer(Expression ~ 1 + (1 | Sex) + (1 | `Gene Symbol`) + (1 | Tissue),
                                            data = CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1_stromal)

VarCorr(fit_stromal_cells_all_solid_tissues_by_sex)



####investigate role of sex-mediated variance for T cells across all solid tissues, excluding T cell exhaustion genes
CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1_T <- subset(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1, `Cell Type` == "T cell")

#take out genes not reflective of T cell biology (i.e., stromal specific genes)
CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1_T <- subset(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1_T, Class != "ECM_density")

##take out sex-specific tissues:
CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1_T <- subset(CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1_T, Tissue != "prostate gland" & Tissue != "breast" & Tissue != "testis"
                                                                           & Tissue != "fallopian tube" & Tissue != "uterus" & Tissue != "ovary")


##Variance by sex is real but much smaller than biological differences across genes and/or tissues.
fit_T_cells_all_solid_tissues_by_sex <-  lmer(Expression ~ 1 + (1 | Sex) + (1 | `Gene Symbol`) + (1 | Tissue),
                                      data = CELLxGENE_gene_expression_092525_by_sex_solid_tissues_b1_T)

VarCorr(fit_T_cells_all_solid_tissues_by_sex)
##hence, it is asssumed safe to aggregate each combination of Tissue, Cell Type and Gene Symbol by sex as well as publication.
###this is the default by which CELLxGENE data is downloaded when selecting genes, tissues and cell types of interest.


###inspect relative expressions in the aggregated data (DEFAULT DOWNLOAD FROM CELLXGENE):
##Note that in the script below we do not just calculate the relative (weighted & normalised) expressions as shown in Table 1 of the manuscript,
##but we also explore cross-tissue correlations as suggested in Table 1, as well as tropisms hinted at within the manuscript.

CELLxGENE_gene_expression_092525_T_cell_exhaustion$Class <- "T_cell_exhaustion"
CELLxGENE_gene_expression_092525_beta_oxidation$Class <- "beta_oxidation"
CELLxGENE_gene_expression_092525_glycolysis_and_glucose_uptake$Class <- "glycolysis_and_glucose_uptake"
CELLxGENE_gene_expression_092525_lactate_oxidation$Class <- "lactate_oxidation"
CELLxGENE_gene_expression_092525_OXPHOS$Class <- "OXPHOS"
CELLxGENE_gene_expression_092525_TCA_cycle$Class <- "TCA_cycle"
CELLxGENE_gene_expression_092525_pH_homeostasis$Class <- "pH_homeostasis"
CELLxGENE_gene_expression_092525_ECM_density$Class <- "ECM_density"


CELLxGENE_gene_expression_092525 <- rbind(CELLxGENE_gene_expression_092525_T_cell_exhaustion, CELLxGENE_gene_expression_092525_beta_oxidation, 
                                          CELLxGENE_gene_expression_092525_glycolysis_and_glucose_uptake, CELLxGENE_gene_expression_092525_lactate_oxidation, 
                                          CELLxGENE_gene_expression_092525_OXPHOS, CELLxGENE_gene_expression_092525_TCA_cycle, 
                                          CELLxGENE_gene_expression_092525_pH_homeostasis, CELLxGENE_gene_expression_092525_ECM_density)


CELLxGENE_gene_expression_092525_solid_tissues <- subset(CELLxGENE_gene_expression_092525, Tissue != "blood" 
                                                                        & Tissue != "bone marrow" & Tissue != "lymph node" 
                                                                        & Tissue != "milk" & Tissue != "pleural fluid" & Tissue != "saliva" & Tissue != "yolk sac" & Tissue != "embryo" 
                                                                        & Tissue != "placenta" & Tissue != "vasculature" & Tissue != "mucosa" & Tissue != "lamina propria" & Tissue != "omentum")

CELLxGENE_gene_expression_092525_solid_tissues_a <- subset(CELLxGENE_gene_expression_092525_solid_tissues, `Cell Type` != "aggregated")

##QC 1: enough publications
CELLxGENE_gene_expression_092525_solid_tissues_b <- CELLxGENE_gene_expression_092525_solid_tissues_a %>%
  dplyr::filter(Tissue %in% unique(at_least_2_pubs$Tissue))


genes_multi_class <- CELLxGENE_gene_expression_092525_solid_tissues_b %>%
  dplyr::group_by(`Gene Symbol`) %>%
  dplyr::summarise(n_classes = n_distinct(Class),
            classes   = paste(unique(Class), collapse = ", "),
            .groups = "drop") %>%
  dplyr::filter(n_classes > 1)

genes_multi_class   ##one gene classed in two panels/classes; makes sense as the panels are related


unique(CELLxGENE_gene_expression_092525_solid_tissues_b$Tissue)


# ##if only keeping focus tissues: uncomment if checking across all tissues
# CELLxGENE_gene_expression_092525_solid_tissues_b <- subset(CELLxGENE_gene_expression_092525_solid_tissues_b, Tissue == "colon" | Tissue == "brain" | Tissue == "breast"
#                                                                       | Tissue == "brain" | Tissue == "esophagus" | Tissue == "eye" | Tissue == "fallopian tube" | Tissue == "kidney" | Tissue == "liver"
#                                                                       | Tissue == "lung" | Tissue == "pancreas" | Tissue == "prostate gland" | Tissue == "skin of body" | Tissue == "stomach" | Tissue == "urinary bladder")
# 


##QC 2: at least 10 Cells Expressing Genes:
CELLxGENE_gene_expression_092525_solid_tissues_b1 <- subset(CELLxGENE_gene_expression_092525_solid_tissues_b, `Number of Cells Expressing Genes` >= 10)
min(CELLxGENE_gene_expression_092525_solid_tissues_b1$`Number of Cells Expressing Genes`)

t_cell_data <- subset(CELLxGENE_gene_expression_092525_solid_tissues_b1, `Cell Type` == "T cell")

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
  dplyr::group_by(`Tissue`, `Gene Symbol`, Class) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * `Cell Count`) / sum(`Cell Count`),
    Total_Cell_Count = sum(`Cell Count`),
    .groups = 'drop'
  )


##as explained in the main manuscript, TOX is selected as a representative T cell exhaustion marker because it is a  master transcriptional regulator 
##that enforces and stabilises the exhausted epigenetic program, effectively “locking” T cells into an exhaustion state.
##we keep all other metabolic genes in there to assess TOX co-expression with metabolic genes across tissues
filtered_data_T_exhaustion <- weighted_expression %>%
  dplyr::filter(Class != "ECM_density" & `Gene Symbol` != "SDHA" & !(Class == "T_cell_exhaustion" & `Gene Symbol` != "TOX")) %>% 
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




##if including complete panel:
filtered_data_T_exhaustion_all <- weighted_expression %>%
  dplyr::filter(Class != "ECM_density" & `Gene Symbol` != "SDHA") %>% 
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
tissue_intrinsic_data <- subset(CELLxGENE_gene_expression_092525_solid_tissues_b1, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")


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
  dplyr::group_by(`Tissue`, `Gene Symbol`, Class) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * `Cell Count`) / sum(`Cell Count`),
    Total_Cell_Count = sum(`Cell Count`),
    .groups = 'drop'
  )


weighted_expression_intrinic_1 <- subset(weighted_expression_intrinic, Class != "T_cell_exhaustion")

filtered_data_intrinic_1 <- weighted_expression_intrinic_1 %>%
  dplyr::group_by(`Gene Symbol`) %>%
  dplyr::mutate(
    normalised_Expression =
      Weighted_Expression / max(Weighted_Expression, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

filtered_data_intrinic_1_heatmap <- as.data.table(filtered_data_intrinic_1)

filtered_data_intrinic_1_heatmap <- subset(filtered_data_intrinic_1_heatmap, `Gene Symbol` != "SDHA") ##remove outlier for heatmap

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


##summarise in across genes per gene panel (Class):
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
CELLxGENE_gene_expression_092525_solid_tissues_b <- CELLxGENE_gene_expression_092525_solid_tissues_a %>%
  dplyr::filter(Tissue %in% unique(at_least_2_pubs$Tissue))

CELLxGENE_gene_expression_092525_solid_tissues_b1 <- subset(CELLxGENE_gene_expression_092525_solid_tissues_b, `Number of Cells Expressing Genes` >= 10)
min(CELLxGENE_gene_expression_092525_solid_tissues_b1$`Number of Cells Expressing Genes`)

CELLxGENE_gene_expression_092525_solid_tissues_b1

tissue_intrinsic_data <- subset(CELLxGENE_gene_expression_092525_solid_tissues_b1, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")

tissue_intrinsic_data <- tissue_intrinsic_data %>%
  dplyr::mutate(
    Tissue = case_when(
      Tissue %in% c("bladder organ", "urinary bladder") ~ "bladder",
      TRUE                                              ~ Tissue
    )
  )

weighted_expression_intrinic <- tissue_intrinsic_data %>%
  dplyr::filter(!is.na(Expression), !is.na(`Cell Count`)) %>%
  dplyr::group_by(`Tissue`, `Gene Symbol`, Class) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * `Cell Count`) / sum(`Cell Count`),
    Total_Cell_Count = sum(`Cell Count`),
    .groups = 'drop'
  )


weighted_expression_intrinic_1 <- subset(weighted_expression_intrinic, Class != "T_cell_exhaustion")

filtered_data_intrinic_1 <- weighted_expression_intrinic_1 %>%
  dplyr::group_by(`Gene Symbol`) %>%
  dplyr::mutate(
    normalised_Expression =
      Weighted_Expression / max(Weighted_Expression, na.rm = TRUE)
  ) 


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

###add the T cell exhaustion panel (so beyond just TOX) as a column
all_panels <- intrinsic_wide %>%
  dplyr::inner_join(
    mean_expression_T_exhaustion_all %>%
      dplyr::select(Tissue, T_cell_exhaustion = Weighted_Expression),
    by = "Tissue"
  )

length(unique(all_panels$Tissue))

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
all solid tissues (N = 25) "
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
       aes(x = glycolysis_and_glucose_uptake,
           y = T_cell_exhaustion,
           color = Tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = my_colours) +  
  labs(
    x = "glycolysis_and_glucose_uptake",
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
CELLxGENE_gene_expression_092525_solid_tissues_b1 <- CELLxGENE_gene_expression_092525_solid_tissues_b1 %>%
  dplyr::mutate(
    Tissue = case_when(
      Tissue %in% c("bladder organ", "urinary bladder") ~ "bladder",
      TRUE                                              ~ Tissue
    )
  )

CELLxGENE_gene_expression_092525_solid_tissues_b1 <- subset(CELLxGENE_gene_expression_092525_solid_tissues_b1, Tissue == "colon" | Tissue == "brain" | Tissue == "breast"
                                                                      | Tissue == "esophagus" | Tissue == "eye" | Tissue == "fallopian tube" | Tissue == "kidney" | Tissue == "liver"
                                                                      | Tissue == "lung" | Tissue == "pancreas" | Tissue == "prostate gland" | Tissue == "skin of body" | Tissue == "stomach" | Tissue == "bladder")


t_cell_data <- subset(CELLxGENE_gene_expression_092525_solid_tissues_b1, `Cell Type` == "T cell")


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



##summarise across genes in the panel:
mean_expression_T_exhaustion_all <- filtered_data_T_exhaustion_all %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)   
  )



##TOX only
mean_expression_T_exhaustion <- filtered_data_T_exhaustion %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "TOX",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)   
  )

##Note that the bins shown in Table 1, meant to roughly categorise tissue-level expressions relative to each other, were based on the values calculated here:
mean_expression_T_exhaustion <- mean_expression_T_exhaustion %>%
  dplyr::mutate(
    Weighted_Expression_norm = round(
      Weighted_Expression / max(Weighted_Expression, na.rm = TRUE),
      2
    )
  )


##Next, assess metabolic intensities across tissue intrinsic/stromal cells:
tissue_intrinsic_data <- subset(CELLxGENE_gene_expression_092525_solid_tissues_b1, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")

weighted_expression_intrinic <- tissue_intrinsic_data %>%
  dplyr::filter(!is.na(Expression), !is.na(`Cell Count`)) %>%
  dplyr::group_by(`Tissue`, `Gene Symbol`, Class) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * `Cell Count`) / sum(`Cell Count`),
    Total_Cell_Count = sum(`Cell Count`),
    .groups = 'drop'
  )

weighted_expression_intrinic_1 <- subset(weighted_expression_intrinic, Class != "T_cell_exhaustion")

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

pheatmap(
  heatmap_matrix_intrinic_1,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "blue"))(100),
  main = "Gene expression by fibroblasts + endothelial cells in normal tissue"
)

##summarise in across genes per gene panel (Class):
mean_expression_intrinsic_1 <- filtered_data_intrinic_1 %>%
  dplyr::group_by(Tissue, Class) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)   
  )



str(filtered_data_intrinic_1)


mean_expression_intrinsic_metabolism <- subset(filtered_data_intrinic_1, Class != "T_cell_exhaustion" & Class != "ECM_density") %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)   
  )



mean_expression_intrinsic_2a <- subset(filtered_data_intrinic_1, `Gene Symbol` == "SLC2A1") %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "SLC2A1",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)    
  )

##Note that the bins shown in Table 1, meant to roughly categorise tissue-level expressions relative to each other, were based on the values calculated here:
mean_expression_intrinsic_2a <- mean_expression_intrinsic_2a %>%
  dplyr::mutate(
    Weighted_Expression_norm = round(
      Weighted_Expression / max(Weighted_Expression, na.rm = TRUE),
      2
    )
  )


mean_expression_intrinsic_3 <- subset(filtered_data_intrinic_1, Class == "pH_homeostasis") %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "mean pH homeostasis genes",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)  
  )

##Note that the bins shown in Table 1, meant to roughly categorise tissue-level expressions relative to each other, were based on the values calculated here:
mean_expression_intrinsic_3 <- mean_expression_intrinsic_3 %>%
  dplyr::mutate(
    Weighted_Expression_norm = round(
      Weighted_Expression / max(Weighted_Expression, na.rm = TRUE),
      2
    )
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

##Note that the ratio is dependent on the panels and cell types elected; in the paper, we show an approximation based on ratios deduced when 
##taking the count of stromal cells (endothelial + fibroblasts) expressing the panel of pH homeostasis genes versus the count of T cells expressing the panel of T cell exhaustion genes
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





############VHL###############
###inspecting the linkage between relative SLC2A1 expression and tropism in VHL lesions 

##VHL plot: 
CELLxGENE_gene_expression_092525_1 <- subset(CELLxGENE_gene_expression_092525, `Number of Cells Expressing Genes` >=10)
CELLxGENE_gene_expression_092525_1_GLUT1 <- subset(CELLxGENE_gene_expression_092525_1, `Gene Symbol` == "SLC2A1")

##elect likely cells-of-origin for lesion formation:
plot_df <- CELLxGENE_gene_expression_092525_1_GLUT1 %>%
  dplyr::filter(`Gene Symbol` == "SLC2A1") %>%
  dplyr::filter(
    (`Tissue` == "adrenal gland"  & `Cell Type` == "chromaffin cell") |
      (`Tissue` %in% c("brain", "spinal cord", "eye") &
         `Cell Type` == "connective tissue cell") | 
      (`Tissue` == "liver"         & `Cell Type` == "endothelial cell") |         
      (!`Tissue` %in% c("adrenal gland","brain","spinal cord","eye","liver") &
         `Cell Type` == "progenitor cell")
  )

highlight_tissues <- c("adrenal gland","brain","eye","kidney","pancreas","spinal cord") ##likely sites of lesion formation in VHL disease

plot_df <- plot_df %>%
  dplyr::mutate(
    colour_grp = ifelse(Tissue %in% highlight_tissues, "highlight", "other"),
    Tissue     = reorder(Tissue, Expression)   
  )

ggplot(
  subset(plot_df, `Number of Cells Expressing Genes` > 15),
  aes(
    x = Expression,
    y = Tissue,
    size = log10(`Number of Cells Expressing Genes`),
    colour = colour_grp
  )
) +
  geom_point(alpha = 0.8) +
  scale_colour_manual(
    values = c("highlight" = "coral1", "other" = "black"),
    guide  = "none"
  ) +
  scale_size_continuous(name = "log10(# Cells expressing SLC2A1)") +
  labs(
    x        = "SLC2A1 Expression",
    y        = "Tissue",
    title    = "SLC2A1 (GLUT1) Expression",
    subtitle = "    Brain, spinal cord & eye: connective tissue cell; 
    adrenal gland: chromaffin cell; 
    liver: endothelial cell;
    rest: progenitor cells"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line    = element_line(colour = "black"),
    panel.grid   = element_blank(),
    plot.subtitle = element_text(size = 10)  
  )









