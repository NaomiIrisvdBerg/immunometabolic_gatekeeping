##Single cell RNAseq analysis publicly available CELLxGENE data across all tissues for human, normal
##In this script, we first explore the inter Sex sources of variance, in order to determine 
#whether the default download of CELLxGENE (one aggregated Expression value for each unique combination of Cell Type, Gene Symbol and Tissue)
#can reliably be used to investigate cross-tissue expression patterns for a given subset of genes and cell types.
##We then analyse cross-tissue patterns including solid tissues other than those listed in the manuscript, to explore robustness of the patterns
##From line ~962 onwards, we explore trends exclusively in the 14 focus tissues
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




##investigate number of publications per tissue
Feb_2026_by_pub

unique(Feb_2026_by_pub$`Gene Symbol`)


##We exclude: non-solid tissues (circulating or fluid-derived), embryonic or developmental tissues, or tissue subcompartments rather than whole solid organs
Feb_2026_by_pub_solid_tissues <- subset(Feb_2026_by_pub, Tissue != "blood" 
                                                                        & Tissue != "bone marrow" & Tissue != "lymph node" 
                                                                        & Tissue != "milk" & Tissue != "pleural fluid" & Tissue != "saliva" & Tissue != "yolk sac" & Tissue != "embryo" 
                                                                        & Tissue != "placenta" & Tissue != "vasculature" & Tissue != "mucosa" & Tissue != "lamina propria" & Tissue != "omentum")


Feb_2026_by_pub_solid_tissues_a <- subset(Feb_2026_by_pub_solid_tissues, Publication != "aggregated")
Feb_2026_by_pub_solid_tissues_a <- subset(Feb_2026_by_pub_solid_tissues_a, `Number of Cells Expressing Genes` >=10)
min(Feb_2026_by_pub_solid_tissues_a$`Number of Cells Expressing Genes`)



##Issue with publication-aggregated data (i.e., CELLxGENE_gene_expression_020626) is that some high-quality studies are excluded if we directly filter for cell types of interest,
#as they may not have "T cell", "endothelial cell" or "fibroblast" Cell Types registered as such.
##e.g., for Tissue == breast, Reed et al. (2024) Nat Genet report finer-grained resolutions only ("fibroblast of mammary gland" instead of "fibroblast")


##Because CELLxGENE aggregates studies with heterogeneous cell-type annotation granularity, we ensured that each study contributed once per major cell class, regardless of taxonomic depth.
#To avoid over-representing studies with finer annotation granularity, each publication contributed a single value per parent cell class: we used the coarse label when available, 
#otherwise we collapsed finer labels using cell-count weighting
##so first: we collapse within each publication to one value per parent cell class;
#and then second: we aggregate across publications to get tissue-level expression for a given gene and cell type combination of interest




Feb_2026_by_pub_solid_tissues_a <-
  Feb_2026_by_pub_solid_tissues_a %>%
  dplyr::mutate(
    Parent_Cell_Type = dplyr::case_when(
      grepl("fibroblast", `Cell Type`, ignore.case = TRUE) ~ "fibroblast",
      grepl("endothelial", `Cell Type`, ignore.case = TRUE) ~ "endothelial cell",
      grepl("T cell", `Cell Type`, ignore.case = FALSE) ~ "T cell",
      TRUE ~ NA_character_
    )
  )

table(Feb_2026_by_pub_solid_tissues_a$Parent_Cell_Type,
      useNA = "ifany")



#######
##any T cells shared between publications?
t_cells <- subset(Feb_2026_by_pub_solid_tissues_a, Tissue == "breast" | Tissue == "kidney" | Tissue == "bladder organ" 
                  |  Tissue == "brain"  |  Tissue == "colon"  |  Tissue == "esophagus"  |  Tissue == "eye"  |  Tissue == "fallopian tube"
                  |  Tissue == "liver"  |  Tissue == "lung"  |  Tissue == "pancreas" |  Tissue == "prostate gland" |  Tissue == "skin of body"
                  |  Tissue == "stomach" | Tissue == "urinary bladder") %>%
  dplyr::filter(Parent_Cell_Type == "T cell")

celltype_pub_counts <- t_cells %>%
  dplyr::distinct(Publication, `Cell Type`) %>%
  dplyr::count(`Cell Type`, name = "n_publications")


total_pubs_with_T <- t_cells %>%
  dplyr::distinct(Publication) %>%
  nrow()

shared_celltypes <- celltype_pub_counts %>%
  dplyr::filter(n_publications == total_pubs_with_T)
#########


##check if T cells and fibroblasts are now captured for high-res studies like Reed et al.
Feb_2026_by_pub_solid_tissues_a %>%
  dplyr::filter(Publication == "Reed et al. (2024) Nat Genet",
                Tissue == "breast") %>%
  dplyr::distinct(`Cell Type`, Parent_Cell_Type)


##here we identify our cells of focus for this analysis: stromal cells and T cells
focus <- c("fibroblast", "endothelial cell", "T cell")

dat_focus <- Feb_2026_by_pub_solid_tissues_a %>%
  dplyr::filter(Parent_Cell_Type %in% focus)



##Within each Parent_Cell_Type for each Publication, we can have nested structures:
cell_counts_unique <- dat_focus %>%
  dplyr::distinct(
    Publication,
    Tissue,
    Parent_Cell_Type,
    `Cell Type`,
    `Cell Count`
  )


##we use the function below to identify likely nested structures based on summations of Cell Counts:
detect_nesting <- function(df) {
  
  counts <- df$`Cell Count`
  labels <- df$`Cell Type`
  
  ##we need at least 3 labels total to have: target + at least 2 others
  if (length(counts) < 3) return(tibble())
  
  results <- list()
  
  for (i in seq_along(counts)) {
    
    target <- counts[i]
    others <- counts[-i]
    other_labels <- labels[-i]
    
    ##we ofcourse need at least 2 others to form a sum
    if (length(others) < 2) next
    
    for (k in 2:length(others)) {
      
      combs <- combn(seq_along(others), k, simplify = FALSE)
      
      for (c in combs) {
        if (sum(others[c]) == target) {
          results[[length(results) + 1]] <- tibble(
            parent_population = labels[i],
            subpopulations = paste(other_labels[c], collapse = " + "),
            parent_count = target
          )
        }
      }
    }
  }
  
  bind_rows(results)
}

##here we see that there are hierarchical structures of subpopulations per focus cell type of interest; 
#this we need to account for when aggregating across publications, as we do not want to over-represent cells that are present in more than one cluster
nesting_results <- cell_counts_unique %>%
  dplyr::group_by(Publication, Tissue, Parent_Cell_Type) %>%
  group_modify(~ detect_nesting(.x)) %>%
  dplyr::ungroup()


##focus on highest-level (largest) nesting parent per Publication * Tissue * Parent_Cell_Type combo as the representative Cell Type for our focus cells:
rep_from_nesting <- nesting_results %>%
  dplyr::group_by(Publication, Tissue, Parent_Cell_Type) %>%
  dplyr::slice_max(order_by = parent_count, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    Publication, Tissue, Parent_Cell_Type,
    Rep_Cell_Type  = parent_population,
    Rep_Cell_Count = parent_count,
    Rep_Source     = "nesting_max_parent"
  )

rep_fallback <- cell_counts_unique %>%
  dplyr::group_by(Publication, Tissue, Parent_Cell_Type) %>%
  dplyr::mutate(has_parent_label = any(`Cell Type` == Parent_Cell_Type)) %>%
  dplyr::arrange(dplyr::desc(`Cell Type` == Parent_Cell_Type), dplyr::desc(`Cell Count`)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    Publication, Tissue, Parent_Cell_Type,
    Rep_Cell_Type  = `Cell Type`,
    Rep_Cell_Count = `Cell Count`,
    Rep_Source     = dplyr::if_else(has_parent_label, "explicit_parent", "max_cellcount")
  )


##here, we use the highest detected nesting parent when nesting exists; otherwise fall back on the Cell Type within the Parent Cell Type that has the highest count:
rep_celltype <- rep_fallback %>%
  dplyr::left_join(
    rep_from_nesting,
    by = c("Publication", "Tissue", "Parent_Cell_Type"),
    suffix = c("_fallback", "_nest")
  ) %>%
  dplyr::mutate(
    Used_Nesting = !is.na(Rep_Cell_Type_nest),
    Rep_Cell_Type  = dplyr::if_else(Used_Nesting, Rep_Cell_Type_nest, Rep_Cell_Type_fallback),
    Rep_Cell_Count = dplyr::if_else(Used_Nesting, Rep_Cell_Count_nest, Rep_Cell_Count_fallback),
    Rep_Source     = dplyr::if_else(Used_Nesting, Rep_Source_nest, Rep_Source_fallback)
  ) %>%
  dplyr::select(Publication, Tissue, Parent_Cell_Type,
                Rep_Cell_Type, Rep_Cell_Count, Rep_Source, Used_Nesting)


##now filter our larger dataframe for the representative cell populations picked per publication * cell type of focus combo:
dat_rep <- dat_focus %>%
  dplyr::inner_join(
    rep_celltype %>% dplyr::select(Publication, Tissue, Parent_Cell_Type, Rep_Cell_Type),
    by = c("Publication", "Tissue", "Parent_Cell_Type")
  ) %>%
  dplyr::filter(`Cell Type` == Rep_Cell_Type)

##quick sanity check: should be 1 distinct cell type per Pub * Tissue * Parent:
dat_rep %>%
  dplyr::distinct(Publication, Tissue, Parent_Cell_Type, `Cell Type`) %>%
  dplyr::count(Publication, Tissue, Parent_Cell_Type) %>%
  dplyr::summarise(max_n = max(n), min_n = min(n))



##per Parent_Cell_Type, we want to take the representative cell type defined above to be passed on as the representative group for aggregation by publication

agg_by_pub <- dat_rep %>%
  dplyr::group_by(Publication, Tissue, Parent_Cell_Type, `Gene Symbol`) %>%
  dplyr::summarise(
    Expression        = Expression,
    Expression_Scaled = `Expression, Scaled`,
    ##keep the representative counts (constant within group):
    Cell_Count     = dplyr::first(`Cell Count`),
    Num_Cells_Expr = dplyr::first(`Number of Cells Expressing Genes`),
    
    .groups = "drop"
  )


rep_meta <- rep_celltype %>%
  dplyr::count(Parent_Cell_Type, Rep_Source)

##attach per-publication metadata to agg_by_pub for tracability:
agg_by_pub <- agg_by_pub %>%
  dplyr::left_join(
    rep_celltype %>% dplyr::select(Publication, Tissue, Parent_Cell_Type, Rep_Source, Used_Nesting),
    by = c("Publication", "Tissue", "Parent_Cell_Type")
  )



##Next step aggregation: across publications per tissue to get tissue-level expression:
agg_by_tissue <- agg_by_pub %>%
  dplyr::group_by(Tissue, Parent_Cell_Type, `Gene Symbol`) %>%
  dplyr::summarise(
    Expression        = weighted.mean(Expression,        w = Cell_Count, na.rm = TRUE),
    Expression_Scaled = weighted.mean(Expression_Scaled, w = Cell_Count, na.rm = TRUE),
    Total_Cell_Count     = sum(Cell_Count, na.rm = TRUE),
    Total_Num_Cells_Expr = sum(Num_Cells_Expr, na.rm = TRUE),
    n_publications       = n_distinct(Publication),
    .groups = "drop"
  )


# ##This is if we assign equal weight per study:
# agg_by_tissue_equal <- by_pub %>%
#   dplyr::group_by(Tissue, Parent_Cell_Type, `Gene Symbol`) %>%
#   dplyr::summarise(
#     Expression = mean(Expression, na.rm = TRUE),
#     Expression_Scaled = mean(Expression_Scaled, na.rm = TRUE),
#     Total_Cell_Count = sum(Cell_Count, na.rm = TRUE),
#     Total_Num_Cells_Expr = sum(Num_Cells_Expr, na.rm = TRUE),
#     n_publications = n_distinct(Publication),
#     .groups = "drop"
#   )





##check publications per tissue that have any of the (aggregated) cell types for any of the genes of interest:
check_publications <- agg_by_pub %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(n_publications = n_distinct(Publication), .groups = "drop")


##we keep urinary bladder as we will merge it with bladder organ to represent "bladder" tissue:
at_least_2_pubs <- subset(
  check_publications,
  n_publications > 1 | Tissue == "urinary bladder"
)


agg_by_tissue$`Cell Type` <- agg_by_tissue$Parent_Cell_Type
agg_by_tissue_solid_tissues_b <- agg_by_tissue %>%
  dplyr::filter(Tissue %in% unique(at_least_2_pubs$Tissue))


##to arrive at 3 publications representing the tissue of origin for bladder cancer, we here regard "urinary bladder" and "bladder organ" as "bladder":
agg_by_tissue_solid_tissues_b1 <- 
  agg_by_tissue_solid_tissues_b %>%
  dplyr::mutate(Tissue = case_when(
    Tissue %in% c("urinary bladder", "bladder organ") ~ "bladder",
    TRUE ~ Tissue
  ))


##do the same filtering for the non-publication aggregated dataset:
agg_by_pub_solid_tissues_b <- agg_by_pub %>%
  dplyr::filter(Tissue %in% unique(at_least_2_pubs$Tissue))

##to arrive at 3 publications representing the tissue of origin for bladder cancer, we here regard "urinary bladder" and "bladder organ" as "bladder":
agg_by_pub_solid_tissues_b1 <- 
  agg_by_pub_solid_tissues_b %>%
  dplyr::mutate(Tissue = case_when(
    Tissue %in% c("urinary bladder", "bladder organ") ~ "bladder",
    TRUE ~ Tissue
  ))


CELLxGENE_gene_expression_020626_by_publication_solid_tissues_b <- agg_by_pub_solid_tissues_b
CELLxGENE_gene_expression_020626_by_publication_solid_tissues_b1 <- agg_by_pub_solid_tissues_b1


CELLxGENE_gene_expression_020626_by_publication_solid_tissues_b$`Cell Type` <- CELLxGENE_gene_expression_020626_by_publication_solid_tissues_b$Parent_Cell_Type
CELLxGENE_gene_expression_020626_by_publication_solid_tissues_b1$`Cell Type` <- CELLxGENE_gene_expression_020626_by_publication_solid_tissues_b1$Parent_Cell_Type




###is subsetting fro cell tyoe even necessary???
Feb_2026_by_pub_solid_tissues_b_summarised <- 
  subset(CELLxGENE_gene_expression_020626_by_publication_solid_tissues_b1, `Cell Type` == "endothelial cell" | `Cell Type` == "fibroblast" | `Cell Type` == "T cell")  %>%
  dplyr::group_by(`Gene Symbol`, Tissue, `Cell Type`) %>%
  dplyr::summarise(
    median_expr = median(Expression, na.rm = TRUE),
    mad_expr = mad(Expression, na.rm = TRUE),
    n_studies = sum(!is.na(Expression))
  )


Feb_2026_by_pub_solid_tissues_b_summarised_fibroblast <- subset(Feb_2026_by_pub_solid_tissues_b_summarised, `Cell Type` == "fibroblast")
Feb_2026_by_pub_solid_tissues_b_summarised_endothelial_cell <- subset(Feb_2026_by_pub_solid_tissues_b_summarised, `Cell Type` == "endothelial cell")
Feb_2026_by_pub_solid_tissues_b_summarised_T_cell <- subset(Feb_2026_by_pub_solid_tissues_b_summarised, `Cell Type` == "T cell")

Feb_2026_by_pub_solid_tissues_b_summarised_stromal <- subset(Feb_2026_by_pub_solid_tissues_b_summarised, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")
unique(Feb_2026_by_pub_solid_tissues_b_summarised_stromal$`Cell Type`)


###if only keeping focus tissues: 
Feb_2026_by_pub_solid_tissues_b_summarised_stromal_focus <- subset(Feb_2026_by_pub_solid_tissues_b_summarised_stromal, Tissue == "colon" | Tissue == "brain" | Tissue == "breast"
                                                                                                   | Tissue == "brain" | Tissue == "esophagus" | Tissue == "eye" | Tissue == "fallopian tube" | Tissue == "kidney" | Tissue == "liver"
                                                                                                   | Tissue == "lung" | Tissue == "pancreas" | Tissue == "prostate gland" | Tissue == "skin of body" | Tissue == "stomach" | Tissue == "bladder")

##investigating role of sex in explaining variance;
##note that in this exported csv we do not have access to publication info as well (can only choose one grouping variable) 
#so here we filter directly for the cell types of interest, without being able to perform within-publication aggregation 
#of lower-res cells by larger cell type category:

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







###inspect relative expressions in the aggregated data (our publication-aware aggregation as opposed to the defaulted aggregation):



agg_by_tissue

CELLxGENE_gene_expression_020626_solid_tissues <- subset(agg_by_tissue, Tissue != "blood" 
                                                         & Tissue != "bone marrow" & Tissue != "lymph node" 
                                                         & Tissue != "milk" & Tissue != "pleural fluid" & Tissue != "saliva" & Tissue != "yolk sac" & Tissue != "embryo" 
                                                         & Tissue != "placenta" & Tissue != "vasculature" & Tissue != "mucosa" & Tissue != "lamina propria" & Tissue != "omentum")


##QC 1: enough publications
CELLxGENE_gene_expression_020626_solid_tissues_b <- CELLxGENE_gene_expression_020626_solid_tissues %>%
  dplyr::filter(Tissue %in% unique(at_least_2_pubs$Tissue))


unique(CELLxGENE_gene_expression_020626_solid_tissues_b$Tissue)

##we QC-ed this pretty early on but just a sanity check
min(CELLxGENE_gene_expression_020626_solid_tissues_b$Total_Num_Cells_Expr)

t_cell_data <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b, `Cell Type` == "T cell")

##group urinary bladder and bladder organ together:
t_cell_data <- t_cell_data %>%
  dplyr::mutate(
    Tissue = case_when(
      Tissue %in% c("bladder organ", "urinary bladder") ~ "bladder",
      TRUE                                              ~ Tissue
    )
  )

weighted_expression <- t_cell_data %>%
  dplyr::filter(!is.na(Expression), !is.na(Total_Cell_Count)) %>%
  dplyr::group_by(`Tissue`, `Gene Symbol`) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * Total_Cell_Count) / sum(Total_Cell_Count),
    Total_Cell_Count = sum(Total_Cell_Count),
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




##if including complete set of genes - including metabolic (stress) ones:
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
tissue_intrinsic_data <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")


#uncomment if wanting to group urinary bladder and bladder organ together:
tissue_intrinsic_data <- tissue_intrinsic_data %>%
  dplyr::mutate(
    Tissue = case_when(
      Tissue %in% c("bladder organ", "urinary bladder") ~ "bladder",
      TRUE                                              ~ Tissue
    )
  )


weighted_expression_intrinic <- tissue_intrinsic_data %>%
  dplyr::filter(!is.na(Expression), !is.na(Total_Cell_Count)) %>%
  dplyr::group_by(`Tissue`, `Gene Symbol`) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * Total_Cell_Count) / sum(Total_Cell_Count),
    Total_Cell_Count = sum(Total_Cell_Count),
    .groups = 'drop'
  )


##take out T cell exhaustion genes
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

##upper and lower block show tissues with higher intensities across these metabolic intensity/stress genes
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
CELLxGENE_gene_expression_020626_solid_tissues_b <- agg_by_tissue %>%
  dplyr::filter(Tissue %in% unique(at_least_2_pubs$Tissue))

min(CELLxGENE_gene_expression_020626_solid_tissues_b$Total_Num_Cells_Expr)

tissue_intrinsic_data <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b, `Cell Type` == "fibroblast" | `Cell Type` == "endothelial cell")

tissue_intrinsic_data <- tissue_intrinsic_data %>%
  dplyr::mutate(
    Tissue = case_when(
      Tissue %in% c("bladder organ", "urinary bladder") ~ "bladder",
      TRUE                                              ~ Tissue
    )
  )

weighted_expression_intrinic <- tissue_intrinsic_data %>%
  dplyr::filter(!is.na(Expression), !is.na(Total_Cell_Count)) %>%
  dplyr::group_by(`Tissue`, `Gene Symbol`) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * Total_Cell_Count) / sum(Total_Cell_Count),
    Total_Cell_Count = sum(Total_Cell_Count),
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



##investigate correlations between individual panels across all solid tissues:
ggplot(all_panels,
       aes(x = SLC2A1,
           y = T_cell_exhaustion,
           color = Tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = my_colours) +  
  labs(
    x = "SLC2A1",
    y = "T_cell_exhaustion panel",
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
min(CELLxGENE_gene_expression_020626_solid_tissues_b$Total_Num_Cells_Expr)

unique(CELLxGENE_gene_expression_020626_solid_tissues_b$Tissue)

CELLxGENE_gene_expression_020626_solid_tissues_b1 <- CELLxGENE_gene_expression_020626_solid_tissues_b %>%
  dplyr::mutate(
    Tissue = case_when(
      Tissue %in% c("bladder organ", "urinary bladder") ~ "bladder",
      TRUE                                              ~ Tissue
    )
  )

unique(CELLxGENE_gene_expression_020626_solid_tissues_b1$Tissue)

CELLxGENE_gene_expression_020626_solid_tissues_b1 <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b1, Tissue == "colon" | Tissue == "brain" | Tissue == "breast"
                                                            | Tissue == "esophagus" | Tissue == "eye" | Tissue == "fallopian tube" | Tissue == "kidney" | Tissue == "liver"
                                                            | Tissue == "lung" | Tissue == "pancreas" | Tissue == "prostate gland" | Tissue == "skin of body" | Tissue == "stomach" | Tissue == "bladder")


t_cell_data <- subset(CELLxGENE_gene_expression_020626_solid_tissues_b1, `Cell Type` == "T cell")


weighted_expression <- t_cell_data %>%
  dplyr::filter(!is.na(Expression), !is.na(Total_Cell_Count)) %>%
  dplyr::group_by(`Tissue`, `Gene Symbol`) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * Total_Cell_Count) / sum(Total_Cell_Count),
    Total_Cell_Count = sum(Total_Cell_Count),
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


labs4 <- c("Low", "Intermediate",
           "High", "Very high")

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
  dplyr::filter(!is.na(Expression), !is.na(Total_Cell_Count)) %>%
  dplyr::group_by(`Tissue`, `Gene Symbol`) %>%
  dplyr::summarise(
    Weighted_Expression = sum(Expression * Total_Cell_Count) / sum(Total_Cell_Count),
    Total_Cell_Count = sum(Total_Cell_Count),
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

##note the lower block of this heatmap contains clades of high-metabolism tissues, and the upper block contains clades of lower metab tissues
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



mean_expression_intrinsic_2a <- subset(filtered_data_intrinic_1, `Gene Symbol` == "SLC2A1") %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    `Gene Symbol` = "SLC2A1",
    Weighted_Expression = mean(normalised_Expression, na.rm = TRUE),
    Total_Cell_Count = max(Total_Cell_Count, na.rm = TRUE)    
  )


labs4 <- c("Low", "Intermediate",
           "High", "Very high")

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
  "urinary bladder"
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



