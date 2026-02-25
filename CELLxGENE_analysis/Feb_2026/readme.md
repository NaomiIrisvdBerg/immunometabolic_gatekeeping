This folder contains the analysis of the CELLxGENE database (for genes of interest, extracted February 2026) as reported in "Immunometabolic Gatekeeping: How Tissue Metabolism Conditions Tumor Immunity" 
The first script (260212_cellxgene_main.R) contains the analysis as reported in the manuscript, using CELLxGENE's default cross-publication aggregated 'Expression' values.
The second script (260212_cellxgene_publication_aware_aggregation.R) contains analysis of the same file, but performs publication aware aggregation (i.e., not defaulted aggregation)
This second script addresses a key limitation of directly analysing the CELLxGENE publication-aggregated dataset: 
because studies use heterogeneous cell-type annotation granularity, filtering only for coarse labels (e.g., "T cell", "endothelial cell", "fibroblast") can exclude high-quality studies that report only finer-grained annotations 
(such as "fibroblast of mammary gland"), while simultaneously over-representing studies with many subclusters. 
To resolve this, the second script performs a publication-aware, two-stage aggregation:
First, within each publication and tissue, all fine-grained cell types are mapped to a common parent cell class, and a single representative population is selected per parent class 
by preferentially identifying hierarchical nesting based on cell-count summations (using the highest-level parent when present), 
or otherwise falling back to an explicit coarse label or the largest cell-count cluster. 
This ensures that each study contributes exactly one value per major cell class, regardless of annotation depth. 
Second, these per-publication representative populations are aggregated across publications to generate tissue-level gene expression estimates, 
using cell-count-weighted means and tracking contributing study counts, thereby producing comparable and unbiased tissue-level expression summaries 
while preserving inclusion of finely annotated studies.
Both scripts reproduce the inverse relationship between tissue-level effector (i.e., CD8+ T) cell prognosticity and metabolic intensity proxied by metabolic gene expression.
