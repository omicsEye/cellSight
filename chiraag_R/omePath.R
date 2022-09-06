library(omePath)

all_results <- 
  read.table("~/../Box/snRNA_CellRanger_Wound_nonWound/analysis/Tweedieverse_output_Pdgfrb/all_results.tsv",
             header = TRUE) %>%
  filter(pval < .05)
sig_results <- 
  read.table("~/../Box/snRNA_CellRanger_Wound_nonWound/analysis/Tweedieverse_output_Pdgfrb/significant_results.tsv",
           header = TRUE)

metagen <- read_tsv("data/GO_UNIREF90_MAP.tsv")
msig <- read_tsv("data/c5.all.v7.1.entrez.tsv")
metabo <- read_tsv("data/smpdb_metabolites.tsv")

rownames(all_results) <- all_results$feature |> toupper()

omePath_results <- omePath(all_results,
                           "data/omePath_output",
                           msig, 
                           pathway_col = "Pathway",
                           feature_col = "symbol",
                           input_metadata = NA,
                           meta = NA,
                           case_label = NA,
                           control_label = NA,
                           score_col = 'coef',
                           pval_threshold = 0.05,
                           fdr_threshold = NA,
                           Pathway.Subject = NA,
                           method = 'ks',
                           min_member = 2,
                           do_plot = TRUE)

