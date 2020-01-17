# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))

# ==============================================================================
# Data
# ==============================================================================
# load metadata
cells = fread("depmap-ids.tsv")

# load Achilles data
achilles = fread("Achilles_gene_dependency.csv", sep = ",", header = TRUE)

# load RNAi data
rnai = fread("D2_combined_gene_dep_scores.csv", sep = ",", header = TRUE)

# ==============================================================================
# Analysis
# ==============================================================================
# only keep prostate lines, transpose so that Genes are rows, cell lines are columns
achilles_prostate = dcast(
    melt(achilles[V1 %in% cells$DepMap_ID, .SD], id.vars = "V1", variable.name = "Gene"),
    Gene ~ V1
)

idx = c(1, which(colnames(rnai) %in% cells$CCLE_Name))
rnai_prostate = rnai[, ..idx]

# clean column names
colnames(achilles_prostate) = c(
    "Gene",
    sapply(
        colnames(achilles_prostate)[-1],
        function(i) {return(cells[DepMap_ID == i, Cell_Line])}
    )
)

colnames(rnai_prostate) = c(
    "Gene",
    sapply(
        colnames(rnai_prostate)[-1],
        function(i) {return(cells[CCLE_Name == i, Cell_Line])}
    )
)

# clean gene names
achilles_prostate[, Gene := gsub(" (.*)$", "", Gene)]
rnai_prostate[, Gene := gsub(" (.*)$", "", Gene)]

# ==============================================================================
# Save data
# ==============================================================================
fwrite(achilles_prostate, "depmap-crispr.tsv", sep = "\t", col.names = TRUE)
fwrite(rnai_prostate, "depmap-rnai.tsv", sep = "\t", col.names = TRUE)
