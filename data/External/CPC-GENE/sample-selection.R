# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
# suppressMessages(library("ggplot2"))
suppressMessages(library("UpSetR"))

# ==============================================================================
# Data
# ==============================================================================
wgs = fread("wgs.tsv")
wgs[, WGS := 1]

rna = fread("rna.tsv")
rna[, RNA := 1]

prot = fread("proteome.tsv")
prot[, Proteome := 1]

chip = fread("h3k27ac.tsv")
chip[, H3K27ac := 1]

T2E = wgs[ETS_Consensus == "Positive", SampleID]
nonT2E = wgs[ETS_Consensus == "Negative", SampleID]
naT2E = wgs[is.na(ETS_Consensus), SampleID]

# merge all data together
agg = merge(wgs, rna, all = TRUE)
agg = merge(agg, prot, all = TRUE)
agg = merge(agg, chip, all = TRUE)
agg[is.na(WGS), WGS := 0]
agg[is.na(RNA), RNA := 0]
agg[is.na(RNA), RNA := 0]
agg[is.na(Proteome), Proteome := 0]
agg[is.na(Proteome), Proteome := 0]
agg[is.na(H3K27ac), H3K27ac := 0]

# ==============================================================================
# Functions
# ==============================================================================
is_t2e = function(row, set) {
    data <- row["SampleID"] %in% set
}

# ==============================================================================
# Plots
# ==============================================================================
png("data-integration.png", width = 12, height = 12, units = "cm", res = 300)
p = upset(
    agg[!is.na(ETS_Consensus)],
    query.legend = "top",
    queries = list(
        list(query.name = "T2E+", params = list(T2E), query = is_t2e, active = TRUE)
    )
)
dev.off()
