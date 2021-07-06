# ==============================================================================
# Meta
# ==============================================================================
# Name
# --------------------------------------
# Description: Combine structural similarity coefficient files into a single table
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))

DIR <- list(
	"res" = file.path("..", "..", "results", "2021-06-30_hicrep")
)


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# metadata
meta <- fread("config.tsv", sep = "\t", header = TRUE)
meta <- meta[(Include == "Yes") & (Source == "Primary")]
SAMPLES <- meta[, SampleID]

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Aggregating data")

# create the data table
agg_data <- data.table(
	SampleID = SAMPLES
)

# add a new column for each sample
for (i in 1:length(SAMPLES)) {
	sid <- SAMPLES[i]
	# each column is initially 0 with name `V1`
	agg_data[, V1 := 0]
	# rename the column to match the sample
	colnames(agg_data)[i + 1] <- sid
}

for (i in 1:length(SAMPLES)) {
	s1 <- SAMPLES[i]
	for (j in 1:length(SAMPLES)) {
		s2 <- SAMPLES[j]
		file_name <- file.path(DIR[["res"]], paste(s1, s2, "scc.txt", sep = "."))
		if (file.exists(file_name) && file.info(file_name)$size > 0) {
			# read in the SCC value
			scc_val <- fread(file_name, sep = "\t")
			if (!is.null(scc_val)) {
				# store it in the aggregate data.table
				agg_data[i, j + 1] <- scc_val[1, V1]
			}
		}
	}
}

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
fwrite(
	agg_data,
	file.path(DIR[["res"]], "scc.tsv"),
	sep = "\t",
	col.names = TRUE
)

