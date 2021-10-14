# ==============================================================================
# Meta
# ==============================================================================
# combine-cancer-types
# --------------------------------------
# Description: Combine SV calls from multiple cancer types. Stratify into prostate and non-prostate
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("GenomicRanges"))


# ==============================================================================
# Functions
# ==============================================================================
#' Re-pair SVs based on their IDs and location
#'
#' @param unpaired data.table of unpaired SV breakpoints
#' @return data.table
repair_svs <- function(unpaired) {
	rbindlist(lapply(
		unpaired[, unique(id)],
		function(sv_id) {
			from <- unpaired[(id == sv_id) & (location = "from"), .SD]
			to <- unpaired[(id == sv_id) & (location = "to"), .SD]
			data.table(
				"chr_from" = from[, chr],
				"chr_from_bkpt" = from[, pos],
				"chr_from_strand" = from[, strand],
				"chr_to" = to[, chr],
				"chr_to_bkpt" = to[, pos],
				"chr_to_strand" = to[, strand],
				"sv_id" = sv_id,
				"variant_type" = from[, variant_type]
			)
		}
	))
}

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load metadata
meta <- fread("config.tsv", sep = "\t", header = TRUE)

svs_hg19_paired <- rbindlist(mapply(
	function(proj_code, type_code, subtype_code) {
		dt1 <- fread(
			paste("structural_somatic_mutation", proj_code, "tsv", sep = "."),
			sep = "\t",
			header = TRUE
		)
		dt1[, `:=`(
			project_code   = proj_code,
			cancer_type	= type_code,
			cancer_subtype = subtype_code
		)]
		return(dt1)
	},
	meta[, project_code],
	meta[, cancer_type],
	meta[, cancer_subtype],
	SIMPLIFY = FALSE
))

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# add "chr" prefix to chr_from and chr_to and convert strands
svs_hg19_paired[, `:=`(
	chr_from = paste0("chr", chr_from),
	chr_to = paste0("chr", chr_to),
	chr_from_strand = ifelse(chr_from_strand == 1, "+", "-"),
	chr_to_strand = ifelse(chr_to_strand == 1, "+", "-")
)]

# make each breakpoint a row
svs_hg19_unpaired <- data.table(
	chr = c(svs_hg19_paired$chr_from, svs_hg19_paired$chr_to),
	pos = c(svs_hg19_paired$chr_from_bkpt, svs_hg19_paired$chr_to_bkpt),
	strand = c(svs_hg19_paired$chr_from_strand, svs_hg19_paired$chr_to_strand),
	id = rep(svs_hg19_paired$sv_id, 2),
	variant_type = rep(svs_hg19_paired$variant_type, 2),
	cancer_type = rep(svs_hg19_paired$cancer_type, 2),
	location = rep(
		c("from", "to"),
		c(svs_hg19_paired[, .N], svs_hg19_paired[, .N])
	)
)

# split into prostate and non-prostate SVs
non_prostate_hg19_paired <- svs_hg19_paired[cancer_type != "prostate"]
prostate_hg19_paired <- svs_hg19_paired[cancer_type == "prostate"]
non_prostate_hg19_unpaired <- svs_hg19_unpaired[cancer_type != "prostate"]
prostate_hg19_unpaired <- svs_hg19_unpaired[cancer_type == "prostate"]

# extract useful columns
useful_columns <- c(
	"chr_from",
	"chr_from_bkpt",
	"chr_from_strand",
	"chr_to",
	"chr_to_bkpt",
	"chr_to_strand",
	"sv_id",
	"variant_type"
)
non_prostate_hg19_paired <- non_prostate_hg19_paired[, .SD, .SDcols = useful_columns]
prostate_hg19_paired <- prostate_hg19_paired[, .SD, .SDcols = useful_columns]

# drop the cancer_type column
non_prostate_hg19_unpaired[, cancer_type := NULL]
prostate_hg19_unpaired[, cancer_type := NULL]


# convert to GRanges objects for conversion
non_prostate_hg19_gr <- non_prostate_hg19_unpaired[, GRanges(
	seqnames = chr,
	ranges = IRanges(
		start = pos,
		width = 1
	),
	strand = strand,
	id = id,
	variant_type = variant_type,
	location = location
)]
prostate_hg19_gr <- prostate_hg19_unpaired[, GRanges(
	seqnames = chr,
	ranges = IRanges(
		start = pos,
		width = 1
	),
	strand = strand,
	id = id,
	variant_type = variant_type,
	location = location
)]

chain <- import.chain(file.path("..", "hg19ToHg38.over.chain"))


# perform liftover conversion
non_prostate_grl <- liftOver(non_prostate_hg19_gr, chain)
prostate_grl <- liftOver(prostate_hg19_gr, chain)

# some lifted over regions are separated by a few bases, so exclude them
prostate_hg38 <- unlist(prostate_grl)
prostate_hg38_unpaired <- as.data.table(prostate_hg38)

prostate_pairing_check <- dcast(
	prostate_hg38_unpaired[, .N, by = c("id", "location")],
	id ~ location,
	value.var = "N"
)
uniquely_mapping_ids_prostate <- prostate_pairing_check[(from == 1) & (to == 1), id]
prostate_hg38_unpaired <- prostate_hg38_unpaired[id %in% uniquely_mapping_ids_prostate]

non_prostate_hg38 <- unlist(non_prostate_grl)
non_prostate_hg38_unpaired <- as.data.table(non_prostate_hg38)

non_prostate_pairing_check <- dcast(
	non_prostate_hg38_unpaired[, .N, by = c("id", "location")],
	id ~ location,
	value.var = "N"
)
uniquely_mapping_ids_prostate <- non_prostate_pairing_check[(from == 1) & (to == 1), id]
non_prostate_hg38_unpaired <- non_prostate_hg38_unpaired[id %in% uniquely_mapping_ids_prostate]

# re-pair regions based on SV IDs
prostate_hg38_paired <- merge(
	x = prostate_hg38_unpaired[location == "from"],
	y = prostate_hg38_unpaired[location == "to"],
	by = "id",
	suffix = c("_from", "_to")
)
non_prostate_hg38_paired <- merge(
	x = non_prostate_hg38_unpaired[location == "from"],
	y = non_prostate_hg38_unpaired[location == "to"],
	by = "id",
	suffix = c("_from", "_to")
)

# remap column names
prostate_hg38_paired <- prostate_hg38_paired[, .(
	chr_from = seqnames_from,
	chr_from_bkpt = start_from,
	chr_from_strand = strand_from,
	chr_to = seqnames_to,
	chr_to_bkpt = start_to,
	chr_to_strand = strand_to,
	id = id,
	variant_type = variant_type_from
)]
non_prostate_hg38_paired <- non_prostate_hg38_paired[, .(
	chr_from = seqnames_from,
	chr_from_bkpt = start_from,
	chr_from_strand = strand_from,
	chr_to = seqnames_to,
	chr_to_bkpt = start_to,
	chr_to_strand = strand_to,
	id = id,
	variant_type = variant_type_from
)]

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")

# save in tables
fwrite(
	prostate_hg19_paired,
	"sv-breakpoints.prostate.hg19.paired.tsv",
	sep = "\t",
	col.names = TRUE
)
fwrite(
	non_prostate_hg19_paired,
	"sv-breakpoints.non-prostate.hg19.paired.tsv",
	sep = "\t",
	col.names = TRUE
)
fwrite(
	prostate_hg19_unpaired,
	"sv-breakpoints.prostate.hg19.unpaired.tsv",
	sep = "\t",
	col.names = TRUE
)
fwrite(
	non_prostate_hg19_unpaired,
	"sv-breakpoints.non-prostate.hg19.unpaired.tsv",
	sep = "\t",
	col.names = TRUE
)

fwrite(
	prostate_hg38_paired,
	"sv-breakpoints.prostate.hg38.paired.tsv",
	sep = "\t",
	col.names = TRUE
)
fwrite(
	non_prostate_hg38_paired,
	"sv-breakpoints.non-prostate.hg38.paired.tsv",
	sep = "\t",
	col.names = TRUE
)
fwrite(
	prostate_hg38_unpaired,
	"sv-breakpoints.prostate.hg38.unpaired.tsv",
	sep = "\t",
	col.names = TRUE
)
fwrite(
	non_prostate_hg38_unpaired,
	"sv-breakpoints.non-prostate.hg38.unpaired.tsv",
	sep = "\t",
	col.names = TRUE
)
