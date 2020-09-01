# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Resolve TAD calls across window sizes"
    )
    PARSER$add_argument(
        "id",
        type = "character",
        help = "Sample ID of patient to be processed"
    )
    PARSER$add_argument(
        "count",
        type = "character",
        help = "Subsampling count to be considered"
    )
    PARSER$add_argument(
        "res",
        type = "integer",
        help = "Resolution of contact matrix"
    )
    PARSER$add_argument(
        "-p", "--prefix",
        type = "character",
        help = "Prefix for output files (`{prefix}.aggregated-domains.tsv` and `{prefix}.aggregated-boundaries.bedGraph` will be created)",
        default = "output"
    )
    PARSER$add_argument(
        "-i", "--in-dir",
        type = "character",
        help = "Parent directory for TAD calls",
        dest = "in_dir"
    )
    PARSER$add_argument(
        "-m", "--min",
        type = "integer",
        help = "Minimum window size to consider",
        default=3
    )
    PARSER$add_argument(
        "-M", "--max",
        type = "integer",
        help = "Maximum window size to consider",
        default=20
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        id = "PCa13266",
        prefix = "Aggregated-TADs/output",
        count = "300000000",
        res = 40000,
        in_dir = "TADs",
        min = 3,
        max = 20
    )
}

# set of window sizes upon which TADs are called
W = seq(ARGS$min, ARGS$max)

# number of bins to look up-/down-stream for merging together
K = 1


# ==============================================================================
# Functions
# ==============================================================================
#' Determin if all values in a are greater than all values in b
#'
#' @param a a vector of integers
#' @param b a vector on integers
#' @return TRUE/FALSE
all_greater = function(a, b) {
    return(min(a) > max(b))
}

list_to_str <- function(l) {
    ifelse(is.na(l) || length(l) == 0, NA, paste(sort(unlist(l)), collapse = ","))
}

# ==============================================================================
# Data
# ==============================================================================
cat("Loading data\n")

# identify all TAD boundaries from all TAD calls across window size parameters
all_boundaries <- rbindlist(lapply(
    W,
    function(w) {
        tads <- fread(
            file.path(ARGS$in_dir, paste0(ARGS$id, ".", ARGS$count, ".res_", ARGS$res, "bp.window_", w, ".domains.bed")),
            header = FALSE,
            col.names = c("chr", "start", "end", "type")
        )
        # only consider the boundaries, not the TADs themselves
        boundaries <- tads[,
            .(
                chr = c(chr, chr),
                pos = c(start, end),
                left_type = c(shift(type, 1), type),
                right_type = c(type, shift(type, 1, type = "lead"))
            )
        ]
        # remove duplicate boundaries
        boundaries <- unique(boundaries)
        # add window size for later use
        boundaries[, w := w]
        return(boundaries)
    }
))

# aggregate `w` column into a single list based on coordinates
agg_boundaries <- all_boundaries[,
    .(
        w = list(unique(w)),
        left_type = list(unique(left_type)),
        right_type = list(unique(right_type))
    ),
    keyby = c("chr", "pos")
]

cat("\t", agg_boundaries[, .N], " unique TAD boundaries\n", sep = "")

# ==============================================================================
# Analysis
# ==============================================================================
# calculate the "persistence" of each boundary
agg_boundaries[, Persistence := lengths(w)]

# 1. Find nearby boundaries that may need to be merged
# -------------------------------------------------
cat("Identifying similar boundaries called with different window sizes\n")

boundaries_to_merge = list()
# keep track of the number of resolutions that need to happen
n_merges = 0
# whether to skip ahead, if the next boundary needs to be resolved with the current boundary,
# since there's not need to check that one twice
skip = FALSE
for (i in 2:agg_boundaries[, .N - 1]) {
    # skip if this boundary was already marked to be resolved with the previous one
    if (skip) {
        skip = FALSE
        next
    }
    # indices of boundaries that need to be resolved
    idx_to_merge = c()
    if (
        (agg_boundaries[i - 1, pos] >= agg_boundaries[i, pos - K * ARGS$res])
        && (length(intersect(agg_boundaries[i - 1, w], agg_boundaries[i, w])) == 0)
    ) {
        idx_to_merge = c(i - 1, i)
    }
    if (
        (agg_boundaries[i + 1, pos] <= agg_boundaries[i, pos + K * ARGS$res])
        && (length(intersect(agg_boundaries[i + 1, w], agg_boundaries[i, w])) == 0)
    ) {
        idx_to_merge = sort(unique(c(idx_to_merge, i, i + 1)))
        skip = TRUE
    }
    if (length(idx_to_merge) > 0) {
        n_merges = n_merges + 1
        boundaries_to_merge[[n_merges]] = idx_to_merge
        i = i + skip
    }
}

# 2.Ensure that each conflict contains unique row indices
# -------------------------------------------------------

# e.g. at this point, one conflict might contain boundaries [i, i+1]
# and the next conflict might contain boundaries [i + 1, i + 2, i + 3]
# these should all be considered as a single conflict, not as two independent conflicts
ignored_conflicts = c()
conflicting_merges = copy(boundaries_to_merge)
for (i in 2:length(conflicting_merges)) {
    # previous conflicting boundaries
    b_prev = conflicting_merges[[i - 1]]
    # current conflicting boundaries
    b_curr = conflicting_merges[[i]]
    b_merged = sort(union(b_prev, b_curr))
    # if there are boundaries in the previous conflict that are also involved in this conflict
    if (length(intersect(b_curr, b_prev)) > 0) {
        # mark the previous conflict to be ignored
        ignored_conflicts = c(ignored_conflicts, i - 1)
        # bring the previous conflicting boundaries forward
        conflicting_merges[[i]] = b_merged
    }
}

# remove all the conflicts that can be ignored
conflicting_merges = conflicting_merges[-ignored_conflicts]

cat("\t", length(conflicting_merges), " conflicts to resolve\n", sep = "")

# 3. Resolve boundary conflicts
# -----------------------------
cat("Resolving conflicting boundary calls\n")

# keep track of resolved indices that can later be removed from agg_boundaries
resolved_idx_to_remove = c()
agg_copy = copy(agg_boundaries)
for (i in 1:length(conflicting_merges)) {
    # indices of boundaries that are in conflict
    idx <- conflicting_merges[[i]]
    # merged set of window sizes that will eventually be passed to the "resolved boundary"
    w_merged <- list(unique(sort(unlist(agg_copy[idx, w]))))
    # if there are 2 boundaries to merge into 1
    if (length(idx) == 2) {
        # local variables for easier readability of the window sizes for each boundary
        w_1 = agg_copy[idx[1], w][[1]]
        w_2 = agg_copy[idx[2], w][[1]]
        # if the windows sizes are all smaller in one than the other
        # resolve smaller w boundary to the larger w boundary
        if (all_greater(w_1, w_2)) {
            # remove the (idx[2])-th row from `agg_boundaries`
            which_idx_to_keep = 1
            # combine set of `w`s together
        } else if (all_greater(w_2, w_1)) {
            # remove the (idx[1])-th row from `agg_boundaries`
            which_idx_to_keep = 2
        # if the window sizes are not completely uniform
        # resolve to the boundary with the larger persistence
        } else {
            which_idx_to_keep = which.max(agg_copy[idx, Persistence])
        }
    # if there are 3 boundaries to merge, take the middle of the loci
    } else if (length(idx) == 3) {
        which_idx_to_keep = 2
    # if more than 3, take the locus with the largest persistence
    } else {
        which_idx_to_keep = which.max(agg_copy[idx, Persistence])
    }
    kept_index = idx[which_idx_to_keep]
    resolved_idx_to_remove = c(resolved_idx_to_remove, idx[-which_idx_to_keep])
    agg_copy[kept_index, w := list(w_merged)]
}
# recalculate persistences based on resolved merges
agg_copy[, Persistence := lengths(w)]
# remove boundaries that have been resolved to be ignored
agg_copy = agg_copy[-unique(resolved_idx_to_remove), .SD]

cat("\tResolved to", agg_copy[, .N], "unique TAD boundaries\n")
cat("\tBoundary persistence:\n")
print(agg_copy[, summary(Persistence)])

# 4. Construct hierarchical TADs across orders
# ------------------------------------------------------------------
cat("Constructing hierarchical TADs\n")

agg_tads = rbindlist(lapply(
    rev(W),
    function(level) {
        print(level)
        # get all boundaries that are a boundary at this window size
        idx = which(agg_copy[, sapply(w, function(v) level %in% v)])
        w_boundaries = agg_copy[idx, .SD]
        # produce TADs from consecutive boundaries
        w_tads = rbindlist(lapply(
            2:w_boundaries[, .N],
            function(i) {
                prev = w_boundaries[i - 1]
                curr = w_boundaries[i]
                # skip making this TAD if chromosomes aren't equal
                # empty data.table gets removed in rbindlist
                if (curr$chr != prev $chr) {
                    return(data.table())
                }
                # create record from previous and current boundary
                dt = data.table(
                    chr = curr$chr,
                    start = prev$pos,
                    end = curr$pos,
                    persistence_left = prev$Persistence,
                    persistence_right = curr$Persistence,
                    type = list_to_str(intersect(prev$right_type, curr$left_type))
                )
                return(dt)
            }
        ))
        w_tads[, w := level]
        return(w_tads)
    }
))

# convert window column to a comma-separated string, instead of a list, before saving
agg_copy[, w := unlist(lapply(w, list_to_str))]
agg_copy[, left_type := unlist(lapply(left_type, list_to_str))]
agg_copy[, right_type := unlist(lapply(right_type, list_to_str))]

# ==============================================================================
# Save data
# ==============================================================================
cat("Saving data\n")

# save aggregated TADs in bedGraph format, where the value column is the depth
fwrite(
    agg_tads,
    paste0(ARGS$prefix, ".agg-domains.tsv"),
    col.names = TRUE,
    sep = "\t"
)

# save aggregated boundaries with their order and window sizes
fwrite(
    agg_copy,
    paste0(ARGS$prefix, ".agg-boundaries.tsv"),
    col.names = TRUE,
    sep = "\t"
)

cat("Done\n")

