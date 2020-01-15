# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))

# set of window sizes upon which TADs are called
W = seq(3, 30)

# number of bins to look up-/down-stream for merging together
K = 1

# contact matrix resolution
BIN_SIZE = 40000

# diretory containing TAD calls
TAD_DIR = file.path("..", "2019-07-08_TADs", "TADs")

# ==============================================================================
# Functions
# ==============================================================================
#' Determines whether a vector of integers is contiguous or not
#'
#' @param x Vector of integers
#' @return TRUE/FALSE
is_contiguous = function(x) {
    bounds = c(min(x), max(x))
    y = seq(bounds[1], bounds[2])
    if (length(x) == length(y)) {
        return(all(y == x[order(x)]))
    }
    return(FALSE)
}

#' Determin if all values in a are greater than all values in b
#'
#' @param a a vector of integers
#' @param b a vector on integers
#' @return TRUE/FALSE
all_greater = function(a, b) {
    return(min(a) > max(b))
}

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata = fread(file.path("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"))

s = metadata[1, get("Sample ID")]
all_boundaries = rbindlist(lapply(
    W,
    function(w) {
        dt = fread(
            file.path(TAD_DIR, paste0("w_", w), paste0("PCa", s, ".40000bp.domains.bed")),
            header = FALSE,
            col.names = c("chr", "start", "end", "type")
        )
        # only consider the boundaries, not the TADs themselves
        boundaries = melt(
            dt,
            id.vars = "chr",
            measure.vars = c("start", "end"),
            value.name = "pos",
            variable.name = "side"
        )
        # remove duplicate boundaries
        boundaries = boundaries[, unique(pos), by = "chr"]
        # add window size for later use
        boundaries[, w := w]
        colnames(boundaries) = c("chr", "pos", "w")
        return(boundaries)
    }
))


# aggregate `w` column into a single list based on coordinates
agg_boundaries = all_boundaries[, lapply(.SD, list), by = c("chr", "pos")]
agg_boundaries = agg_boundaries[order(chr, pos), .SD]

# ==============================================================================
# Analysis
# ==============================================================================
# calculate the "order" of each boundary
agg_boundaries[, Order := lengths(w)]

# 1. Find nearby boundaries that may need to be merged
# -------------------------------------------------
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
    if (i %% 100 == 0) print(i)
    # indices of boundaries that need to be resolved
    idx_to_merge = c()
    if (
        (agg_boundaries[i - 1, pos] >= agg_boundaries[i, pos - K * BIN_SIZE])
        && (length(intersect(agg_boundaries[i - 1, w], agg_boundaries[i, w])) == 0)
    ) {
        idx_to_merge = c(i - 1, i)
    }
    if (
        (agg_boundaries[i + 1, pos] <= agg_boundaries[i, pos + K * BIN_SIZE])
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
conflicts = copy(boundaries_to_merge)
for (i in 2:length(conflicts)) {
    if (i %% 100 == 0) print (i)
    # previous conflicting boundaries
    b_prev = conflicts[[i - 1]]
    # current conflicting boundaries
    b_curr = conflicts[[i]]
    b_merged = sort(union(b_prev, b_curr))
    # if there are boundaries in the previous conflict that are also involved in this conflict
    if (length(intersect(b_curr, b_prev)) > 0) {
        # mark the previous conflict to be ignored
        ignored_conflicts = c(ignored_conflicts, i - 1)
        # bring the previous conflicting boundaries forward
        conflicts[[i]] = b_merged
    }
}

# remove all the conflicts that can be ignored
conflicts = conflicts[-ignored_conflicts]

# 3. Resolve boundary conflicts
# -----------------------------

# keep track of resolved indices that can later be removed from agg_boundaries
resolved_idx_to_remove = c()
agg_copy = copy(agg_boundaries)
for (i in 1:length(conflicts)) {
    if (i %% 100 == 0) print(i)
    idx = conflicts[[i]]
    # if there are 2 boundaries to merge into 1
    if (length(idx) == 2) {
        # local variables for easier readability of the window sizes for each boundary
        w_1 = agg_copy[idx[1], w][[1]]
        w_2 = agg_copy[idx[2], w][[1]]
        w_merged = list(sort(c(w_1, w_2)))
        # if the windows sizes are all smaller in one than the other
        # resolve smaller w boundary to the larger w boundary
        if (all_greater(w_1, w_2)) {
            # remove the (idx[2])-th row from `agg_boundaries`
            resolved_idx_to_remove = c(resolved_idx_to_remove, idx[2])
            # combine set of `w`s together
            agg_copy[idx[1], w := w_merged]
        } else if (all_greater(w_2, w_1)) {
            # remove the (idx[1])-th row from `agg_boundaries`
            resolved_idx_to_remove = c(resolved_idx_to_remove, idx[1])
            # combine set of `w`s together
            agg_copy[idx[2], w := w_merged]
        # if the window sizes are not completely uniform
        # resolve to the boundary with the larger order
        } else if (agg_copy[idx[1], Order] > agg_copy[idx[2], Order]) {
            agg_copy[idx[1], w := w_merged]
            resolved_idx_to_remove = c(resolved_idx_to_remove, idx[2])
        } else if (agg_copy[idx[1], Order] > agg_copy[idx[2], Order]) {
            agg_copy[idx[2], w := w_merged]
            resolved_idx_to_remove = c(resolved_idx_to_remove, idx[1])
        }
    # if there are 3 boundaries to merge, take the middle of the loci
    } else if (length(idx) == 3) {
        w_1 = agg_copy[idx[1], w][[1]]
        w_2 = agg_copy[idx[2], w][[1]]
        w_3 = agg_copy[idx[3], w][[1]]
        w_merged = list(sort(c(w_1, w_2, w_3)))
        resolved_idx_to_remove = c(resolved_idx_to_remove, idx[1], idx[3])
        agg_copy[idx[2], w := w_merged]
    # if more than 3, take the locus with the largest order
    } else {
        idx_with_max_order = which.max(agg_copy[idx, Order])
        idx_to_keep = idx[idx_with_max_order]
        resolved_idx_to_remove = c(resolved_idx_to_remove, idx[-idx_with_max_order])
        w_merged = unique(sort(unlist(agg_copy[idx, w])))
        agg_copy[idx_to_keep, w := w_merged]
    }
}
# remove boundaries that have been resolved to be ignored
agg_copy = agg_copy[-unique(resolved_idx_to_remove), .SD]
# recalculate orders based on resolved merges
agg_copy[, Order := lengths(w)]


# ==============================================================================
# Plots
# ==============================================================================
