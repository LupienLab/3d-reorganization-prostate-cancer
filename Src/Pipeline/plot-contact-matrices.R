# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("Matrix"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dryhic"))
suppressMessages(library("dplyr"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Plot Hi-C contact maps"
    )
    PARSER$add_argument(
        "raw",
        type = "character",
        help = "Sparse matrix file with raw read counts"
    )
    PARSER$add_argument(
        "norm",
        type = "character",
        help = "Sparse matrix file with normalized read counts"
    )
    PARSER$add_argument(
        "bin_signal",
        type = "character",
        help = "TSV file with normalized bin signal information"
    )
    PARSER$add_argument(
        "domains",
        type = "character",
        help = "TSV file with called domain information"
    )
    PARSER$add_argument(
        "res",
        type = "integer",
        help = "Bin size resolution"
    )
    PARSER$add_argument(
        "-p", "--prefix",
        type = "character",
        help = "Output file prefix",
        default = "contacts"
    )
    PARSER$add_argument(
        "-c", "--chr",
        type = "character",
        help = "Specific chromosome to plot. Otherwise, a plot will be made for each."
    )
    ARGS <- PARSER$parse_args()
}

CHRS = paste0("chr", c(1:22, "X", "Y"))

# update prefix to include resolution
ARGS$prefix = paste0(ARGS$prefix, ".", as.integer(ARGS$res), "bp")

# ==============================================================================
# Functions
# ==============================================================================
#' Plot HiC contact matrix
#'
#' This function takes a HiC contact matrix and creates heatmap-like prepresentation of it
#' @import magrittr
#' @import Matrix
#' @param mat HiC contact matrix (it could be the output of \code{\link{get_contacts_matrix}})
#' @param coord A vector of size two with the start and end coordinates of the desired region to plot
#' @param tads A data.frame where each row is a TAD, the columns are from.id and to.id
#' @param resolution Resolution (bin size) in bp
#' @param transformation Transformation to apply to the contacts prior to the heatmap representation
#' @param color Vector of colors
#' @param sym Logical indicating if the color scale should be symmetrical (around 0)
#' @param trim Numeric proportion of the data that should be 'flattened' on both extremes
#' @param rotate Logical indicating if the matrix should be rotated
#' @param unit_x_axis Numeric unit for the X axis
#' @param label_x_axis Character X axis label
#' @param na.col Color to depict \code{NA} values
#' @param ... Further arguments to be passed to \code{\link{plot}}
#' @export
#' @examples
#' plot(0)

plot_matrix <- function(mat, coord, tads = NULL, resolution,
                        transformation = logfinite,
                        color = colorRampPalette(c("white", "red"))(100),
                        sym = FALSE, trim = .01, rotate = FALSE,
                        unit_x_axis = 1e6,
                        label_x_axis = "Genomic Position / Mbp",
                        na.col = "white",
                        ...){

    # settings

    options(scipen = 999)
    
    # prepare matrix

    rownames(mat) <- colnames(mat) <- rownames(mat) %>% gsub("^.*:", "", .)
    
    lims <- seq(floor(coord[1] / resolution),
                ceiling(coord[2] / resolution)) * resolution
    
    i <- rownames(mat) %in% lims
    
    mat <- as.matrix(mat[i, i])
    mat <- mat[match(lims, rownames(mat)), match(lims, rownames(mat))]
    rownames(mat) <- colnames(mat) <- lims

    if(rotate) mat <- mat[nrow(mat):1,]
    
    # prepare axis info and parameters

    guides <- pretty(x = rownames(mat) %>% as.numeric)

    guides_pos <- data.frame(y = 1:nrow(mat), x = rownames(mat) %>% as.numeric) %>%
        lm(y ~ x, .) %>%
        predict(newdata = data.frame(x = guides))
    
    par(mar = c(4, 0, 0, 0), pty = "s")

    # trimm

    if(trim > 0){

        trim <- as.matrix(mat) %>% c %>% quantile(c(trim / 2, 1 - trim / 2), na.rm = T)
        mat[mat < trim[1]] <- trim[1]
        mat[mat > trim[2]] <- trim[2]

    }
    
    # prepare range of colors

    x <- as.matrix(mat) %>% transformation %>% as.matrix
    if(sym){

        upper <- max(abs(x), na.rm = T)
        lower <- - upper
        
    }else{
        
        lower <- min(x, na.rm = T)
        upper <- max(x, na.rm = T)
        
    }

    # transform scores into colors

    if (max(x, na.rm = T) == min(x, na.rm = T)) {
        x[] <- color[round(length(color) / 2)]
    } else if (!is.finite(lower) || !is.finite(upper) || is.na(lower) || is.na(upper)) {
        x[] <- color[round(length(color) / 2)]
    } else {
        x[] <- color[cut(
                c(x),
                seq(lower, upper, len = length(color) + 1),
                include = T
            )]
    }

    x[is.na(x)] <- na.col
    
    # get limits of genomic region

    range_pos <- as.numeric(rownames(x)) %>% range

    # get matrix dimensions
    
    nr <- nrow(x)
    nc <- ncol(x)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d

    # plot void region

    plot(NA, type="n",
         xlim=c(0, nr), ylim=c(0, nc),
         xlab=label_x_axis,
         ylab="",
         asp=1, axes = F, cex.lab = 1.5, ...)

    # add heatmap

    rasterImage(as.raster(unclass(x)),
                xleft = 0, xright = nc, ybottom = 0, ytop = nr,
                interpolate = FALSE)
    axis(1, at = guides_pos,
         labels = guides / unit_x_axis, cex.axis = 1.5)

    # add TAD annotations

    if (!is.null(tads)) {
        # convert row and column indices to coordinates for mapping
        tad_lines = rbindlist(apply(
            tads,
            1,
            function(row) {
                dt = data.table(
                    x = c(row[1], row[1], row[2]),
                    y = c(row[1], row[2], row[2])
                )
                if (!rotate) {
                    # using nr - X since (0, 0) is the bottom left x coordinate
                    # in the map, but the matrix is plotted with the 1st row/col
                    # being in the top left
                    dt[, y := nr - y]
                }
                return(dt)
            }
        ))

        lines(x = tad_lines$x, y = tad_lines$y)
    }
    
    invisible()
    
}


# ==============================================================================
# Data
# ==============================================================================
bins = fread(ARGS$bin_signal, sep = "\t", header = TRUE)
domains = fread(ARGS$domains, sep = "\t", header = TRUE)
mtx_raw = readMM(ARGS$raw)
mtx_norm = readMM(ARGS$norm)

bin_names = bins[, paste0(chr, ':', from.coord)]
dimnames(mtx_raw) = list(bin_names, bin_names)
dimnames(mtx_norm) = list(bin_names, bin_names)


# ==============================================================================
# Plots
# ==============================================================================
if (!is.null(ARGS$chr)) {
    chrs_to_plot = ARGs$chr
} else {
    chrs_to_plot = CHRS
}

cat("Plotting\n")
for (plot_chr in chrs_to_plot) {
    cat("  ", plot_chr, "\n")
    bin_idx = bins[, which(chr == plot_chr)]
    coord_lims = c(
        bins[bin_idx, ][1, from.coord],
        bins[bin_idx, ][length(bin_idx), from.coord] + ARGS$res
    )
    domains_tuples = domains[chr == plot_chr & tag == "domain", .(from.id, to.id)]
    pdf(
        paste(ARGS$prefix, plot_chr, "raw", "pdf", sep = "."),
        height = 12,
        width = 12
    )
    plot_matrix(
        mat = as.matrix(mtx_raw[bin_idx, bin_idx]),
        coord = coord_lims,
        resolution = ARGS$res
    )
    dev.off()

    pdf(
        paste(ARGS$prefix, plot_chr, "ice-corrected", "pdf", sep = "."),
        height = 12,
        width = 12
    )
    plot_matrix(
        mat = as.matrix(mtx_norm[bin_idx, bin_idx]),
        coord = coord_lims,
        resolution = ARGS$res
    )
    dev.off()

    pdf(
        paste(ARGS$prefix, plot_chr, "ice-corrected-with-TADs", "pdf", sep = "."),
        height = 12,
        width = 12
    )
    plot_matrix(
        mat = as.matrix(mtx_norm[bin_idx, bin_idx]),
        coord = coord_lims,
        resolution = ARGS$res,
        tads = domains_tuples
    )
    dev.off()
}


