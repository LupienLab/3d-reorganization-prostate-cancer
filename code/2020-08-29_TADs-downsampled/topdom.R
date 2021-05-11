TopDom <- function(matrix.data=NULL, bins=NULL, matrix.file=NULL, window.size, outFile=NULL, statFilter=TRUE, sparse = FALSE, verbose = TRUE)
{
    if (verbose) {
        cat("#########################################################################\n")
        cat("Step 0 : Load Data\n")
        cat("#########################################################################\n")
    }
    if ((is.null(matrix.data) || is.null(bins)) && is.null(matrix.file)) {
        stop("`matrix.data` and `bins`, or `matrix.file` must be defined.")
    } else if (!is.null(matrix.file)) {
        matdf = read.table(matrix.file, header=F)
        
        if( ncol(matdf) - nrow(matdf) == 3) {
            colnames(matdf) <- c("chr", "from.coord", "to.coord")
        } else if( ncol(matdf) - nrow(matdf) ==4 ) {
            colnames(matdf) <- c("id", "chr", "from.coord", "to.coord")
        } else {
            cat("Unknown Type of matrix file")
            return(0)
        }    
        bins <- data.frame(
            id=1:n_bins, 
            chr=matdf[, "chr"], 
            from.coord=matdf[, "from.coord"], 
            to.coord=matdf[, "to.coord"]
        )
        matrix.data <- as.matrix( matdf[, (ncol(matdf) - nrow(matdf)+1 ):ncol(matdf)] )
    }
    n_bins = nrow(bins)
    mean.cf <- rep(0, n_bins)
    pvalue <- rep(1, n_bins)

    local.ext = rep(-0.5, n_bins)
    window.size = as.numeric(window.size)

    if (verbose) {
        cat("-- Done!\n")
        cat("Step 0 : Done !!\n")

        cat("#########################################################################\n")
        cat("Step 1 : Generating binSignals by computing bin-level contact frequencies\n")
        cat("#########################################################################\n")
    }
    ptm <- proc.time()
    for(i in 1:n_bins) {
        diamond = Get.Diamond.Matrix(mat.data=matrix.data, i=i, size=window.size)
        mean.cf[i] = mean(diamond)
    }
    
    eltm = proc.time() - ptm
    if (verbose) {
        cat(paste("Step 1 Running Time : ", eltm[3]), "\n")
        cat("Step 1 : Done !!\n")
        
        cat("#########################################################################\n")
        cat("Step 2 : Detect TD boundaries based on binSignals\n")
        cat("#########################################################################\n")
    }
    
    ptm = proc.time()

    gap.idx = Which.Gap.Region(matrix.data=matrix.data, window.size)
    
    proc.regions = Which.process.region(rmv.idx=gap.idx, n_bins=n_bins, min.size=3)
    
    for( i in 1:nrow(proc.regions)) {
        start = proc.regions[i, "start"]
        end = proc.regions[i, "end"]
        
        if (verbose) cat(paste("Process Regions from ", start, "to", end), "\n")
            
        local.ext[start:end] = Detect.Local.Extreme(x=mean.cf[start:end])
    }
    
    eltm = proc.time() - ptm
    if (verbose) {
        cat(paste("Step 2 Running Time : ", eltm[3]), "\n")
        cat("Step 2 : Done !!\n")
        
    }
    if(statFilter)
    {
        if (verbose) {
            cat("#########################################################################\n")
            cat("Step 3 : Statistical Filtering of false positive TD boundaries\n")
            cat("#########################################################################\n")
            cat("-- Matrix Scaling....\n")
        }
    
        ptm = proc.time()
        # SPARSIFY #
        scale.matrix.data = matrix.data
        for( i in 1:(2*window.size) )
        {
            scale.matrix.data[ seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins) ] = scale( matrix.data[ seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins) ] )
        }
        
        if (verbose) cat("-- Compute p-values by Wilcox Ranksum Test\n")
        for( i in 1:nrow(proc.regions))
        {
            start = proc.regions[i, "start"]
            end = proc.regions[i, "end"]
            
            if (verbose) cat(paste("Process Regions from ", start, "to", end), "\n")
            
            pvalue[start:end] <- Get.Pvalue(matrix.data=scale.matrix.data[start:end, start:end], size=window.size, scale=1)
        }
        if (verbose) {
            cat("-- Done!\n")
            cat("-- Filtering False Positives\n")
        }
        
        local.ext[intersect( union(which( local.ext==-1), which(local.ext==-1)), which(pvalue<0.05))] = -2
        local.ext[which(local.ext==-1)] = 0
        local.ext[which(local.ext==-2)] = -1
        if (verbose) cat("-- Done!\n")
        
        eltm = proc.time() - ptm
        if (verbose) {
            cat(paste("Step 3 Running Time : ", eltm[3]), "\n")
            cat("Step 3 : Done!\n")
        }
    } else pvalue = 0

    domains = Convert.Bin.To.Domain.TMP(
        bins=bins, 
        signal.idx=which(local.ext==-1), 
        gap.idx=which(local.ext==-0.5), 
        pvalues=pvalue, 
        pvalue.cut=0.05
    )
    
    bins = cbind(
        bins, 
        local.ext = local.ext,
        mean.cf = mean.cf, 
        pvalue = pvalue
    )
    
    bedform = domains[, c("chr", "from.coord", "to.coord", "tag")]
    colnames(bedform) = c("chrom", "chromStart", "chromEnd", "name")
    
    if( !is.null(outFile) ) {
        if (verbose) {
            cat("#########################################################################\n")
            cat("Writing Files\n")
            cat("#########################################################################\n")
        }
        # output files
        outBinSignal =    paste(outFile, ".binSignal", sep="")
        outDomain = paste(outFile, ".domain", sep="")
        outBed = paste(outFile, ".bed", sep="")
        if (verbose) {
            cat(paste("binSignal File :", outBinSignal) )
            cat(paste("Domain File :", outDomain) )
            cat(paste("Bed File : ", outBed))
        }
        # save files
        write.table(bins, file=outBinSignal, quote=F, row.names=F, col.names=T, sep="\t")
        write.table( domains, file=outDomain, quote=F, row.names=F, col.names=T, sep="\t")
        write.table( bedform, file=outBed, quote=F, row.names=F, col.names=F, sep="\t")
    }
    
    if (verbose) {
        cat("Done!!\n")
        cat("Job Complete !\n")
    }
    
    return(list(binSignal=bins, domain=domains, bed=bedform))
}

#' Get.Diamond.Matrix
#' @param mat.data : N by N matrix, where each element indicate contact frequency
#' @param i :integer, bin index
#' @param size : integer, window size to expand from bin
#' @return : matrix.
Get.Diamond.Matrix <- function(mat.data, i, size)
{
    # SPARSIFY #
    n_bins = nrow( mat.data )
    if(i==n_bins) return(NA)
    
    lowerbound = max( 1, i-size+1 )
    upperbound = min( i+size, n_bins)
    
    return( mat.data[lowerbound:i, (i+1):upperbound] )
}

#' Which.process.region
#' @param rmv.idx : vector of idx, remove index vector
#' @param n_bins : total number of bins
#' @param min.size : minimum size of bins
#' @return : data.frame of proc.regions
Which.process.region <- function(rmv.idx, n_bins, min.size=3)
{
    gap.idx = rmv.idx
    
    proc.regions = data.frame(start=numeric(0), end=numeric(0))
    proc.set = setdiff(1:n_bins, gap.idx)
    n_proc.set = length(proc.set)
    
    i=1
    while(i < n_proc.set )
    {
        start = proc.set[i]
        j = i+1
        
        while(j <= n_proc.set)
        {
            if( proc.set[j] - proc.set[j-1] <= 1) j = j + 1
            else {
                proc.regions = rbind(proc.regions, c(start=start, end=proc.set[j-1]) )
                i = j
                break
            }
        }
        
        if(j >= n_proc.set ) {
            proc.regions = rbind(proc.regions, c(start=start, end=proc.set[j-1]) )
            break
        }
    }
    
    colnames(proc.regions) = c("start", "end")
    proc.regions <- proc.regions[ which( abs(proc.regions[,"end"] - proc.regions[, "start"]) >= min.size ), ]
    
    return(proc.regions)
}


#' Which.Gap.Region
#' @breif version 0.0.2 used
#' @param matrix.data : n by n matrix
#' @return gap index
Which.Gap.Region <- function(matrix.data, w)
{
    # SPARSIFY #
    n_bins = nrow(matrix.data)
    gap = rep(0, n_bins)
    
    for(i in 1:n_bins)
    {
        if( sum( matrix.data[i, max(1, i-w):min(i+w, n_bins)], na.rm = TRUE ) == 0 ) gap[i]=-0.5
    }
    
    idx = which(gap == -0.5)
    return(idx)
}

#' Detect.Local.Extreme
#' @param x : original signal to find local minima
#' @return vector of local extrme, -1 if the index is local minimum, 1 if the index is local maxima, 0 otherwise.
Detect.Local.Extreme <- function(x)
{
    n_bins = length(x)
    ret = rep(0, n_bins)
    x[is.na(x)]=0
    
    if(n_bins <= 3)
    {
        ret[which.min(x)]=-1
        ret[which.max(x)]=1
        
        return(ret)
    }
    # Norm##################################################3
    new.point = Data.Norm(x=1:n_bins, y=x)
    x=new.point$y
    ##################################################
    cp = Change.Point(x=1:n_bins, y=x)
    
    if( length(cp$cp) <= 2 ) return(ret)
    if( length(cp$cp) == n_bins) return(ret)
    for(i in 2:(length(cp$cp)-1))
    {
        if( x[cp$cp[i]] >= x[cp$cp[i]-1] && x[cp$cp[i]] >= x[cp$cp[i]+1] ) ret[cp$cp[i]] = 1
        else if(x[cp$cp[i]] < x[cp$cp[i]-1] && x[cp$cp[i]] < x[cp$cp[i]+1]) ret[cp$cp[i]] = -1
        
        min.val = min( x[ cp$cp[i-1] ], x[ cp$cp[i] ] )
        max.val = max( x[ cp$cp[i-1] ], x[ cp$cp[i] ] )
        
        if( min( x[cp$cp[i-1]:cp$cp[i]] ) < min.val ) ret[ cp$cp[i-1] - 1 + which.min( x[cp$cp[i-1]:cp$cp[i]] ) ] = -1
        if( max( x[cp$cp[i-1]:cp$cp[i]] ) > max.val ) ret[ cp$cp[i-1] - 1 + which.max( x[cp$cp[i-1]:cp$cp[i]] ) ] = 1
    }

    return(ret)
}

#' Data.Norm
#' @param x : x axis vector
#' @param x : y axis vector
#' @return list of normalized x and y
Data.Norm <- function(x, y)
{
    ret.x = rep(0, length(x))
    ret.y = rep(0, length(y))
    
    ret.x[1] = x[1]
    ret.y[1] = y[1]
    
    diff.x = diff(x)
    diff.y = diff(y)
    
    scale.x = 1 / mean( abs(diff(x) ) )
    scale.y = 1 / mean( abs( diff(y) ) )
    
    #print(scale.x)
    #print(scale.y)
    
    for(i in 2:length(x))
    {
        ret.x[i] = ret.x[i-1] + (diff.x[i-1]*scale.x)
        ret.y[i] = ret.y[i-1] + (diff.y[i-1]*scale.y)
    }
    
    return(list(x=ret.x, y=ret.y))
}

#' Change.Point
#' @param x : x axis vector
#' @param x : y axis vector
#' @return change point index in x vector,
#' Note that the first and the last point will be always change point
Change.Point <- function( x, y )
{
    if( length(x) != length(y)) 
    {
        cat("ERROR : The length of x and y should be the same")
        return(0)
    }
    
    n_bins <- length(x)
    Fv <- rep(NA, n_bins)
    Ev <- rep(NA, n_bins)
    cp <- 1
    
    i=1
    Fv[1]=0
    while( i < n_bins )
    {
        j=i+1
        Fv[j] = sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 )
        
        while(j<n_bins)
        {
            j=j+1
            k=(i+1):(j-1)
            Ev[j] = ( sum( abs( (y[j]-y[i] )*x[k] - (x[j] -x[i])*y[k] - (x[i]*y[j]) + (x[j]*y[i]) ) ) / sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 ) )
            Fv[j] = sqrt( (x[j]-x[i])^2 + (y[j] - y[i])^2 ) - ( sum( abs( (y[j]-y[i] )*x[k] - (x[j] -x[i])*y[k] - (x[i]*y[j]) + (x[j]*y[i]) ) ) / sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 ) )
            
            #################################################
            #Not Original Code
            if( is.na(Fv[j]) || is.na(Fv[j-1]) ) {
                j = j-1
                cp <- c(cp, j)
                break
            }
            ####################################################3
            if(Fv[j] < Fv[j-1] ) {
                j = j - 1
                cp <- c(cp, j )
                break
            }
        }
        i=j
    }
    
    cp <- c(cp, n_bins)
    
    return(list(cp=cp, objF=Fv, errF=Ev))
}

#' Get.Pvalue
#' @param matrix.data : matrix
#' @param size : size to extend
#' @param scale : scale parameter if necessary. deprecated parameter
#' @return computed p-value vector
Get.Pvalue <- function( matrix.data, size, scale=1 )
{
    n_bins = nrow(matrix.data)
    pvalue <- rep(1, n_bins)
    
    for( i in 1:(n_bins-1) )
    {
        # SPARSIFY #
        dia = as.vector( Get.Diamond.Matrix2(matrix.data, i, size=size) )
        ups = as.vector( Get.Upstream.Triangle(matrix.data, i, size=size) )
        downs = as.vector( Get.Downstream.Triangle(matrix.data, i, size=size) )

        wil.test = wilcox.test(x=dia*scale, y=c(ups, downs), alternative="less", exact=F)
        pvalue[i] = wil.test$p.value    
        
        #cat(paste(i, "=", wil.test$p.value) )    
    }
    
    pvalue[ is.na(pvalue) ] = 1
    return(pvalue)
}

#' Get.Upstream.Triangle
#' @param mat.data : matrix data
#' @param i : bin index
#' @param size : size of window to extend
#' @return upstream triangle matrix
Get.Upstream.Triangle <- function(mat.data, i, size)
{
    n_bins = nrow(mat.data)
    
    lower = max(1, i-size)
    tmp.mat = mat.data[lower:i, lower:i]
    return( tmp.mat[ upper.tri( tmp.mat, diag=F ) ] )
}

#' Get.Downstream.Triangle
#' @param mat.data : matrix data
#' @param i : bin index
#' @param size : size of window to extend
#' @return downstream triangle matrix
Get.Downstream.Triangle <- function(mat.data, i, size)
{
    n_bins = nrow(mat.data)
    if(i==n_bins) return(NA)
    
    upperbound = min(i+size, n_bins)
    tmp.mat = mat.data[(i+1):upperbound, (i+1):upperbound]
    return( tmp.mat[ upper.tri( tmp.mat, diag=F ) ] )
}

#' Get.Diamond.Matrix2
#' @param mat.data : matrix data
#' @param i : bin index
#' @param size : size of window to extend
#' @return diamond matrix
Get.Diamond.Matrix2 <- function(mat.data, i, size)
{
    n_bins = nrow(mat.data)
    new.mat = matrix(rep(NA, size*size), nrow=size, ncol=size)
    
    for(k in 1:size)
    {
        if(i-(k-1) >= 1 && i < n_bins)
        {
            lower = min(i+1, n_bins)
            upper = min(i+size, n_bins)
            
            new.mat[size-(k-1), 1:(upper-lower+1)] = mat.data[i-(k-1), lower:upper]
        }
    }
    
    return(new.mat)
}
    
#' Convert.Bin.To.Domain
#' @param bins : bin information
#' @param signal.idx : signal index
#' @param signal.idx : gap index
#' @param pvalues : pvalue vector
#' @param pvalue.cut : pvalue threshold
#' @return dataframe storing domain information
Convert.Bin.To.Domain <- function(bins, signal.idx, gap.idx, pvalues=NULL, pvalue.cut=NULL)
{
    n_bins = nrow(bins)
    ret = data.frame(chr=character(0), from.id=numeric(0), from.coord=numeric(0), to.id=numeric(0), to.coord=numeric(0), tag=character(0), size=numeric(0))
    levels( x=ret[, "tag"] ) = c("domain", "gap", "boundary")
    
    rmv.idx = setdiff(1:n_bins, gap.idx)
    proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
    from.coord = bins[proc.region[, "start"], "from.coord"]
    n_procs = nrow(proc.region)
    gap = data.frame(
        chr=rep(bins$chr[1], n_procs),
        from.id=rep(0, n_procs),
        from.coord=from.coord,
        to.id=rep(0, n_procs),
        to.coord=rep(0, n_procs),
        tag=rep("gap", n_procs),
        size=rep(0, n_procs),
        stringsAsFactors=F
    )
    
    rmv.idx = union(signal.idx, gap.idx)
    proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
    n_procs = nrow(proc.region)
    from.coord = bins[proc.region[, "start"], "from.coord"]
    domain = data.frame(
        chr=rep(bins$chr[1], n_procs),
        from.id=rep(0, n_procs),
        from.coord=from.coord,
        to.id=rep(0, n_procs),
        to.coord=rep(0, n_procs),
        tag=rep("domain", n_procs),
        size=rep(0, n_procs),
        stringsAsFactors=F
    )
    
    rmv.idx = setdiff(1:n_bins, signal.idx)
    proc.region = as.data.frame( Which.process.region(rmv.idx, n_bins, min.size=1) )
    n_procs = nrow(proc.region)
    if(n_procs>0)
    {
        from.coord = bins[proc.region[, "start"]+1, "from.coord"]    
        boundary = data.frame(
            chr=rep(bins$chr[1], n_procs),
            from.id=rep(0, n_procs),
            from.coord=from.coord,
            to.id=rep(0, n_procs),
            to.coord=rep(0, n_procs),
            tag=rep("boundary", n_procs),
            size=rep(0, n_procs),
            stringsAsFactors=F
        )
        ret = rbind(ret, boundary)
    }
    
    ret = rbind(gap, domain)
    ret = ret[order(ret[,3]), ]
    
    ret[, "to.coord"] = c(ret[2:nrow(ret), "from.coord"], bins[n_bins, "to.coord"])
    ret[, "from.id"] = match( ret[, "from.coord"], bins[, "from.coord"] )
    ret[, "to.id"] = match(ret[, "to.coord"], bins[, "to.coord"])
    ret[, "size"] = ret[,"to.coord"]-ret[,"from.coord"]
    
    if(!is.null(pvalues) && !is.null(pvalue.cut))
    {
        for(i in 1:nrow(ret))
        {
            if(ret[i, "tag"]=="domain")
            {
                domain.bins.idx = ret[i, "from.id"]:ret[i, "to.id"]
                p.value.constr = which( pvalues[ domain.bins.idx ] < pvalue.cut )
                
                if( length(domain.bins.idx) == length(p.value.constr)) ret[i, "tag"] = "boundary"
            }
        }
    }
    
    new.bdr.set = data.frame(
        chr=character(0),
        from.id=numeric(0),
        from.coord=numeric(0),
        to.id=numeric(0),
        to.coord=numeric(0),
        tag=character(0),
        size=numeric(0)
    )
    stack.bdr = data.frame(
        chr=character(0),
        from.id=numeric(0),
        from.coord=numeric(0),
        to.id=numeric(0),
        to.coord=numeric(0),
        tag=character(0),
        size=numeric(0)
    )
    
    i=1
    while(i <= nrow(ret))
    {
        if( ret[i, "tag"] == "boundary" )
        {
            stack.bdr = rbind(stack.bdr, ret[i, ])
        } else if(nrow(stack.bdr)>0) {
            new.bdr = data.frame(
                chr=bins[1, "chr"], 
                from.id = min( stack.bdr[, "from.id"]), 
                from.coord=min(stack.bdr[, "from.coord"]),
                to.id = max( stack.bdr[, "to.id"]),
                to.coord=max(stack.bdr[, "to.coord"]),
                tag="boundary",
                size=max(stack.bdr[, "to.coord"]) - min(stack.bdr[, "from.coord"])
            )
            new.bdr.set = rbind(new.bdr.set, new.bdr)
            stack.bdr = data.frame(
                chr=character(0),
                from.id=numeric(0),
                from.coord=numeric(0),
                to.id=numeric(0),
                to.coord=numeric(0),
                tag=character(0),
                size=numeric(0)
            )
        }
        i = i + 1
    }
    
    ret = rbind( ret[ ret[, "tag"]!="boundary", ], new.bdr.set )
    ret = ret[order(ret[, "to.coord"]), ]
    
    return(ret)
}


#' Convert.Bin.To.Domain
#' @param bins : bin information
#' @param signal.idx : signal index
#' @param signal.idx : gap index
#' @param pvalues : pvalue vector
#' @param pvalue.cut : pvalue threshold
#' @return dataframe storing domain information
Convert.Bin.To.Domain.TMP <- function(bins, signal.idx, gap.idx, pvalues=NULL, pvalue.cut=NULL)
{
    n_bins = nrow(bins)
    ret = data.frame(
        chr=character(0),
        from.id=numeric(0),
        from.coord=numeric(0),
        to.id=numeric(0),
        to.coord=numeric(0),
        tag=character(0),
        size=numeric(0)
    )
    levels( x=ret[, "tag"] ) = c("domain", "gap", "boundary")
    
    rmv.idx = setdiff(1:n_bins, gap.idx)
    proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
    from.coord = bins[proc.region[, "start"], "from.coord"]
    n_procs = nrow(proc.region)
    gap = data.frame(
        chr=rep(bins$chr[1], n_procs),
        from.id=rep(0, n_procs),
        from.coord=from.coord,
        to.id=rep(0, n_procs),
        to.coord=rep(0, n_procs),
        tag=rep("gap", n_procs),
        size=rep(0, n_procs),
        stringsAsFactors=F
    )
    
    rmv.idx = union(signal.idx, gap.idx)
    proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
    n_procs = nrow(proc.region)
    from.coord = bins[proc.region[, "start"], "from.coord"]
    domain = data.frame(
        chr=rep(bins$chr[1], n_procs),
        from.id=rep(0, n_procs),
        from.coord=from.coord,
        to.id=rep(0, n_procs),
        to.coord=rep(0, n_procs),
        tag=rep("domain", n_procs),
        size=rep(0, n_procs),
        stringsAsFactors=F
    )
    
    rmv.idx = setdiff(1:n_bins, signal.idx)
    proc.region = as.data.frame( Which.process.region(rmv.idx, n_bins, min.size=1) )
    n_procs = nrow(proc.region)
    if(n_procs>0)
    {
        # from.coord = bins[proc.region[, "start"]+1, "from.coord"]
        from.coord = bins$from.coord[proc.region$start + 1]
        boundary = data.frame(
            chr=rep(bins$chr[1], n_procs),
            from.id=rep(0, n_procs),
            from.coord=from.coord,
            to.id=rep(0, n_procs),
            to.coord=rep(0, n_procs),
            tag=rep("boundary", n_procs),
            size=rep(0, n_procs),
            stringsAsFactors=F
        )
        ret = rbind(ret, boundary)
    }

    ret = rbind(gap, domain)
    ret = ret[order(ret[,3]), ]
    
    ret$to.coord = c(
        ret$from.coord[2:nrow(ret)],
        bins$to.coord[n_bins]
    )
    ret$from.id = match(ret$from.coord, bins$from.coord)
    ret$to.id = match(ret$to.coord, bins$to.coord)
    ret$size = ret$to.coord - ret$from.coord
    
    if(!is.null(pvalues) && !is.null(pvalue.cut))
    {
        for(i in 1:nrow(ret))
        {
            if(ret$tag[i] == "domain")
            {
                domain.bins.idx = ret$from.id[i]:ret$to.id[i]
                p.value.constr = which( pvalues[ domain.bins.idx ] < pvalue.cut )
                
                if( length(domain.bins.idx) == length(p.value.constr)) ret$tag[i] = "boundary"
            }
        }
    }
    
    return(ret)
}
