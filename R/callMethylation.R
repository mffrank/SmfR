# Functions to call methylation events from alignment files using QuasR

#' Call methylation events on all cytosines
#' @param proj qProject object
#' @param range GRanges object with ONE!!!! range
#' @param sample name of the Genome alignment
#' @importFrom QuasR qMeth
#' @import data.table
#' @import Matrix
#' @import GenomicRanges
#' @return list with GRanges of Cytosine positions and Matrix with methylation events
#' \itemize{
#'   \item GRanges object describing the locations of called cytosines with some metadata for each site:
#'   \describe{
#'     \item{strand}{of the cytosine (strands can later be collapsed)}
#'     \item{met_rate}{Number of methylated sites divided by all sites that are called}
#'   }
#'   \item Matrix with the methylation calls for each read and genomic position. Uncalled=NA
#' }
#' @export

getCMethMatrix<-function(proj,range,sample){
  Cs = qMeth(proj[sample], query=range,mode="allC",reportLevel="alignment")
  # use data.table to get a 1,0 matrix of methylation profiles
  # make the data.table object
  # allmet = data.table(meth=Cs[[sample]]$meth,
  #                   aid=Cs[[sample]]$aid,
  #                   cid=Cs[[sample]]$Cid,
  #                   strand=Cs[[sample]]$strand)
  # dtm = dcast(allmet, aid~cid+strand, value.var="meth")
  # Make a sparse matrix to save on memory
  aid = factor(Cs[[sample]]$aid)
  cid = factor(Cs[[sample]]$Cid)
  if(length(aid)==0|length(cid)==0) return(list(met_locs = GRanges(), met_mat = NULL))
  allmet = Matrix::sparseMatrix(i = as.numeric(aid),
                              j = as.numeric(cid),
                              x = Cs[[sample]]$meth + 1,
                              dimnames = list(levels(aid), as.character(levels(cid))),
                              index1 = TRUE)
  # ronames = dtm$aid
  # dtm[,aid := NULL] # remove unwanted row
  # CpGm = as.matrix(dtm)
  # Get the strand and localization of each methylation event
  # strands = gsub("^.*_","",colnames(CpGm))
  # locs = as.numeric(gsub("_.*","",colnames(CpGm)))
  # Every C should have a uniqe strand,
  # therefore we can just get the first strand we find for each C position
  locs = as.numeric(levels(cid))
  strands = Cs[[sample]]$strand[match(locs, Cs[[sample]]$Cid)]
  # Cleanup methylation matrix
  # colnames(CpGm) <- locs
  # rownames(CpGm)=ronames
  # Produce GRanges object with descriptive statistics
  ran = GRanges(seqnames = seqnames(range),
                ranges = IRanges(start = locs, width = 1),
                strand = strands)
                # coverage = diff(allmet@p),
                # met_rate = colMeans_drop0(allmet) - 1)
  res <- list(met_locs = ran, met_mat = allmet)
  res <- calculateMetStats(res)
  return(res)
}

calculateMetStats <- function(met){
  met$met_locs$coverage <- diff(met$met_mat@p)
  met$met_locs$met_rate <- colMeans_drop0(met$met_mat) - 1
  return(met)
}

collapseStrands <- function(met, collapseCpG = TRUE, collapseGPC = TRUE){
  met_locs <- met$met_locs
  # Collapse CpG
  if(collapseCpG){
    if(!('cg_motif' %in% names(met$met_locs@elementMetadata))){
      met_locs <- addType(met_locs, motifs = c(cg_motif=DNAString("NCG")))
    }
    cg = met_locs$cg_motif
    met <- .collapseStrands(met_locs, met$met_mat, filter = cg, shift = -1)
  }

  met_locs <- met$met_locs
  # Collapse GpC
  if(collapseGPC){
    if(!('gc_motif' %in% names(met$met_locs@elementMetadata))){
      met_locs <- addType(met_locs, motifs = c(gc_motif=DNAString("GCN")))
    }
    gc <- met_locs$gc_motif
    met <- .collapseStrands(met_locs, met$met_mat, filter = gc, shift=1)
  }

  #filt <- cg | gc
  #met <- .collapseStrands(met_locs, met$met_mat, filter = filt)
  return(met)
}

.collapseStrands <- function(met_locs, met_mat, filter, shift = 1){
  #shift <- match.arg(shift, several.ok = F)

  minus_strand <- ((strand(met_locs) == '-') & filter)
  # Move all minus strands to the left
  start(met_locs[minus_strand]) <- start(met_locs[minus_strand]) + shift
  end(met_locs[minus_strand]) <- end(met_locs[minus_strand]) + shift
  strand(met_locs[minus_strand]) <- '+'
  # Remove the sites that have a detected event on the plus strand
  if(shift == -1){
    collapse <- c(FALSE, diff(start(met_locs)) == 0)
  }else{
    collapse <- c(diff(start(met_locs)) == 0, FALSE)
  }
  met_locs <- met_locs[!collapse]

  # Correct the sites in the met_mat
  # for(c in which(collapse)){
  #   met_mat[,c-1] <- met_mat[,c-1] + met_mat[,c]
  # }
  # ## Fast version (using sparse matrix tricks and the fact that one read can only be on one strand)
  # if(shift == -1){
  #   for(co  in which(collapse)){
  #     # Shift all values in collapse strands one column to the left
  #     p = met_mat@p[(co-1):(co+1)] # P values that are needed
  #     met_mat@p[co] <- p[3] # Correct p for the now empty collapsed column
  #     met_mat@i[(p[1]+1):p[3]] <- sort(met_mat@i[(p[1]+1):p[3]]) # sort the i values of the two columns
  #   }
  #   met_mat <- met_mat[,!collapse, drop = F]
  # }else if(shift == 1){
  #   for(co  in which(collapse)){
  #     # Shift all values in collapse strands one column to the right
  #     p = met_mat@p[(co-1):(co+1)] # P values that are needed
  #     met_mat@p[co] <- p[1] # Correct p for the now empty collapsed column
  #     met_mat@i[(p[1]+1):p[3]] <- sort(met_mat@i[(p[1]+1):p[3]]) # sort the i values of the two columns
  #   }
  #   met_mat <- met_mat[,-(which(collapse) - 1), drop = F]
  # }
  co <- which(collapse)
  if(shift == -1){
    met_mat@p[co] <- met_mat@p[co - shift]
  }else if(shift == 1){
    met_mat@p[co+1] <- met_mat@p[co]
  }
  # Construct the matrix again to sort row indices correctly (construction function handles this)
  # see ?Matrix:::.validateCsparse
  met_mat <- sparseMatrix(i = met_mat@i,
                          p = met_mat@p,
                          x = met_mat@x,
                          dimnames = dimnames(met_mat),
                          dims = dim(met_mat),
                          index1 = F)
  # Matrix::.validateCsparse(met_mat, sort.if.needed = TRUE) # This handles sorting too and might be faster
  met_mat <- met_mat[,-co, drop = F]
  colnames(met_mat) <- as.character(start(met_locs))

  # Recalculate statistics
  res <- list(met_locs = met_locs, met_mat = met_mat)
  res <- calculateMetStats(res)
  res$met_locs <- addContext(res$met_locs, width = width(res$met_locs$context[1]))
  return(res)
}
#' Find the neighboring bases of a cytosine in a Granges object
#'
#' @param met_locs GRanges object indicating the mehtylated Cytosines
#' @param genome Bsgenome, Reference Genome of the alignment
#'
#' @return GRanges object with mehtylated Cytosines and added context metadata
#' @importFrom Biostrings getSeq
#' @import ggplot2
#' @export
#'
#' @examples
#' ml <- GRanges(seqnames = 'chr1',
#'               IRanges(start=c(3508273, 3508274, 3508275), width = 1),
#'               strand = c("-", "-", "+"))
#' addContext(ml)
#'
addContext <- function(met_locs, genome = Mmusculus, width = 5){
  met_locs$context=Biostrings::getSeq(genome,
                          resize(met_locs, width = width, fix='center'),
                          as.character=F)
  return(met_locs)
}

#' Filter reads using nonCpG/GpC sites as a refrerence
#'
#' @param met Methylation data object
#' @param conv_rate_cutoff the minimum conversion rate per read to be kept
#' @param plot boolean, whether to plot the distribution of methylation for CpG,GpC and nonCpGGpC sites
#' @param keep_noinfo boolean, whether to keep reads that have no non CpG/GpC sites.
#' (For those reads the conversion rate cannot be estimated)
#'
#' @return met subset for reads with sufficient conversion rate
#' @import ggplot2
#' @export
#'
filterByConvRate <- function(met, conv_rate_cutoff = 0.8, plot = TRUE, keep_noinfo = TRUE){
  met_locs <- met$met_locs
  if(!('cg_motif' %in% names(met$met_locs@elementMetadata))){
    met_locs <- addType(met_locs, motifs = c(cg_motif=DNAString("NCG")))
  }
  if(!('gc_motif' %in% names(met$met_locs@elementMetadata))){
    met_locs <- addType(met_locs, motifs = c(gc_motif=DNAString("GCN")))
  }
  cggc_mask <- met_locs$cg_motif | met_locs$gc_motif
  conv_rate <- rowMeans_drop0(met$met_mat[,!cggc_mask, drop=F]) - 1
  if(plot){
    plot_dt <- data.frame(non_GC_CG = conv_rate,
                          CG = rowMeans_drop0(met$met_mat[,met_locs$cg_motif, drop=F]) - 1,
                          GC = rowMeans_drop0(met$met_mat[,met_locs$gc_motif, drop=F]) - 1)
    plot_dt <- melt(plot_dt, id.vars=NULL, variable.name = "Methylation_Type", value.name = "Methylation_Rate")
    pl <- ggplot2::ggplot(plot_dt[!is.na(plot_dt$Methylation_Rate),], aes(x=Methylation_Rate, fill=Methylation_Type)) +
      geom_histogram(position = "dodge", binwidth = 0.05) +
      geom_vline(xintercept = 1 - conv_rate_cutoff, color = 'red') +
      theme_bw()
    print(pl)
  }
  if(keep_noinfo){
    conv_rate[is.na(conv_rate)] = 0
  }else{
    conv_rate[is.na(conv_rate)] = 1
  }
  ## Filter
  met$met_mat <- met$met_mat[conv_rate < (1 - conv_rate_cutoff),,drop=F]
  ## Recalculate statistics
  # Check if some sites have 0 reads now and remove those that do
  nnz = diff(met$met_mat@p)
  met$met_mat <- met$met_mat[,nnz > 0]
  met$met_locs <- met$met_locs[nnz > 0]
  # Recalculate the site methylation rates
  met <- calculateMetStats(met)
  return(met)
}

#' Add a type column to GRanges object of a methylation call by greping contexts
#'
#' @param met_locs GRanges object indicating the mehtylated Cytosines
#' @param motifs Vector of named DNAstrings that contain the motifs to extract
#'
#' @return GRanges object with added metadata columns indicating the type
#' @export
#'
#' @examples
addType <- function(met_locs,
                    motifs = c(acc_motif=DNAString("DGCHN"),
                               endmet_motif=DNAString("NWCGW"))){
  motif_l <- sapply(motifs, length)
  if(any(motif_l[1] != motif_l)){
    stop('All motifs must be the same length. You can add N to extend them')
  }
  if(motif_l[1]%%2 == 0){
    stop('Motifs must be the same width in both directions around C. You can add N characters to extend them')
  }
  if(!any(names(met_locs@elementMetadata) == 'context') | (length(met_locs$context[[1]]) < motif_l[1])){
    met_locs <- addContext(met_locs, width = motif_l[1])
    context <- met_locs$context
  } else if ((length(met_locs@elementMetadata$context[1]) > motif_l[1])){
    # Subset the context to the right width (ensures that there are no sporadic matches)
    context <- subseq(met_locs$context, width = motif_l[1], start = length(met_locs$context[[1]]) - motif_l[1])
  }else{
    context <- met_locs$context
  }
  if(any(names(motifs)=="")) stop("Some motifs are not named. Provide a name.")
  # Match all motifs to the C contexts
  mot_match <- sapply(motifs, vcountPattern, context, fixed=F)
  if(any(rowSums(mot_match)>1)) warning("Some sites match to multiple motifs")
  mot_match <- mot_match > 0
  met_locs@elementMetadata <- cbind(met_locs@elementMetadata, mot_match)
  return(met_locs)
}

subsetByMotifs <- function(met,motifs = c('acc_motif', 'endmet_motif'), discard_aditional_reads = TRUE){
  motmask <- apply(met$met_locs@elementMetadata[, motifs, drop = F], 1, any)
  met$met_locs <- met$met_locs[motmask]
  met$met_mat <- met$met_mat[, motmask]
  # Check if some reads can be thrown away
  if(discard_aditional_reads){
    nnz <- tabulate(met$met_mat@i + 1, nbins = nrow(met$met_mat))
    met$met_mat <- met$met_mat[nnz>0,,drop = F]
  }
  return(met)
}

subsetByGrange <- function(met, range){
  ov <- findOverlaps(met$met_locs, range)
  if(length(ov)==0) return(NULL)

  met$met_locs <- met$met_locs[queryHits(ov)]
  met$met_mat <- met$met_mat[, queryHits(ov), drop = F]
  met$met_mat <- met$met_mat[tabulate(met$met_mat@i + 1, nbins = nrow(met$met_mat))>0,,drop = F]
  met$subs_range <- range

  return(met)
}
#splitByMotifs <- function(met, motifs = c('acc_motif', 'endmet_motif'))
