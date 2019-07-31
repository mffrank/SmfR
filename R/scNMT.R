
recoverStrand <- function(met_locs){
  met_locs <- addContext(met_locs, width = 1)
  minus_strand <- vcountPattern("C", met_locs$context, fixed = F) == 0
  strand(met_locs) <- '+'
  strand(met_locs)[minus_strand] <- "-"
  met_locs <- addContext(met_locs)
  return(met_locs)
}

collapseStrandsScNMT <- function(met, collapseCpG = TRUE, collapseGPC = TRUE){
  # This is a bit different than for smf since we can have info on the plus and minus strand in the same cell.
  # For convenience we'll take the info from the plus strand for now if both are present
  met_locs <- met$met_locs
  # Collapse CpG
  if(collapseCpG){
    if(!('cg_motif' %in% names(met$met_locs@elementMetadata))){
      met_locs <- addType(met_locs, motifs = c(cg_motif=DNAString("NCG")))
    }
    cg = met_locs$cg_motif
    met <- .collapseStrandsScNMT(met_locs, met$met_mat, filter = cg, shift = -1)
  }

  met_locs <- met$met_locs
  # Collapse GpC
  if(collapseGPC){
    if(!('gc_motif' %in% names(met$met_locs@elementMetadata))){
      met_locs <- addType(met_locs, motifs = c(gc_motif=DNAString("GCN")))
    }
    gc <- met_locs$gc_motif
    met <- .collapseStrandsScNMT(met_locs, met$met_mat, filter = gc, shift=1)
  }

  met <- calculateMetStats(met)
  met$met_locs <- addContext(met$met_locs, width = width(met$met_locs$context[1]))
  return(met)

}

.collapseStrandsScNMT <- function(met_locs, met_mat, filter, shift = 1){

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

  co <- which(collapse)
  # This is too slow
  # for(co in which(collapse)[1:100]){
  #   a = met_mat[,co]
  #   b = met_mat[,co+shift]
  #   res = a
  #   disagree <- a != b
  #   onlyB <- disagree & a == 0
  #   onlyA <- disagree & b == 0
  #   res[disagree] <- 0
  #   res[onlyA] <- a[onlyA]
  #   res[onlyB] <- b[onlyB]
  #   met_mat[,co] <- res
  # }
  # Collapse the columns with a trick
  met_mat@x[met_mat@x == 2] <- 3
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
  # Now the entries in the same i,j field got added (i.e. 1,3=4;0,1=1; etc)
  # We can assign the desired states
  met_mat@x[met_mat@x == 2] <- 1 # Comes from 1,1
  met_mat@x[met_mat@x == 3] <- 2 # Comes from 0,3
  met_mat[met_mat == 4] <- 0 # Comes from 1,3
  met_mat@x[met_mat@x == 6] <- 2 # Comes from 3,3

  met_mat <- met_mat[,-co, drop = F]
  colnames(met_mat) <- as.character(start(met_locs))

  res <- list(met_locs = met_locs, met_mat = met_mat)
  return(res)

}
