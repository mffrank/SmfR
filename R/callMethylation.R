# Functions to call methylation events from alignment files using QuasR

#' Call methylation events on all cytosines
#' @param proj qProject object
#' @param range GRanges object with ONE!!!! range
#' @param sample name of the Genome alignment

getCMethMatrix<-function(proj,range,sample){
  Cs=qMeth(proj[sample], query=range,mode="allC",reportLevel="alignment")
  # use data.table to get a 1,0 matrix of methylation profiles
  all.cids=unique(Cs[[sample]]$Cid) # get all possible C locations
  # make the data.table object
  dt=data.table(meth=Cs[[sample]]$meth ,aid=Cs[[sample]]$aid ,cid=Cs[[sample]]$Cid)
  dtm=dcast(dt, aid~cid, value.var="meth")
  ronames=dtm$aid
  dtm[,aid:=NULL] # remove unwanted row
  CpGm=as.matrix(dtm)
  rownames(CpGm)=ronames
  return(CpGm)
}
