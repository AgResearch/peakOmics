normalizeRunorder <-
function(dmat, run.seq, log2=FALSE){
  if(log2) warnings("fitting log-linear")
  dmat=as.matrix(dmat)
  if(is.null(colnames(dmat)))
    stop("provide a matrix with colnames")
  if(nrow(dmat)!=length(run.seq))
    stop("data matrix should be in the dimension of sample x ions");
  dn=matrix(0, nrow(dmat), ncol(dmat))
  colnames(dn)=colnames(dmat)
  for(j in 1:ncol(dmat)){
    ion.d=dmat[,j]
    if(any(is.na(ion.d))){
       ion.d[which(is.na(ion.d))]=min(ion.d, na.rm=T)
    }
    if(log2) ion.d=log2(ion.d) #0?, for when log2 applied
    lm1=lm(ion.d ~ run.seq)
    bsln=fitted(lm1)
    margin=max(bsln)-bsln
    baseoff=ion.d+margin
    dn[,j]=baseoff
  }
  dn
}
