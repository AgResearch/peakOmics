deIsotoping=function(peakData, peakId, ppm.tol=20, rt.tol=0.2, cor.tol=0.9){
  if(ncol(peakData)!=nrow(peakId))
    stop("Provide data in such format: peakData is in sample x peak, and 
         peakId is a data.frame with required columns 'mz' and 'rt'")
  d.cor=cor(peakData, use="complete.obs")
  d.cor[lower.tri(d.cor, diag=TRUE) ]=0

  MZ=peakId$mz; RT=peakId$rt
  idxM=which(d.cor>cor.tol, arr.ind=TRUE)
  #make a matrix
  p1.MZ=MZ[idxM[,1]]; p2.MZ=MZ[idxM[,2]]
  p1.RT=RT[idxM[,1]]; p2.RT=RT[idxM[,2]]
  MZ.diff=abs(p1.MZ-p2.MZ)
  RT.diff=abs(p1.RT-p2.RT)
  #add ppm to isoData
  ppm=MZ.diff/p1.MZ*1e6         #by p2.MZ
  isoData=cbind(idxM, p1.MZ, p2.MZ, p1.RT, p2.RT, MZ.diff, ppm, RT.diff)
  #*** same peak
  #if mz is judged by ppm, the same peak should be treated differently
  samePeak.idx=which(ppm < ppm.tol & RT.diff<rt.tol)
  cat("---- Groups of same peaks ---\n")
  print(isoData[samePeak.idx, ])

  #the same peak -- take strong one
  isoData[samePeak.idx,]
  idx1=isoData[samePeak.idx, 1:2]
  wk.idx=rep(NA, nrow(idx1))
  for(i in 1:nrow(idx1)){
    criterion1=median(peakData[,idx1[i,1]])>median(peakData[,idx1[i,2]])
    wk.idx[i]=ifelse(criterion1, idx1[i,1], idx1[i, 2])
  }
  length(wk.idx)                     #remove it, peak idx
  #keep the peaknames to be removed
  #samePeak.wk=colnames(d1)[wk.idx]

  #***isotopic peaks [M+1]
  isoPeak.idx=which(abs(MZ.diff-1.007276)/p1.MZ*1e6 < ppm.tol & RT.diff<rt.tol)
  isoData[isoPeak.idx,]
  idx2=isoData[isoPeak.idx, 1:2]
  #need to keep both monomass idx and weaker idx
  mono.idx=iso.wk.idx=rep(NA, nrow(idx2))
  #pair.idx=NULL;
  for(i in 1:nrow(idx2)){
    m1=mean(peakData[,idx2[i,1]]);
    m2=mean(peakData[,idx2[i,2]]);
    if(m1>m2 & m1/m2 > 2){
      mono.idx[i]=idx2[i,1]  #keep stronger
      iso.wk.idx[i]=idx2[i,2]   #keep the weaker
      #pair.idx=rbind(pair.idx, idx2[i,c(1,2)])  #and order it
    }else if(m1<m2 & m2/m1 > 2){
      mono.idx[i]=idx2[i,2]
      iso.wk.idx[i]=idx2[i,1]
      #pair.idx=rbind(pair.idx, idx2[i,c(2,1)])
    }else{
      cat("Check: possible isotope pairs but its ratio is < the criterion: \n")
      cat(colnames(peakData)[idx2[i,1]], "and", colnames(peakData)[idx2[i,2]], "\n")
      cat("The ratio of the pair", max(c(m1,m2))/min(c(m1,m2)), "\n\n")
      #!!! should retain both? Here simply take the first one
      mono.idx[i]=idx2[i,1]
      iso.wk.idx[i]=idx2[i,2] #
      #pair.idx=rbind(pair.idx, idx2[i,])
    }
  }
  mono.idx                                #keep it
  iso.wk.idx                              #remove it

  #give the identified pairs, search M+2 in d1 matrix
  M1=colnames(peakData)[iso.wk.idx];
  searchList=peakId[iso.wk.idx,]
  M2.idx=rep(0, length(iso.wk.idx))
  for(i in 1:length(M1)){
    mz1=searchList$mz[i]
    rt1=searchList$rt[i]
    mz1.diff=peakId$mz-mz1  #search those > mz1
    rt1.diff=peakId$rt-rt1
    #M2 - M1 > (1 +/- delta)
    delta=mz1*ppm.tol*1e-6
    #w=which( mz1.diff>(1.007276-delta) & mz1.diff<(1.007276+delta))
    idx=which( mz1.diff>(1.007276-delta) & mz1.diff<(1.007276+delta) & abs(rt1.diff)<rt.tol)
    #cat("to search pairs of: ", M1[i], "\n")
    #check peak intensity
    M1.k=match(M1[i], colnames(peakData))
    M1.int=mean(peakData[,M1.k])

    if(length(idx)==1){
      M2.int=mean(peakData[,idx])
      if((M1.int/M2.int) >2)
        M2.idx[i]=idx
    }else if(length(idx)>1){
      min.k=which.min(colMeans(peakData[,idx]))
      #print(min.k); print(w)
      M2.int=mean(peakData[,idx[min.k]])
      if((M1.int/M2.int)>2) #[M+1]/[M+2]
        M2.idx[i]=idx[min.k]
      #print(idx[min.k])
    }else{
      M2.idx[i]=0 #not found
    }
  }
  tmp.idx=cbind(mono.idx,  iso.wk.idx, M2.idx)
  tmp.idx
  #barplot(colMeans(d1[,tmp.idx[176,]]))

  #!! update "peakId" with isopeaks detected added to help annotation
  isoCluster=rep("", nrow(peakId))
  Mplus2.name=rep("", length(M2.idx))
  for(i in 1:length(M2.idx)){
    if(M2.idx[i]!=0) Mplus2.name[i]=peakId$pkname[M2.idx[i]]
  }
  Mplus2.name
  isoNames=paste(peakId$pkname[iso.wk.idx], Mplus2.name, sep="/")
  isoNames
  isoCluster[mono.idx]=isoNames
  peakId=data.frame(peakId, isoPeaks=isoCluster)

  #!!remove both the sample peak and weak isotopic
  remove.idx=c(wk.idx, iso.wk.idx, M2.idx[M2.idx!=0])  #OK with duplicated idx
  length(unique(remove.idx))

  d2=peakData[,-remove.idx]                      #d2:  dIso data
  dim(d2)
  peakId2=peakId[-remove.idx,]             #860 peaks retained
  cat("Checking -- ", ncol(d2)==nrow(peakId2), "\n")
  cat("Checking peak names -- ", all(peakId2$pkname==colnames(d2)), "\n")

  list(dat=d2, peakId=peakId2)
}#>>>
