
detection.p.val <- function(x, tau = 0.015,calls=TRUE,alpha1=getAlpha1(cleancdfname(cdfName(x))),alpha2=getAlpha2(cleancdfname(cdfName(x))),ignore.saturated=TRUE) {
  if(alpha1 < 0)      {stop("alpha1 must be  > 0 "); }
  if(alpha1 > alpha2) {stop("alpha2 must be  > alpha1 "); }
  if(alpha2 > 1)      {stop("alpha2 must be  <1 "); }

  cat("Getting probe level data...\n");
  pms <-as.matrix(pm(x));
  mms <-as.matrix(mm(x));

  # Saturation:
  # shouldn't be a problem with new scanners or those that have had an engineer visit
  if(ignore.saturated) { sat <- 46000; }
  else { sat <- -1; }
  
  pns <- probeNames(x);
  unique.pns <- unique(pns);
  cat("Computing p-values\n");
  p<-sapply(1:length(pms[1,]),function(x) { 
    .C("DetectionPValue",as.double(pms[,x]),as.double(mms[,x]),as.character(pns),as.integer(length(mms[,x])),
	as.double(tau),as.double(sat),dpval=double(length(unique.pns)),length(unique.pns),PACKAGE="simpleaffy")$dpval;
  });
  rownames(p) <- unique.pns;
  colnames(p) <- paste(sampleNames(x),".detection.p.val",sep="");
  if(!calls) { 
    l <- list(detection.p.values=p); 
  }
  else       {
    cat("Doing PMA Calls\n");

    calls <- sapply(p,function(y) { if(y < alpha1) { return("P") } else { if(y < alpha2) { return("M") } else { return("A") }}});
    calls <- matrix(calls,nrow=nrow(p),ncol=ncol(p));
    colnames(calls) <- paste(sampleNames(x),".present",sep="");
    rownames(calls) <- rownames(p)
    l<- list(pval=p,call=calls);
    return(l); 
  }

}     
