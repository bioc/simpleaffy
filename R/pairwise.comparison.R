"pairwise.comparison" <- function(x,group,members=NULL,logged=T,spots=NULL,a.order=NULL,b.order=NULL,method="unlogged") {
  if(is.null(members)) {
    pd <- unique(as.character(pData(x)[,group]))
    if(is.null(pd)) {
      stop(paste("Can't find a group called",group));
    }      
    if(length(pd) != 2)  {
      stop("There must be exactly two groups for a pairwise comparison. Please specify which groups you want to compare.");
    }
    members <- pd;
  }
  if(!is.null(spots)) {
    pmac <- detection.p.val(spots);
    results <- get.fold.change.and.t.test(x,group,members,logged=logged,return.exprs="samples",a.order=a.order,b.order=b.order,method=method);
    results <- c(results,pmac);
  }
  else {
    results <- get.fold.change.and.t.test(x,group,members,logged=logged,return.exprs="samples",a.order=a.order,b.order=b.order,method=method);
  }
  return(results); 
}
