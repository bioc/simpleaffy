"pairwise.comparison" <- function(x,group,members,spots=NULL) {
  if(!is.null(spots)) {
    pmac <- detection.p.val(spots);
    results <- get.fold.change.and.t.test(x,group,members,return.exprs="samples");
    results <- c(results,pmac);
  }
  else {
    results <- get.fold.change.and.t.test(x,group,members,return.exprs="samples");
  }
  return(results); 
}
