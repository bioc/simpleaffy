library("methods")
library("affy")
# holds the results of a pairwise comparison
setClass("PairComp",representation(means="matrix",fc="numeric",tt="numeric",calls="matrix",group="character",members="character"))

#accessor methods
setGeneric("means", function(object) standardGeneric("means"))
setMethod("means","PairComp",function(object) object@means)

setGeneric("fc", function(object) standardGeneric("fc"))
setMethod("fc","PairComp",function(object) object@fc)

setGeneric("tt", function(object) standardGeneric("tt"))
setMethod("tt","PairComp",function(object) object@tt)

setGeneric("calls", function(object) standardGeneric("calls"))
setMethod("calls","PairComp",function(object) object@calls)

setGeneric("group", function(object) standardGeneric("group"))
setMethod("group","PairComp",function(object) object@group)

setGeneric("members", function(object) standardGeneric("members"))
setMethod("members","PairComp",function(object) object@members)


##subseting. can only happen by gene

setMethod("[", "PairComp", function(x,i,j,...,drop=FALSE) {
  if(nrow(calls(x))>0) {calls <- calls(x)[i,,...,drop=FALSE];}
  else { calls <- matrix(); }
  y <- new ("PairComp",means=means(x)[i,,...,drop=FALSE],fc=fc(x)[i,...,drop=FALSE],tt=tt(x)[i,...,drop=FALSE],calls=calls,group=group(x),members=members(x))
  return(y)
})

setReplaceMethod("[", "PairComp", function(x, i,j,...,value) {
  stop("operation not supported");
})



"get.array.subset.exprset" <-
function(x,group,members) {
  pd <- pData(x);
  grp <- pd[,colnames(pd) == group];
  return(x[,is.element(grp,members)]);
}

"get.array.subset.affybatch" <-
function(x,group,members) {
  pd <- pData(x);
  grp <- pd[,colnames(pd) == group];
  return(x[is.element(grp,members),]);
}

setGeneric("get.array.subset", function(x,group,members) standardGeneric("get.array.subset"))
setMethod("get.array.subset","AffyBatch",get.array.subset.affybatch);
setMethod("get.array.subset","exprSet",get.array.subset.exprset);

"get.fold.change.and.t.test" <- function(x,group,members,logged = TRUE, a.order=NULL,b.order=NULL,method=c("unlogged","logged","median")) {
  
  a.samples <- exprs(get.array.subset(x,group,members[1]));
  b.samples <- exprs(get.array.subset(x,group,members[2]));
 
  pw <- FALSE;

  if(!is.null(a.order)) { 
    a.samples <- a.samples[,a.order];
    if(!is.null(b.order)) { 
      b.samples <- b.samples[,b.order]; 
      pw <- TRUE;
 
    }
    else {
     stop("Both a.order and b.order must be specified for a paired t-test");
    }
  }
  method <- match.arg(method)
  m <- switch(method,
              logged   = 2,
              unlogged = 1,
              median   = 3);

  a.samples.array <- as.double(t(a.samples));
  b.samples.array <- as.double(t(b.samples));

  nacol <- as.integer(length(a.samples[1,]));
  ngene <- as.integer(length(a.samples[,1]));
  nbcol <- as.integer(length(b.samples[1,]));

  if(class(logged) != "logical") stop("Parameter 'logged' should be TRUE or FALSE")
  if((nacol == 1) | (nbcol == 1))  warning("There was only one sample in one (or both) of your sample groups. Not computing t-tests - instead, returning 0.0 for p-scores...");

  c.res <- .C("FCMandTT",a.samples.array,b.samples.array,nacol,nbcol,ngene,as.logical(logged),pw,as.integer(m),ma = double(ngene),mb = double(ngene),fc = double(ngene),tt = double(ngene),PACKAGE="simpleaffy")
  means <- cbind(c.res$ma,c.res$mb);
  colnames(means) <- members;
  fc <- c.res$fc;
  tt <- c.res$tt;
  names(fc)       <- rownames(a.samples);
  rownames(means) <- rownames(a.samples);
  names(tt)       <- rownames(a.samples);
  
  return(new("PairComp",fc=fc,tt=tt,means=means,group=group,members=members))
}





"pairwise.comparison" <- function(x,group,members=NULL,spots=NULL,a.order=NULL,b.order=NULL,method="unlogged",logged=TRUE) {
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
    results <- get.fold.change.and.t.test(x,group,members,logged=logged,a.order=a.order,b.order=b.order,method=method);
    calls(results) <- pmac;
  }
  else {
    results <- get.fold.change.and.t.test(x,group,members,logged=logged,a.order=a.order,b.order=b.order,method=method);
  }
  return(results); 
}


pairwise.filter <- function(object,x,min.exp=log2(100),min.exp.no=0,min.present.no=3,fc=1.0,tt=0.001) {

  if(class(object) != "PairComp") { stop("Can only filter an object of class 'PairComp'"); }
  if(class(x) != "exprSet") { stop("Can only filter using class 'exprSet' for parameter 'x'"); }
  pass.fc              <- (abs(fc(object)) > fc);
  pass.tt              <- (tt(object) < tt );

  samples <- exprs(get.array.subset(x,group(object),members(object)));

  no.chips             <- length(colnames(samples));

  intensity.thresh     <- array(sapply(samples[,1:no.chips],function(x) { if(x > min.exp) { 1 } else { 0 } } ),dim=dim(samples));

  min.intensity.thresh <- rowSums(intensity.thresh)

  
  pass.intensity <- (min.intensity.thresh >= min.exp.no);

  if(nrow(calls(object))>0) {
    present.thresh      <- array(sapply(calls(object)[,1:no.chips],function(x) { if(x == "P") { 1 } else { 0 } } ),dim=dim(calls(object)));
    min.present.thresh  <- rowSums(present.thresh);
    pass.present   <- (min.present.thresh >= min.present.no);

    return(object[(pass.fc & pass.tt & pass.intensity & pass.present),]);
  }
  else {
    return(object[(pass.fc & pass.tt & pass.intensity),]);
  }
}
