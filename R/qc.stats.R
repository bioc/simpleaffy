library("methods")

#holds the results of a pairwise comparison
setClass("QCStats",representation(scale.factors="numeric",target="numeric",percent.present="numeric",average.background="numeric",minimum.background="numeric",maximum.background="numeric",spikes="matrix",qc.probes="matrix"));

#accessor methods
setGeneric("sfs", function(object) standardGeneric("sfs"))
setMethod("sfs","QCStats",function(object) object@scale.factors)

setGeneric("target", function(object) standardGeneric("target"))
setMethod("target","QCStats",function(object) object@target)

setGeneric("percent.present", function(object) standardGeneric("percent.present"))
setMethod("percent.present","QCStats",function(object) object@percent.present)

setGeneric("avbg", function(object) standardGeneric("avbg"))
setMethod("avbg","QCStats",function(object) object@average.background)

setGeneric("minbg", function(object) standardGeneric("minbg"))
setMethod("minbg","QCStats",function(object) object@minimum.background)

setGeneric("maxbg", function(object) standardGeneric("maxbg"))
setMethod("maxbg","QCStats",function(object) object@maximum.background)

setGeneric("spikeInProbes", function(object) standardGeneric("spikeInProbes"))
setMethod("spikeInProbes","QCStats",function(object) object@spikes)

setGeneric("qcProbes", function(object) standardGeneric("qcProbes"))
setMethod("qcProbes","QCStats",function(object) object@qc.probes)


#for ratios
.namegrep3 <- function(stems,all) {
   sapply(stems,function(stem) {
     grep(paste(stem,"[_-]3.?_?.?_at$",sep=""),all,value=T)
   });
}

.namegrepM <- function(stems,all) {
   sapply(stems,function(stem) {
     grep(paste(stem,"[_-]M.?_?.?_at$",sep=""),all,value=T)
   });
}
.namegrep5 <- function(stems,all) {
   sapply(stems,function(stem) {
     grep(paste(stem,"[_-]5.?_?.?_at$",sep=""),all,value=T)
   });
}

.getRatios <- function(x) {
   vals <- x@qc.probes;
   unique.names <- colnames(vals)
   unique.names <- sub("[_-]5.?_?.?_at$","",unique.names,perl=T);
   unique.names <- sub("[_-]3.?_?.?_at$","",unique.names,perl=T);
   unique.names <- sub("[_-]M.?_?.?_at$","",unique.names,perl=T);
   unique.names <- unique(unique.names);
   p3 <- .namegrep3(unique.names,colnames(vals))
   p5 <- .namegrep5(unique.names,colnames(vals));
   pM <- .namegrepM(unique.names,colnames(vals));
   res1 <- rbind(c(),(vals[,p3] - vals[,p5]))

   colnames(res1) <- paste(unique.names,".3'/5'",sep="")

   res2 <- rbind(c(),(vals[,pM] - vals[,p5]))
   colnames(res2) <- paste(unique.names,".M /5'",sep="")
   r <- cbind(res1,res2)
   return(r)

}

setGeneric("ratios", function(object) standardGeneric("ratios"))
setMethod("ratios","QCStats",function(object) .getRatios(object))

#Create a new environment containing the probenames and parameters required by the qc functions
#These data are stored in the data directory of the package as tab delimited files

.createQCEnvironment <- function() {
  if(!exists(".qcEnv")) {
    .qcEnv <- new.env()
     data(alpha,envir =.qcEnv )
     data(spikes,envir =.qcEnv )
     data(qc.probes,envir =.qcEnv )
     assign(".qcEnv",.qcEnv,envir=globalenv())
  }
}


getTao <- function(name) {
  0.015;
}

getAlpha1 <- function(name) {
  as.character(get("alpha",envir=.qcEnv)[name,"alpha1"])
}

getAlpha2 <- function(name) {
  as.character(get("alpha",envir=.qcEnv)[name,"alpha2"])
}

getActin3 <- function(name) {
  as.character(get("qc.probes",envir=.qcEnv)[name,"actin3"])
}

getActinM <- function(name) {
  as.character(get("qc.probes",envir=.qcEnv)[name,"actinM"])
}

getActin5 <- function(name) {
  as.character(get("qc.probes",envir=.qcEnv)[name,"actin5"])
}

getGapdh3 <- function(name) {
  as.character(get("qc.probes",envir=.qcEnv)[name,"gapdh3"])
}

getGapdhM <- function(name) {
  as.character(get("qc.probes",envir=.qcEnv)[name,"gapdhM"])
}

getGapdh5 <- function(name) {
  as.character(get("qc.probes",envir=.qcEnv)[name,"gapdh5"])
}

getAllQCProbes <- function(name) {
  r <- as.matrix(get("qc.probes",envir=.qcEnv)[name,])
  n <- names(r)
  r <- r[!is.na(r)]
  names(r) <- n
  return(r)
}

getBioB <- function(name) {
  as.character(get("spikes",envir=.qcEnv)[name,"biob"])
}


getBioC <- function(name) {
  as.character(get("spikes",envir=.qcEnv)[name,"bioc"])
}


getBioD <- function(name) {
  as.character(get("spikes",envir=.qcEnv)[name,"biod"])
}


getCreX <- function(name) {
  as.character(get("spikes",envir=.qcEnv)[name,"crex"])
}


getAllSpikeProbes <- function(name) {
  r <- as.matrix(get("spikes",envir=.qcEnv)[name,])
  n <- names(r)
  r <- r[!is.na(r)]
  names(r) <- n
  return(r)
}


haveQCParams <- function(name) {
  name %in% rownames(get("alpha",envir=.qcEnv))
}



qc.affy <-function(unnormalised,normalised=NULL,tau=0.015,logged=TRUE,cdfn=cleancdfname(cdfName(unnormalised))) {

  if(is.null(normalised)) {getAllSpikeProbes("hgu133acdf")

    normalised <- call.exprs(unnormalised,"mas5");
  }

  if(!haveQCParams(cleancdfname(cdfName(unnormalised)))) {
	stop(paste("I'm sorry, I do not know about chip type:",cleancdfname(cdfName(unnormalised))))
  }

  x <- exprs(normalised);

  det <- detection.p.val(unnormalised,tau=tau,alpha1=getAlpha1(cdfn),alpha2=getAlpha2(cdfn));

  dpv<-apply(det$call,2,function(x) { 
            x[x!="P"] <- 0;
	    x[x=="P"] <- 1;
            x<-as.numeric(x);			       
            return(100 * sum(x)/length(x));
  });

  sfs    <- normalised@description@preprocessing$sfs;
  target <- normalised@description@preprocessing$tgt;

  if(!logged) { x <- log2(x); }

  bgsts<-.bg.stats(unnormalised)$zonebg

  meanbg <- apply(bgsts,1,mean);

  minbg  <- apply(bgsts,1,min);

  maxbg  <- apply(bgsts,1,max);

  stdvbg <- sqrt(apply(bgsts,1,var));




  #get the probenames for the QC probes for this chip
  qc.probenames <- getAllQCProbes(cdfn);

  qc.probe.vals <- rbind(c(),(sapply(qc.probenames, function(y) {x[y,]})))
  rownames(qc.probe.vals) <- colnames(x);

  spike.probenames <- getAllSpikeProbes(cdfn);
  spike.vals <- rbind(c(),(sapply(spike.probenames, function(y) {x[y,]})))

  rownames(spike.vals) <- colnames(x);
  return(new("QCStats",scale.factors=sfs,target=target,percent.present=dpv,average.background=meanbg,minimum.background=minbg,maximum.background=maxbg,
              spikes=spike.vals,qc.probes=qc.probe.vals));
}


setGeneric("qc", function(unnormalised,...) standardGeneric("qc"))
setMethod("qc","AffyBatch",function(unnormalised,...) qc.affy(unnormalised,...)) 

setGeneric("getQCParams", function(x) standardGeneric("getQCParams") )

setMethod("getQCParams","AffyBatch",function(x) {
  n <- cleancdfname(cdfName(x))

  res <- list(getGapdh3(n),getGapdhM(n),getGapdh5(n),getActin3(n),getActinM(n),getActin5(n),getBioB(n),getBioC(n),getBioD(n),getCreX(n),getAlpha1(n),getAlpha2(n),getTao(n))
  names(res) <- c("Gapdh3","GapdhM","Gapdh5","Actin3","ActinM","Actin5","BioB","BioC","BioD","CreX","Alpha1","Alpha2","Tao")
  return(res)
})


.bg.stats <- function(unnormalised, grid=c(4,4)) {
pms         <- unlist(pmindex(unnormalised))
mms         <- unlist(mmindex(unnormalised))
all         <- c(pms,mms)
intensities <- exprs(unnormalised)
rws <- nrow(unnormalised)
cls <- ncol(unnormalised)
zonebg <- c();
zonesd <- c();
for(no in 1:length(unnormalised)){
  this.array <- intensities[,no];
  result <- .C("bgmas",as.integer(as.vector(all)),as.integer(length(all)),
       as.double(as.vector(this.array)),as.integer(length(this.array)),
       as.integer(rws),
       as.integer(cls),
       as.integer(grid[1]),as.integer(grid[2]),
       zonebg=double(grid[1] * grid[2]),
       zonesd=double(grid[1] * grid[2]),corrected=double(length(this.array)),PACKAGE="simpleaffy");
  zonesd <- rbind(zonesd, result$zonesd);
  zonebg <- rbind(zonebg, result$zonebg);
  }
  colnames(zonesd) <- paste("zone",1:16,sep=".");
  colnames(zonebg) <- paste("zone",1:16,sep=".");
  rownames(zonesd) <- sampleNames(unnormalised);
  rownames(zonebg) <- sampleNames(unnormalised);
  return(list(zonebg=zonebg,zonesd=zonesd))
}


