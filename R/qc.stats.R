library("methods")

#holds the results of a pairwise comparison
setClass("QCStats",representation(scale.factors="numeric",target="numeric",percent.present="numeric",average.background="numeric",minimum.background="numeric",maximum.background="numeric",spikes="matrix",qc.probes="matrix",bioBCalls="character"));

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
  get("alpha",envir=.qcEnv)[name,"alpha1"]
}

getAlpha2 <- function(name) {
  get("alpha",envir=.qcEnv)[name,"alpha2"]
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

  bb <- getBioB(cdfn)

  if(!is.na(bb)) {
    biobcalls <- det$call[bb,]
  }
  else {
     biobcalls <- NULL
  }

  return(new("QCStats",scale.factors=sfs,target=target,percent.present=dpv,average.background=meanbg,minimum.background=minbg,maximum.background=maxbg,
              spikes=spike.vals,qc.probes=qc.probe.vals,bioBCalls=biobcalls));
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



.plot.qc.stats2<-function(x,fc.line.col,sf.ok.region,chip.label.col,sf.thresh,gdh.thresh,ba.thresh,present.thresh,bg.thresh,label,main,usemid,...) {

  sfs    <- log2(sfs(x))

  n      <- length(sfs)

  meansf <- mean(sfs)

  dpv <- percent.present(x)
  dpv <- (round(100*dpv))/100;

  abg <- avbg(x)
  abg <- (round(100*abg))/100;
	
  sfs <- sfs + 6.0;
  if(is.null(label)) { label <- 1:n }
  col=c("red","green");
  d1 <- 0.0;
  d2 <- 0.0;
  d3 <- 0.0;

  for(i in 1:n) {
    for(j in 1:n) { 
      d1 <- max(abs(sfs[i] - sfs[j]),d1);
      d2 <- max(abs(dpv[i] - dpv[j]),d1);
      d3 <- max(abs(abg[i] - abg[j]),d3);
    }
  }
  i<-1:101;
  arc <- 2 * pi / 100;
  x1 <- (7.5 +  meansf) * cos(i * arc);
  y1 <- (7.5 +  meansf) * sin(i * arc);
  x2 <- (4.5 +  meansf) * cos(i * arc);
  y2 <- (4.5 +  meansf) * sin(i * arc);

  plot(x1[1],y1[1],pch=".",xlim=range(-12,12),ylim=range(-12,12),xaxt="n",yaxt="n",xlab="",ylab="",col=sf.ok.region,main=main,...);

  polygon(c(x1,rev(x2),x1[1]),c(y1,rev(y2),y1[1]),col=sf.ok.region,border=sf.ok.region);

  x1 <- 3 * cos(i * arc);
  y1 <- 3 * sin(i * arc);
  points(x1,y1,pch=".",col=fc.line.col);
  x1 <- 6 * cos(i * arc);
  y1 <- 6 * sin(i * arc);
  points(x1,y1,pch=".",col=fc.line.col);
  x1 <- 9 * cos(i * arc);
  y1 <- 9 * sin(i * arc);
  points(x1,y1,pch=".",col=fc.line.col);
  text(0,3,"-3",col=fc.line.col)
  text(0,6,"0",col=fc.line.col)
  text(0,9,"+3",col=fc.line.col)
  arc <- 2 * pi / n;
  rats <- ratios(x);
  if(!usemid) {
    gdh <- rats[,2];
    ba  <- rats[,1];
  }
  else {  
    gdh <- rats[,4];
    ba  <- rats[,3];
  }
  bb  <- x@bioBCalls

  for(i in 1:n) {
    if(d1 > sf.thresh) { col = "red" } else {col="blue"}
     x1 <- sfs[i] * cos(i * arc);
     y1 <- sfs[i] * sin(i * arc);
     x2 <- 6 * cos(i * arc);
     y2 <- 6 * sin(i * arc);
     lines(c(x2,x1),c(y2,y1),col=col);
     points(x1,y1,col=col,pch=20);
     text(x1,y1,col=chip.label.col,label=label[i],adj=0.2 * c(cos(i * arc),sin(i * arc)));
     x2 <- (6 + gdh[i]) * cos(i * arc);
     y2 <- (6 + gdh[i]) * sin(i * arc);
     if(gdh[i] > gdh.thresh) { col = "red" } else {col="blue"}	
     points(x2,y2,pch=1,col=col);
     x2 <- (6 + ba[i]) * cos(i * arc);
     y2 <- (6 + ba[i]) * sin(i * arc);
     if(ba[i] > ba.thresh) { col = "red" } else {col="blue"}	
     points(x2,y2,pch=2,col=col);

     if(d2 > present.thresh) { col = "red" } else {col="blue"}
     x2 <- (9 * cos(i * arc));
     y2 <- (9 * sin(i * arc));
     text(x2,y2,label=paste(dpv[i],"%",sep=""),col=col);

     if(d3 > bg.thresh) { col = "red" } else {col="blue"}
     x2 <- (11 * cos(i * arc));
     y2 <- (11 * sin(i * arc));
     text(x2,y2,label=abg[i],col=col);

     if(bb[i]!="P") {
       x2 <- (12 * cos(i * arc));
       y2 <- (12 * sin(i * arc));
       text(x2,y2,label="bioB",col="red");
     }
  }
  if(!usemid) {
    legend(-11,12,pch=2:1,colnames(rats)[1:2])
  }
  else {
    legend(-11,12,pch=2:1,colnames(rats)[3:4])
  }
}

plot.qc.stats<-function(x,fc.line.col="black",sf.ok.region="light blue",chip.label.col="black",sf.thresh = 3.0,gdh.thresh = 1.25,ba.thresh = 3.0,present.thresh=10,bg.thresh=20,label=NULL,main="QC Stats",usemid=F,spread=c(-8,8),type="l",...) {
  if(type=="c") { 
    .plot.qc.stats2(x,fc.line.col,sf.ok.region,chip.label.col,sf.thresh,gdh.thresh,ba.thresh,present.thresh,bg.thresh,label,main,usemid,...) 
    return()
  }
  old.par <- par()
  par(mai=c(0,0,0,0))
  sfs    <- log2(sfs(x))

  n      <- length(sfs)

  meansf <- mean(sfs)

  dpv <- percent.present(x)
  dpv <- (round(100*dpv))/100;

  abg <- avbg(x)
  abg <- (round(100*abg))/100;
	
  if(is.null(label)) { label <- names(maxbg(x)) }
  d1 <- 0.0;
  d2 <- 0.0;
  d3 <- 0.0;

  for(i in 1:n) {
    for(j in 1:n) { 
      d1 <- max(abs(sfs[i] - sfs[j]),d1);
      d2 <- max(abs(dpv[i] - dpv[j]),d1);
      d3 <- max(abs(abg[i] - abg[j]),d3);
    }
  }

  # set up plotting area - a column for array names next to a column for the QC

  m <- matrix(c(4,2,1,3) ,nrow=2,ncol=2)
  layout(m,c(1,2),c(0.1,1))
  # the title
  if(is.null(main)) { main="" }
  plot(0,0,xlim=range(0,1),ylim=range(0,1),type="n",yaxs="i",xaxt="n",yaxt="n",bty="n")
  text(0.5,0.5,labels=main,adj=0,cex=2)

  # write out the array names

  cex <- 1.0
  while((maxwidth <- max(sapply(label,function(a) { strwidth(a,units="f",cex)} ))) >= 1) cex <- cex/maxwidth

  plot(0,0,xlim=range(0,1),ylim=range(-1,n),type="n",yaxs="i",xaxt="n",yaxt="n",bty="n")
  text(1,(1:n)-0.5,labels=label,adj=1)

  plot(0,0,xlim=spread,ylim=c(-1,n),type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",bty="n")

  x1 <- (1.5 +  meansf)
  y1 <- 0
  x2 <- (-1.5 +  meansf)
  y2 <- n

  polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=sf.ok.region,border=sf.ok.region);
  lines(c(0,0),c(0,n),lty=1,col=fc.line.col)
  lines(c(-1,-1),c(0,n),lty=2,col="grey")
  lines(c(-2,-2),c(0,n),lty=2,col="grey")
  lines(c(-3,-3),c(0,n),lty=2,col=fc.line.col)
  lines(c(1,1),c(0,n),lty=2,col="grey")
  lines(c(2,2),c(0,n),lty=2,col="grey")
  lines(c(3,3),c(0,n),lty=2,col=fc.line.col)
  text(3,-1,"3",pos=3,col=fc.line.col)
  text(2,-1,"2",pos=3,col=fc.line.col)
  text(1,-1,"1",pos=3,col=fc.line.col)
  text(-3,-1,"-3",pos=3,col=fc.line.col)
  text(-2,-1,"-2",pos=3,col=fc.line.col)
  text(-1,-1,"-1",pos=3,col=fc.line.col)
  text(0,-1,"0",pos=3,col=fc.line.col)

  rats <- ratios(x);
  if(!usemid) {
    gdh <- rats[,2];
    ba  <- rats[,1];
  }
  else {
    gdh <- rats[,4];
    ba  <- rats[,3];
  }

  bb  <- x@bioBCalls

  for(i in 1:n) {
    x1<-spread[1]
    x2<-spread[2]
    y1<-i-1;
    y2<-i-1;
    lines(c(x1,x2),c(y1,y2),lty=2,col="light grey")
    if(d1 > sf.thresh) { col = "red" } else {col="blue"}
     x1 <- sfs[i]
     y1 <- i-0.25
     lines(c(0,x1),c(y1,y1),col=col);

     points(x1,y1,col=col,pch=20);
     x2 <- gdh[i]
     y2 <- i-0.5;
     if(gdh[i] > gdh.thresh) { col = "red" } else {col="blue"}	
     points(x2,y2,pch=1,col=col);

     x2 <- ba[i];
     y2 <- i-0.5;
     if(ba[i] > ba.thresh) { col = "red" } else {col="blue"}	
     points(x2,y2,pch=2,col=col);

     if(d2 > present.thresh) { col = "red" } else {col="blue"}
     x2 <- spread[1]
     y2 <- i-0.25
     dpvs<-paste(dpv[i],"%",sep="")
     text(x2,y2,label=dpvs,col=col,pos=4);
     if(d3 > bg.thresh) { col = "red" } else {col="blue"}
     x2 <- spread[1]
     y2 <- i-0.75
     text(x2,y2,label=abg[i],col=col,pos=4);
     if(bb[i]!="P") {
       x2 <- (12 * cos(i * arc));
       y2 <- (12 * sin(i * arc));
       text(0,i-1,label="bioB",col="red");
     }

  }
  plot(0,0,xlim=range(0,1),ylim=range(0,1),type="n",yaxs="i",xaxt="n",yaxt="n",bty="n")
  if(!usemid) {
    points(0.25,0.25,pch=1)
    text(0.3,0.25,colnames(rats)[2],pos=4)
    points(0.25,0.5,pch=2)
    text(0.3,0.5,colnames(rats)[1],pos=4)
  }
  else {
    points(0.25,0.25,pch=1)
    text(0.3,0.25,colnames(rats)[4],pos=4)
    points(0.25,0.5,pch=2)
    text(0.3,0.5,colnames(rats)[3],pos=4)
  }
  
  ow <- options("warn")$warn
  options(warn=-1)
  par(old.par)
  options(warn=ow)
}


setMethod("plot","QCStats",function(x,y) plot.qc.stats(x,...))
