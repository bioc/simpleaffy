qc.stats <-function(unnormalised,normalised=NULL,logged=T,tau=0.015,alpha1=0.04,alpha2=0.06) {
  if(is.null(normalised)) {
    normalised <- call.exprs(unnormalised,"mas5");
  }
  x <- exprs(normalised);
  det <- detection.p.val(unnormalised,tau=tau,alpha1=alpha1,alpha2=alpha2);
  dpv<-apply(det$call,2,function(x) {
            x[x!="P"] <- 0;
	    x[x=="P"] <- 1;
            x<-as.numeric(x);
            return(100 * sum(x)/length(x));
  });
  sfs <- normalised@description@preprocessing$sfs;
  sfs2<-log2(sfs);
  sfs3<-sapply(sfs2,function(a){ if( a <0 ) { -2^-a} else { 2^a}})
  if(!logged) { x <- log2(x); }
  n <- cleancdfname(cdfName(unnormalised));
  if(n %in% c("hgu133acdf", "hgu133atagcdf", "hgu133bcdf", "hgu133plus2cdf", "hgu133a_2cdf", "hgu95acdf",
              "hgu95av2cdf", "hgu95bcdf", "hgu95ccdf", "hgu95dcdf", "hgu95ecdf")) {
    rats <- 2^cbind((x["AFFX-HUMGAPDH/M33197_3_at",] - x["AFFX-HUMGAPDH/M33197_5_at",]),
                    (x["AFFX-HUMGAPDH/M33197_3_at",] - x["AFFX-HUMGAPDH/M33197_M_at",]),
                    (x["AFFX-HSAC07/X00351_3_at",]  - x["AFFX-HSAC07/X00351_5_at",]),
                    (x["AFFX-HSAC07/X00351_3_at",]  - x["AFFX-HSAC07/X00351_M_at",]));

  }
  else {
    if (n %in% c("mgu74acdf", "mgu74av2cdf", "mgu74bcdf", "mgu74bv2cdf", "mgu74ccdf", "mgu74cv2cdf", "moe430acdf", "moe430bcdf", "moe430_2cdf", "moe430a_2cdf")) {
      rats <- 2^cbind((x["AFFX-GapdhMur/M32599_3_at",] - x["AFFX-GapdhMur/M32599_5_at",]),
                      (x["AFFX-GapdhMur/M32599_3_at",] - x["AFFX-GapdhMur/M32599_M_at",]),
                      (x["AFFX-b-ActinMur/M12481_3_at",]  - x["AFFX-b-ActinMur/M12481_5_at",]),
                     (x["AFFX-b-ActinMur/M12481_3_at",]  - x["AFFX-b-ActinMur/M12481_M_at",]));
    }
    else {
      if (n %in% c("rae230acdf")) {
        rats <- 2^cbind((x["AFFX_Rat_GAPDH_3_at",] - x["AFFX_Rat_GAPDH_5_at",]),
                        (x["AFFX_Rat_GAPDH_3_at",] - x["AFFX_Rat_GAPDH_M_at",]),
                        (x["AFFX_Rat_beta-actin_3_at",]  - x["AFFX_Rat_beta-actin_5_at",]),
                       (x["AFFX_Rat_beta-actin_3_at",]  - x["AFFX_Rat_beta-actin_5_at",]));
      }
      else { stop(paste("Sorry - I'm afraid I don't know about the spike probes on '",n,"' arrays.")) }
    }
  }


  bgsts<-bg.stats(unnormalised)$zonebg
  backgrounds<-apply(bgsts,1,mean);
  backgrounds<-cbind(backgrounds,apply(bgsts,1,min));
  backgrounds<-cbind(backgrounds,apply(bgsts,1,max));
  backgrounds<-cbind(backgrounds,sqrt(apply(bgsts,1,var)));
  if(n %in% c("hgu133acdf", "hgu133atagcdf", "hgu133bcdf","hgu133a_2cdf", "hgu133plus2cdf", "moe430acdf", "moe430bcdf", "moe430_2cdf", "moe430a_2cdf" )) {
    spikes <- cbind(det$call["AFFX-r2-Ec-bioB-3_at",],
                    det$call["AFFX-r2-Ec-bioC-3_at",],
                    det$call["AFFX-r2-Ec-bioD-3_at",],
                    det$call["AFFX-r2-P1-cre-3_at",]);
  }
  else  {
     if (n %in% c("mgu74acdf", "mgu74av2cdf", "mgu74bcdf", "mgu74bv2cdf", "mgu74ccdf", "mgu74cv2cdf","hgu95acdf",
                  "hgu95av2cdf", "hgu95bcdf", "hgu95ccdf", "hgu95dcdf", "hgu95ecdf","rae230acdf")) {
       spikes <- cbind(det$call["AFFX-BioB-3_at",],
                       det$call["AFFX-BioC-3_at",],
                       det$call["AFFX-BioDn-3_at",],
                      det$call["AFFX-CreX-3_at",]);
     }
  }

  spikes[spikes != "P"] <- 0;
  spikes[spikes == "P"] <- 1;
  spikes <- matrix(as.numeric(spikes),nrow=nrow(spikes))
  rats <- cbind(rats,as.vector(sfs),as.vector(sfs2),as.vector(sfs3),normalised@description@preprocessing$tgt,dpv,backgrounds,spikes);
  colnames(rats) <- c("GAPDH.3'/5'","GAPDH.3'/M","beta.actin.3'/5'","beta.actin.3'/M",
                      "s.f.","log2(s.f.)","fc(s.f.)","tgt",
                      "%present",
                      "avbg","minbg","maxbg","stdevbg",
                      "bioB","bioC","bioD","creX");
  return(rats);
}


bg.stats <- function(unnormalised, grid=c(4,4)) {
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
       zonesd=double(grid[1] *
       grid[2]),corrected=double(length(this.array)), PACKAGE="simpleaffy");
  zonesd <- rbind(zonesd, result$zonesd);
  zonebg <- rbind(zonebg, result$zonebg);
  }
  colnames(zonesd) <- paste("zone",1:16,sep=".");
  colnames(zonebg) <- paste("zone",1:16,sep=".");
  rownames(zonesd) <- sampleNames(unnormalised);
  rownames(zonebg) <- sampleNames(unnormalised);
  return(list(zonebg=zonebg,zonesd=zonesd))
}

plot.qc.stats<-function(x,unnormalised,fc.line.col="black",sf.ok.region="light blue",chip.label.col="black",sf.thresh = 3.0,gdh.thresh = 3.0,ba.thresh = 3.0,present.thresh=10,bg.thresh=20,label=NULL,...) {
  n <- length(rownames(x));
  sfs <- x[,"fc(s.f.)"]
  meansf <- mean(sfs)

  dpv <- x[,"%present"]
  dpv <- (round(100*dpv))/100;
  abg <- round(x[,"avbg"])

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
  plot(x1[1],y1[1],pch=".",xlim=range(-12,12),ylim=range(-12,12),xaxt="n",yaxt="n",xlab="",ylab="",col=sf.ok.region,...);
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
  for(i in 1:n) {
    if(d1 > sf.thresh) { col = "red" } else {col="blue"}
     x1 <- sfs[i] * cos(i * arc);
     y1 <- sfs[i] * sin(i * arc);
     x2 <- 6 * cos(i * arc);
     y2 <- 6 * sin(i * arc);
     lines(c(x2,x1),c(y2,y1),col=col);
     points(x1,y1,col=col,pch=20);
     text(x1,y1,col=chip.label.col,label=label[i],adj=0.2 * c(cos(i * arc),sin(i * arc)));
     gdh <- log2(x[i,1]);
     ba  <- log2(x[i,3]);
     x2 <- (6 + gdh) * cos(i * arc);
     y2 <- (6 + gdh) * sin(i * arc);
     if(gdh > gdh.thresh) { col = "red" } else {col="blue"}
     points(x2,y2,pch=1,col=col);
     x2 <- (6 + ba) * cos(i * arc);
     y2 <- (6 + ba) * sin(i * arc);
     if(ba > ba.thresh) { col = "red" } else {col="blue"}
     points(x2,y2,pch=2,col=col);

     if(d2 > present.thresh) { col = "red" } else {col="blue"}
     x2 <- (9 * cos(i * arc));
     y2 <- (9 * sin(i * arc));
     text(x2,y2,label=paste(dpv[i],"%",sep=""),col=col);

     if(d3 > bg.thresh) { col = "red" } else {col="blue"}
     x2 <- (10 * cos(i * arc));
     y2 <- (10 * sin(i * arc));
     text(x2,y2,label=abg[i],col=col);

     if(x[i,"bioB"]!="1") {
       x2 <- (11 * cos(i * arc));
       y2 <- (11 * sin(i * arc));
       text(x2,y2,label="bioB",col="red");
     }
  }
  legend(-10,10,pch=1:2,c("GAPDH","Beta Actin"))
}
