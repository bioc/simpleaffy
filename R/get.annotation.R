get.annotation <- function(x,cdfname) {
  library(cdfname,character.only=T);
  symb  <- unlist(multiget(x,env=get(paste(cdfname,"SYMBOL",sep=""))),use.names=F);
  desc  <- unlist(multiget(x,env=get(paste(cdfname,"GENENAME",sep=""))),use.names=F);
  accno <- unlist(multiget(x,env=get(paste(cdfname,"ACCNUM",sep=""))),use.names=F);
  uni   <- unlist(multiget(x,env=get(paste(cdfname,"UNIGENE",sep=""))),use.names=F);
  acc.lnk <- paste("=HYPERLINK(\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=nucleotide&term=",accno,"\",\"",accno,"\")",sep="");
  uni.lnk <- paste("=HYPERLINK(\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=unigene&term=",uni,"&dopt=unigene\",\"",uni,"\")",sep="")	
  res <- cbind(symb,acc.lnk,uni.lnk,desc);
  colnames(res) <- c("gene name","accession","unigene","description");
  return(res);
}

write.annotation <- function(summary,file="results/annotation.table.xls") {
  write.table(summary,file=file,sep="\t",quote=F,col.names=NA)
}

results.summary <- function(results,cdfname) {
  res <- cbind(results$means,results$fc,sapply(results$fc,function(x) { if(x <0) { -1 * 2 ^ (-1 * x) } else { 2^x } }),results$tt,get.annotation(names(results$fc),cdfname));
  cns <- colnames(res);
  cns[3]<-"log2(fc)";
  cns[4]<-"fc";
  cns[5]<-"tt";
  colnames(res) <- cns;
  return(res);
}
