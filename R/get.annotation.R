get.annotation <- function(x,cdfname) {
  library(cdfname,character.only=T);
  symb  <- unlist(multiget(x,env=get(paste(cdfname,"SYMBOL",sep=""))),use.names=F);
  desc  <- unlist(multiget(x,env=get(paste(cdfname,"GENENAME",sep=""))),use.names=F);
  accno <- unlist(multiget(x,env=get(paste(cdfname,"ACCNUM",sep=""))),use.names=F);
  uni   <- unlist(multiget(x,env=get(paste(cdfname,"UNIGENE",sep=""))),use.names=F);
  acc.lnk <- paste("=HYPERLINK(\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=nucleotide&term=",accno,"\",\"",accno,"\")",sep="");
  uni.lnk <- paste("=HYPERLINK(\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=unigene&term=",uni,"&dopt=unigene\",\"",uni,"\")",sep="")	
  return(cbind(symb,acc.lnk,uni.lnk,desc));
}

write.annotation <- function(summary,file="results/annotation.table.xls") {
  write.table(x,file=file,sep="\t",quote=F)
}

results.summary<- function(results,cdfname) {
  cbind(x$means,x$fc,x$tt,get.annotation(names(x$fc),cdfname))
}
