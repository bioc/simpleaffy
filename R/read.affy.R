"read.affy" <-
function(covdesc="covdesc",...) {
  samples <- read.phenoData(covdesc);
  files.to.read <- rownames(pData(samples));
  eset <- ReadAffy(filenames=files.to.read,...);
  newPhenoData <- cbind(pData(eset),pData(samples)[rownames(pData(eset)),]);
  colnames(newPhenoData) <- c(colnames(pData(eset)),colnames(pData(samples)));
  tmp <- as.list(colnames(newPhenoData));
  names(tmp) <- colnames(newPhenoData);
  newPhenoData <- new("phenoData",pData=newPhenoData,varLabels=tmp)

  eset@phenoData <- newPhenoData;	
  return(eset);
}
