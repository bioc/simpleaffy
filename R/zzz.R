.First.lib <-  function(lib,pkg,where) {
  library.dynam("simpleaffy",pkg,lib);
  require(affy,quietly=TRUE);
  require(methods,quietly=TRUE);
  where <- match(paste("package:", pkg, sep=""), search());
  .initClassesAndMethods(where);
  cacheMetaData(as.environment(where));
  cat("Welcome to 'simpleaffy' V 2.10                                \n");
  cat("      Produced by The Paterson Institute for Cancer Research \n");
  cat("      and funded by CANCER RESEARCH UK.                      \n");
  cat("      http://bioinformatics.picr.man.ac.uk/simpleaffy        \n");
  cat("      mailto: microarray@picr.man.ac.uk                      \n");
  .createQCEnvironment()
}

.initClassesAndMethods <- function(where) {

}
