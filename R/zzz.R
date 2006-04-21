.First.lib <-  function(lib,pkg,where) {
  library.dynam("simpleaffy",pkg,lib);
  require(affy,quietly=TRUE);
  require(methods,quietly=TRUE);
  where <- match(paste("package:", pkg, sep=""), search());
  .initClassesAndMethods(where);
  cacheMetaData(as.environment(where));
  .createQCEnvironment()
}

.initClassesAndMethods <- function(where) {

}
