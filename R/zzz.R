.onLoad <-  function(libname, pkgname) {
    require("methods", quietly=TRUE)
    require("utils", quietly=TRUE)
    require("affy", quietly=TRUE)
    require("genefilter", quietly=TRUE)
    require("Biobase", quietly=TRUE)
    library.dynam("simpleaffy", pkgname, libname)
    .initClassesAndMethods()
}

.initClassesAndMethods <- function() {
   .createQCEnvironment()
}
