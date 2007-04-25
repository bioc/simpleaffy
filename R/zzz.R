.onLoad <-  function(libname, pkgname) {
    require("methods", quietly=TRUE)
    library.dynam("simpleaffy", pkgname, libname)
    .initClassesAndMethods()
}

.initClassesAndMethods <- function() {
   .createQCEnvironment()
}
