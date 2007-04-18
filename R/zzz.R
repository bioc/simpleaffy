.First.lib <-  function(libname, pkgname) {
    require("methods", quietly=TRUE)
    library.dynam("simpleaffy", pkgname, libname)
}
.onLoad <- .First.lib
