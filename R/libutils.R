inlineCxxPlugin <- function(...) {
    ismacos <- Sys.info()[["sysname"]] == "Darwin"
    openmpflag <- if (ismacos) "" else "$(SHLIB_OPENMP_CFLAGS)"
    plugin <- Rcpp::Rcpp.plugin.maker(
                  include.before = "#include <LDshrink.h>",
                  libs           = paste(
                      openmpflag,
                      RcppParallel::RcppParallelLibs(),
                      "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)"
                  ),
                  package        = "LDshrink"
              )
    settings <- plugin()
    settings$env$PKG_CPPFLAGS <- paste("-I../inst/include", openmpflag)
    if (!ismacos) settings$env$USE_CXX11 <- "yes"
    settings
}
