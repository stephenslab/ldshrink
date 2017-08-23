libpath <- function() {
    cat(sprintf("%s/LDshrink/libs/LDshrink%s", installed.packages()["LDshrink", "LibPath"][1], .Platform$dynlib.ext))
}
incpath <- function() {
    cat(sprintf("-I%s/LDshrink/include/", installed.packages()["LDshrink", "LibPath"][1], .Platform$dynlib.ext))
}
