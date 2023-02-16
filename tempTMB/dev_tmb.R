#' Develop the other TMB models here!
compile("besag.cpp")
dyn.load(dynlib("besag"))
