#' Develop the other TMB models here!
library(TMB)

mw <- sf::st_as_sf(mw)

constant_tmb(mw, its = 1000)
iid_tmb(mw, its = 1000)
besag_tmb(mw, its = 1000)

#' BYM2
compile("bym2.cpp")
dyn.load(dynlib("bym2"))

#' Develop package function

#' Test works using packaged function
bym2_tmb(mw, its = 1000)

#' FCK
compile("fck.cpp")
dyn.load(dynlib("fck"))

#' Develop package function

#' Test works using packaged function

#' FIK
compile("fik.cpp")
dyn.load(dynlib("fik"))

#' Develop package function

#' Test works using packaged function

#' CK
compile("ck.cpp")
dyn.load(dynlib("ck"))

#' Develop package function

#' Test works using packaged function

#' IK
compile("ik.cpp")
dyn.load(dynlib("ik"))

#' Develop package function

#' Test works using packaged function
