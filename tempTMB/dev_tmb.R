#' Develop the other TMB models here!

load("data/mw.rda")
mw <- sf::st_as_sf(mw)

constant_tmb(mw, its = 1000)
iid_tmb(mw, its = 1000)
besag_tmb(mw, its = 1000)
bym2_tmb(mw, its = 1000)
fck_tmb(mw, its = 1000)
fik_tmb(mw, its = 1000)

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
