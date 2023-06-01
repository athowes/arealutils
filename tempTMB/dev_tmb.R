#' Develop the other TMB models here!

load("data/mw.rda")
mw <- sf::st_as_sf(mw)

constant_tmb(mw, its = 1000)
iid_tmb(mw, its = 1000)
besag_tmb(mw, its = 1000)
bym2_tmb(mw, its = 1000)
fck_tmb(mw, its = 1000)
fik_tmb(mw, its = 1000)
# ck_tmb(mw, its = 1000)
# ik_tmb(mw, its = 1000)

#' CK
TMB::compile("ck.cpp")
dyn.load(TMB::dynlib("ck"))

#' IK
TMB::compile("ik.cpp")
dyn.load(TMB::dynlib("ik"))
