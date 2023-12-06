#' Develop the other TMB models here!

load("data/mw.rda")
mw <- sf::st_as_sf(mw)

constant_tmb(mw, its = 1000)
iid_tmb(mw, its = 1000)
besag_tmb(mw, its = 1000)
bym2_tmb(mw, its = 1000)
fck_tmb(mw, its = 1000)
fik_tmb(mw, its = 1000)
system.time({ck_tmb(mw, its = 1000)})
system.time({ik_tmb(mw, its = 1000)})

#' Test that the cross-validation enabling code works
cv_test <- iid_tmb(mw, its = 1000, ii_mis = 0)
summary(cv_test)
