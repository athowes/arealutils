#' Develop the tmbstan models here!

load("data/mw.rda")
mw <- sf::st_as_sf(mw)

constant_tmbstan(mw)
iid_tmbstan(mw)
besag_tmbstan(mw)
bym2_tmbstan(mw)
fck_tmbstan(mw)
fik_tmbstan(mw)
system.time({ck_tmbstan(mw)})
system.time({ik_tmbstan(mw)})

#' Test that the cross-validation enabling code works
cv_test <- iid_tmbstan(mw, ii_mis = 0)
summary(cv_test)
cv_test
