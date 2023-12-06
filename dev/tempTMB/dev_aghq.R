#' Develop the aghq models here!

load("data/mw.rda")
mw <- sf::st_as_sf(mw)

system.time({constant_aghq(mw)})
system.time({iid_aghq(mw)})
system.time({besag_aghq(mw)})
system.time({bym2_aghq(mw)})
system.time({fck_aghq(mw)})
system.time({fik_aghq(mw)})
system.time({ck_aghq(mw)})
system.time({ik_aghq(mw)})

#' Test that the cross-validation enabling code works
cv_test <- iid_aghq(mw, ii_mis = 0)
summary(cv_test)
