#' Develop the aghq models here!

load("data/mw.rda")
mw <- sf::st_as_sf(mw)

constant_aghq(mw, its = 1000)
iid_aghq(mw, its = 1000)
besag_aghq(mw, its = 1000)
bym2_aghq(mw, its = 1000)
fck_aghq(mw, its = 1000)
fik_aghq(mw, its = 1000)
system.time({ck_aghq(mw, its = 1000)})
system.time({ik_aghq(mw, its = 1000)})
