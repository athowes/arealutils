#' Develop the aghq models here!

load("data/mw.rda")
mw <- sf::st_as_sf(mw)

constant_aghq(mw)
iid_aghq(mw)
besag_aghq(mw)
bym2_aghq(mw)
fck_aghq(mw)
fik_aghq(mw)
system.time({ck_aghq(mw)})
system.time({ik_aghq(mw)})
