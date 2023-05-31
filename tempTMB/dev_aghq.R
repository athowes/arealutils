#' Develop the aghq models here!
mw <- sf::st_as_sf(mw)

constant_aghq(mw, its = 1000)
iid_aghq(mw, its = 1000)
besag_aghq(mw, its = 1000)
bym2_aghq(mw, its = 1000)
fck_aghq(mw, its = 1000)
