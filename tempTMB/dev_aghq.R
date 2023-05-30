#' Develop the aghq models here!
mw <- sf::st_as_sf(mw)

constant_aghq(mw)
iid_aghq(mw)
besag_aghq(mw)
