#' Develop the tmbstan models here!
mw <- sf::st_as_sf(mw)

constant_tmbstan(mw)
iid_tmbstan(mw)
besag_tmbstan(mw)
bym2_tmbstan(mw)
fck_tmbstan(mw)
fik_tmbstan(mw)
# ck_tmbstan(mw)
# fk_tmbstan(mw)