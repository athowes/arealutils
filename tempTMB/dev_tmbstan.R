#' Develop the tmbstan models here!
mw <- sf::st_as_sf(mw)

#' Constant
compile("constant.cpp")
dyn.load(dynlib("constant"))

#' Develop package function

#' Test works using packaged function
constant_tmbstan(mw)
