#' Develop the other TMB models here!
library(TMB)

#' IID
compile("iid.cpp")
dyn.load(dynlib("iid"))
#' Fit model here...

#' Besag
compile("besag.cpp")
dyn.load(dynlib("besag"))
#' Fit model here...

#' BYM2
compile("bym2.cpp")
dyn.load(dynlib("bym2"))
#' Fit model here...

#' FCK
compile("fck.cpp")
dyn.load(dynlib("fck"))
#' Fit model here...

#' FIK
compile("fik.cpp")
dyn.load(dynlib("fik"))
#' Fit model here...

#' CK
compile("ck.cpp")
dyn.load(dynlib("ck"))
#' Fit model here...

#' IK
compile("ik.cpp")
dyn.load(dynlib("ik"))
#' Fit model here...
