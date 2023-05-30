// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_arealutils_TMBExports
#include <TMB.hpp>
#include "besag.hpp"
#include "bym2.hpp"
#include "constant.hpp"
#include "iid.hpp"
#include "mvn_covariance.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "besag") {
    return besag(this);
  } else if(model == "bym2") {
    return bym2(this);
  } else if(model == "constant") {
    return constant(this);
  } else if(model == "iid") {
    return iid(this);
  } else if(model == "mvn_covariance") {
    return mvn_covariance(this);
  } else {
    Rf_error("Unknown model.");
  }
  return 0;
}
