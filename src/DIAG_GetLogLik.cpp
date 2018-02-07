#include <RcppArmadillo.h>
#include "DIAG_diagnostics.h"

//Function that computes the log-likelihood for spBDwDM model----------------------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec GetLogLik(Rcpp::List DatObj_List, Rcpp::List Para_List,
                       int NKeep) {

  //Convert Rcpp::Lists to C++ structs
  datobjDIAG DatObj = ConvertDatObjDIAG(DatObj_List);
  paraDIAG Para = ConvertParaDIAG(Para_List);

  //Set data object
  int FamilyInd = DatObj.FamilyInd;

  //Compute log-likelihood
  arma::colvec LogLik(NKeep);
  if (FamilyInd == 0) LogLik = NormalLogLik(DatObj, Para, NKeep);
  if (FamilyInd == 1) LogLik = ProbitLogLik(DatObj, Para, NKeep);
  if (FamilyInd == 2) LogLik = TobitLogLik(DatObj, Para, NKeep);

  //Return log-likelihood
  return LogLik;

}



//Function that computes the log-likelihood for STBDwDM model for the mean parameters------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double GetLogLikMean(Rcpp::List DatObj_List, Rcpp::List Para_List) {

  //Convert Rcpp::Lists to C++ structs
  datobjDIAG DatObj = ConvertDatObjDIAG(DatObj_List);
  paraDIAG Para = ConvertParaDIAG(Para_List);

  //Set data object
  int FamilyInd = DatObj.FamilyInd;

  //Compute log-likelihood
  double LogLikMean;
  if (FamilyInd == 0) LogLikMean = NormalLogLikMean(DatObj, Para);
  if (FamilyInd == 1) LogLikMean = ProbitLogLikMean(DatObj, Para);
  if (FamilyInd == 2) LogLikMean = TobitLogLikMean(DatObj, Para);

  //Return log-likelihood
  return LogLikMean;

}



//Function that computes the log-likelihood for spBDwDM model----------------------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec GetLogLik_lmc(Rcpp::List DatObj_List, Rcpp::List Para_List,
                       int NKeep) {

  //Convert Rcpp::Lists to C++ structs
  datobjDIAG DatObj = ConvertDatObjDIAG(DatObj_List);
  paraDIAG Para = ConvertParaDIAG(Para_List);

  //Set data object
  int FamilyInd = DatObj.FamilyInd;

  //Compute log-likelihood
  arma::colvec LogLik(NKeep);
  if (FamilyInd == 0) LogLik = NormalLogLik_lmc(DatObj, Para, NKeep);
  if (FamilyInd == 1) LogLik = ProbitLogLik_lmc(DatObj, Para, NKeep);
  if (FamilyInd == 2) LogLik = TobitLogLik_lmc(DatObj, Para, NKeep);

  //Return log-likelihood
  return LogLik;

}



//Function that computes the log-likelihood for STBDwDM model for the mean parameters------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double GetLogLikMean_lmc(Rcpp::List DatObj_List, Rcpp::List Para_List) {

  //Convert Rcpp::Lists to C++ structs
  datobjDIAG DatObj = ConvertDatObjDIAG(DatObj_List);
  paraDIAG Para = ConvertParaDIAG(Para_List);

  //Set data object
  int FamilyInd = DatObj.FamilyInd;

  //Compute log-likelihood
  double LogLikMean;
  if (FamilyInd == 0) LogLikMean = NormalLogLikMean_lmc(DatObj, Para);
  if (FamilyInd == 1) LogLikMean = ProbitLogLikMean_lmc(DatObj, Para);
  if (FamilyInd == 2) LogLikMean = TobitLogLikMean_lmc(DatObj, Para);

  //Return log-likelihood
  return LogLikMean;

}



//Function that computes the log-likelihood for spBDwDM model----------------------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec GetLogLik_novar(Rcpp::List DatObj_List, Rcpp::List Para_List,
                           int NKeep) {

  //Convert Rcpp::Lists to C++ structs
  datobjDIAG DatObj = ConvertDatObjDIAG(DatObj_List);
  paraDIAG Para = ConvertParaDIAG(Para_List);

  //Set data object
  int FamilyInd = DatObj.FamilyInd;

  //Compute log-likelihood
  arma::colvec LogLik(NKeep);
  if (FamilyInd == 0) LogLik = NormalLogLik_novar(DatObj, Para, NKeep);
  if (FamilyInd == 1) LogLik = ProbitLogLik_novar(DatObj, Para, NKeep);
  if (FamilyInd == 2) LogLik = TobitLogLik_novar(DatObj, Para, NKeep);

  //Return log-likelihood
  return LogLik;

}



//Function that computes the log-likelihood for STBDwDM model for the mean parameters------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double GetLogLikMean_novar(Rcpp::List DatObj_List, Rcpp::List Para_List) {

  //Convert Rcpp::Lists to C++ structs
  datobjDIAG DatObj = ConvertDatObjDIAG(DatObj_List);
  paraDIAG Para = ConvertParaDIAG(Para_List);

  //Set data object
  int FamilyInd = DatObj.FamilyInd;

  //Compute log-likelihood
  double LogLikMean;
  if (FamilyInd == 0) LogLikMean = NormalLogLikMean_novar(DatObj, Para);
  if (FamilyInd == 1) LogLikMean = ProbitLogLikMean_novar(DatObj, Para);
  if (FamilyInd == 2) LogLikMean = TobitLogLikMean_novar(DatObj, Para);

  //Return log-likelihood
  return LogLikMean;

}
