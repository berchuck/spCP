#include <RcppArmadillo.h>
#include "MCMC_spCP.h"

//Initiate burn-in progress bar--------------------------------------------------------------------------
void BeginBurnInProgress(mcmcobj McmcObj, bool Interactive) {

  //Set MCMC object
  int BarLength = McmcObj.BarLength;

  //Initialize burn-in bar
  if (Interactive) {
    Rcpp::Rcout << std::fixed << "Burn-in progress:  |";
    for (int i = 0; i < BarLength - 1; i++) Rcpp::Rcout << std::fixed << " ";
    Rcpp::Rcout << std::fixed <<  "|" << std::fixed;
  }
  if (!Interactive) {
    Rcpp::Rcout << std::fixed << "Burn-in progress:  0%..  ";
  }

}



//Function to pilot adapt tuning parameter--------------------------------------------------------------
double PilotAdaptFunc(double TuningParameter, double AcceptancePct) {

  //Adjust tuning parameter using scaling based on size of acceptance rate
  if (AcceptancePct >= 0.90) TuningParameter *= 1.3;
  if ( (AcceptancePct >= 0.75 ) & (AcceptancePct < 0.90 ) ) TuningParameter *= 1.2;
  if ( (AcceptancePct >= 0.45 ) & (AcceptancePct < 0.75 ) ) TuningParameter *= 1.1;
  if ( (AcceptancePct <= 0.25 ) & (AcceptancePct > 0.15 ) ) TuningParameter *= 0.9;
  if ( (AcceptancePct <= 0.15 ) & (AcceptancePct > 0.10 ) ) TuningParameter *= 0.8;
  if (AcceptancePct <= 0.10) TuningParameter *= 0.7;
  return TuningParameter;

}



//Function for implementing pilot adaptation in MCMC sampler--------------------------------------------
metrobj PilotAdaptation(datobj DatObj, metrobj MetrObj, mcmcobj McmcObj) {

  //Set data objects
  int M = DatObj.M;

  //Set Metropolis objects
  arma::vec MetropLambda0Vec = MetrObj.MetropLambda0Vec;
  arma::vec AcceptanceLambda0Vec = MetrObj.AcceptanceLambda0Vec;
  arma::vec MetropLambda1Vec = MetrObj.MetropLambda1Vec;
  arma::vec AcceptanceLambda1Vec = MetrObj.AcceptanceLambda1Vec;
  arma::vec MetropEtaVec = MetrObj.MetropEtaVec;
  arma::vec AcceptanceEtaVec = MetrObj.AcceptanceEtaVec;
  double MetropAlpha = MetrObj.MetropAlpha;
  double AcceptanceAlpha = MetrObj.AcceptanceAlpha;

  //Set MCMC objects
  int PilotAdaptDenominator = McmcObj.PilotAdaptDenominator;

  //Get acceptance percentages
  arma::vec PctLambda0Vec = AcceptanceLambda0Vec / double(PilotAdaptDenominator);
  arma::vec PctLambda1Vec = AcceptanceLambda1Vec / double(PilotAdaptDenominator);
  arma::vec PctEtaVec = AcceptanceEtaVec / double(PilotAdaptDenominator);
  double PctAlpha = AcceptanceAlpha / double(PilotAdaptDenominator);

  //Update Tuning Parameter
  for (int i = 0; i < M; i++) MetropLambda0Vec(i) = PilotAdaptFunc(MetropLambda0Vec(i), PctLambda0Vec(i));
  for (int i = 0; i < M; i++) MetropLambda1Vec(i) = PilotAdaptFunc(MetropLambda1Vec(i), PctLambda1Vec(i));
  for (int i = 0; i < M; i++) MetropEtaVec(i) = PilotAdaptFunc(MetropEtaVec(i), PctEtaVec(i));
  MetropAlpha = PilotAdaptFunc(MetropAlpha, PctAlpha);
  MetrObj.MetropLambda0Vec = MetropLambda0Vec;
  MetrObj.MetropLambda1Vec = MetropLambda1Vec;
  MetrObj.MetropEtaVec = MetropEtaVec;
  MetrObj.MetropAlpha = MetropAlpha;

  //Zero the acceptance counters
  AcceptanceLambda0Vec.zeros();
  AcceptanceLambda1Vec.zeros();
  AcceptanceEtaVec.zeros();
  AcceptanceAlpha = 0;
  MetrObj.AcceptanceLambda0Vec = AcceptanceLambda0Vec;
  MetrObj.AcceptanceLambda1Vec = AcceptanceLambda1Vec;
  MetrObj.AcceptanceEtaVec = AcceptanceEtaVec;
  MetrObj.AcceptanceAlpha = AcceptanceAlpha;
  return MetrObj;

}



//Function for implementing pilot adaptation in MCMC sampler--------------------------------------------
metrobj_lmc PilotAdaptation_lmc(datobj_lmc DatObj, metrobj_lmc MetrObj, mcmcobj McmcObj) {

  //Set data objects
  int M = DatObj.M;

  //Set Metropolis objects
  arma::vec MetropBeta0Vec = MetrObj.MetropBeta0Vec;
  arma::vec AcceptanceBeta0Vec = MetrObj.AcceptanceBeta0Vec;
  arma::vec MetropBeta1Vec = MetrObj.MetropBeta1Vec;
  arma::vec AcceptanceBeta1Vec = MetrObj.AcceptanceBeta1Vec;
  arma::vec MetropLambda0Vec = MetrObj.MetropLambda0Vec;
  arma::vec AcceptanceLambda0Vec = MetrObj.AcceptanceLambda0Vec;
  arma::vec MetropLambda1Vec = MetrObj.MetropLambda1Vec;
  arma::vec AcceptanceLambda1Vec = MetrObj.AcceptanceLambda1Vec;
  arma::vec MetropEtaVec = MetrObj.MetropEtaVec;
  arma::vec AcceptanceEtaVec = MetrObj.AcceptanceEtaVec;
  arma::vec MetropSigma = MetrObj.MetropSigma;
  arma::vec AcceptanceSigma = MetrObj.AcceptanceSigma;
  arma::vec MetropAlpha = MetrObj.MetropAlpha;
  arma::vec AcceptanceAlpha = MetrObj.AcceptanceAlpha;

  //Set MCMC objects
  int PilotAdaptDenominator = McmcObj.PilotAdaptDenominator;

  //Get acceptance percentages
  arma::vec PctBeta0Vec = AcceptanceBeta0Vec / double(PilotAdaptDenominator);
  arma::vec PctBeta1Vec = AcceptanceBeta1Vec / double(PilotAdaptDenominator);
  arma::vec PctLambda0Vec = AcceptanceLambda0Vec / double(PilotAdaptDenominator);
  arma::vec PctLambda1Vec = AcceptanceLambda1Vec / double(PilotAdaptDenominator);
  arma::vec PctEtaVec = AcceptanceEtaVec / double(PilotAdaptDenominator);
  arma::vec PctSigma = AcceptanceSigma / double(PilotAdaptDenominator);
  arma::vec PctAlpha = AcceptanceAlpha / double(PilotAdaptDenominator);

  //Update Tuning Parameter
  for (int i = 0; i < M; i++) MetropBeta0Vec(i) = PilotAdaptFunc(MetropBeta0Vec(i), PctBeta0Vec(i));
  for (int i = 0; i < M; i++) MetropBeta1Vec(i) = PilotAdaptFunc(MetropBeta1Vec(i), PctBeta1Vec(i));
  for (int i = 0; i < M; i++) MetropLambda0Vec(i) = PilotAdaptFunc(MetropLambda0Vec(i), PctLambda0Vec(i));
  for (int i = 0; i < M; i++) MetropLambda1Vec(i) = PilotAdaptFunc(MetropLambda1Vec(i), PctLambda1Vec(i));
  for (int i = 0; i < M; i++) MetropEtaVec(i) = PilotAdaptFunc(MetropEtaVec(i), PctEtaVec(i));
  for (int i = 0; i < 15; i++) MetropSigma(i) = PilotAdaptFunc(MetropSigma(i), PctSigma(i));
  for (int i = 0; i < 5; i++) MetropAlpha(i) = PilotAdaptFunc(MetropAlpha(i), PctAlpha(i));
  MetrObj.MetropBeta0Vec = MetropBeta0Vec;
  MetrObj.MetropBeta1Vec = MetropBeta1Vec;
  MetrObj.MetropLambda0Vec = MetropLambda0Vec;
  MetrObj.MetropLambda1Vec = MetropLambda1Vec;
  MetrObj.MetropEtaVec = MetropEtaVec;
  MetrObj.MetropSigma = MetropSigma;
  MetrObj.MetropAlpha = MetropAlpha;

  //Zero the acceptance counters
  AcceptanceBeta0Vec.zeros();
  AcceptanceBeta1Vec.zeros();
  AcceptanceLambda0Vec.zeros();
  AcceptanceLambda1Vec.zeros();
  AcceptanceEtaVec.zeros();
  AcceptanceSigma.zeros();
  AcceptanceAlpha.zeros();
  MetrObj.AcceptanceBeta0Vec = AcceptanceBeta0Vec;
  MetrObj.AcceptanceBeta1Vec = AcceptanceBeta1Vec;
  MetrObj.AcceptanceLambda0Vec = AcceptanceLambda0Vec;
  MetrObj.AcceptanceLambda1Vec = AcceptanceLambda1Vec;
  MetrObj.AcceptanceEtaVec = AcceptanceEtaVec;
  MetrObj.AcceptanceSigma = AcceptanceSigma;
  MetrObj.AcceptanceAlpha = AcceptanceAlpha;
  return MetrObj;

}



//Function for implementing pilot adaptation in MCMC sampler--------------------------------------------
metrobj_novar PilotAdaptation_novar(datobj_novar DatObj, metrobj_novar MetrObj, mcmcobj McmcObj) {

  //Set data objects
  int M = DatObj.M;

  //Set Metropolis objects
  arma::vec MetropLambdaVec = MetrObj.MetropLambdaVec;
  arma::vec AcceptanceLambdaVec = MetrObj.AcceptanceLambdaVec;
  arma::vec MetropEtaVec = MetrObj.MetropEtaVec;
  arma::vec AcceptanceEtaVec = MetrObj.AcceptanceEtaVec;
  double MetropAlpha = MetrObj.MetropAlpha;
  double AcceptanceAlpha = MetrObj.AcceptanceAlpha;

  //Set MCMC objects
  int PilotAdaptDenominator = McmcObj.PilotAdaptDenominator;

  //Get acceptance percentages
  arma::vec PctLambdaVec = AcceptanceLambdaVec / double(PilotAdaptDenominator);
  arma::vec PctEtaVec = AcceptanceEtaVec / double(PilotAdaptDenominator);
  double PctAlpha = AcceptanceAlpha / double(PilotAdaptDenominator);

  //Update Tuning Parameter
  for (int i = 0; i < M; i++) MetropLambdaVec(i) = PilotAdaptFunc(MetropLambdaVec(i), PctLambdaVec(i));
  for (int i = 0; i < M; i++) MetropEtaVec(i) = PilotAdaptFunc(MetropEtaVec(i), PctEtaVec(i));
  MetropAlpha = PilotAdaptFunc(MetropAlpha, PctAlpha);
  MetrObj.MetropLambdaVec = MetropLambdaVec;
  MetrObj.MetropEtaVec = MetropEtaVec;
  MetrObj.MetropAlpha = MetropAlpha;

  //Zero the acceptance counters
  AcceptanceLambdaVec.zeros();
  AcceptanceEtaVec.zeros();
  AcceptanceAlpha = 0;
  MetrObj.AcceptanceLambdaVec = AcceptanceLambdaVec;
  MetrObj.AcceptanceEtaVec = AcceptanceEtaVec;
  MetrObj.AcceptanceAlpha = AcceptanceAlpha;
  return MetrObj;

}



//Output Metropolis object for summary-------------------------------------------------------------------
Rcpp::List OutputMetrObj(metrobj MetrObj) {

  Rcpp::List Out = Rcpp::List::create(Rcpp::Named("AcceptanceLambda0Vec") = MetrObj.AcceptanceLambda0Vec,
                                      Rcpp::Named("MetropLambda0Vec") = MetrObj.MetropLambda0Vec,
                                      Rcpp::Named("AcceptanceLambda1Vec") = MetrObj.AcceptanceLambda1Vec,
                                      Rcpp::Named("MetropLambda1Vec") = MetrObj.MetropLambda1Vec,
                                      Rcpp::Named("AcceptanceEtaVec") = MetrObj.AcceptanceEtaVec,
                                      Rcpp::Named("MetropEtaVec") = MetrObj.MetropEtaVec,
                                      Rcpp::Named("AcceptanceAlpha") = MetrObj.AcceptanceAlpha,
                                      Rcpp::Named("MetropAlpha") = MetrObj.MetropAlpha);
  return Out;

}



//Output Metropolis object for summary-------------------------------------------------------------------
Rcpp::List OutputMetrObj_lmc(metrobj_lmc MetrObj) {

  Rcpp::List Out = Rcpp::List::create(Rcpp::Named("AcceptanceBeta0Vec") = MetrObj.AcceptanceBeta0Vec,
                                      Rcpp::Named("MetropBeta0Vec") = MetrObj.MetropBeta0Vec,
                                      Rcpp::Named("AcceptanceBeta1Vec") = MetrObj.AcceptanceBeta1Vec,
                                      Rcpp::Named("MetropBeta1Vec") = MetrObj.MetropBeta1Vec,
                                      Rcpp::Named("AcceptanceLambda0Vec") = MetrObj.AcceptanceLambda0Vec,
                                      Rcpp::Named("MetropLambda0Vec") = MetrObj.MetropLambda0Vec,
                                      Rcpp::Named("AcceptanceLambda1Vec") = MetrObj.AcceptanceLambda1Vec,
                                      Rcpp::Named("MetropLambda1Vec") = MetrObj.MetropLambda1Vec,
                                      Rcpp::Named("AcceptanceEtaVec") = MetrObj.AcceptanceEtaVec,
                                      Rcpp::Named("MetropEtaVec") = MetrObj.MetropEtaVec,
                                      Rcpp::Named("AcceptanceSigma") = MetrObj.AcceptanceSigma,
                                      Rcpp::Named("MetropSigma") = MetrObj.MetropSigma,
                                      Rcpp::Named("AcceptanceAlpha") = MetrObj.AcceptanceAlpha,
                                      Rcpp::Named("MetropAlpha") = MetrObj.MetropAlpha);
  return Out;

}



//Output Metropolis object for summary-------------------------------------------------------------------
Rcpp::List OutputMetrObj_novar(metrobj_novar MetrObj) {

  Rcpp::List Out = Rcpp::List::create(Rcpp::Named("AcceptanceLambdaVec") = MetrObj.AcceptanceLambdaVec,
                                      Rcpp::Named("MetropLambdaVec") = MetrObj.MetropLambdaVec,
                                      Rcpp::Named("AcceptanceEtaVec") = MetrObj.AcceptanceEtaVec,
                                      Rcpp::Named("MetropEtaVec") = MetrObj.MetropEtaVec,
                                      Rcpp::Named("AcceptanceAlpha") = MetrObj.AcceptanceAlpha,
                                      Rcpp::Named("MetropAlpha") = MetrObj.MetropAlpha);
  return Out;

}



//Initiate burn-in progress bar-------------------------------------------------------------------------------------
void SamplerProgress(int s, mcmcobj McmcObj) {

  //Set MCMC object
  int NSims = McmcObj.NSims;
  int NBurn = McmcObj.NBurn;

  //Add a new percentage
  Rcpp::Rcout.precision(0);
  Rcpp::Rcout << std::fixed << 100 * (s - NBurn) / NSims << "%..  ";

}



//Function for storing raw MCMC samples to to an object in memory
arma::colvec StoreSamples(datobj DatObj, para Para) {

  //Set data object
  int M = DatObj.M;

  //Set parameter objects
  double Alpha = Para.Alpha;
  arma::vec Delta = Para.Delta;
  arma::mat Sigma = Para.Sigma;
  arma::vec Phi = Para.Phi;

  //Save raw samples
  int counter = 0;
  arma::colvec col(5 + 15 + 1 + 5 * M);
  col(0) = Alpha;
  for (int i = 0; i < 5; i++) col(i + 1) = Delta(i);
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j <= i; j++) {
      col(6 + counter) = Sigma(i, j);
      counter++;
    }
  }
  for (int i = 0; i < (5 * M); i++) col(i + 21) = Phi(i);
  return col;
}



//Function for storing raw MCMC samples to to an object in memory
arma::colvec StoreSamples_lmc(datobj_lmc DatObj, para_lmc Para) {

  //Set data object
  int M = DatObj.M;

  //Set parameter objects
  arma::vec Alpha = Para.Alpha;
  arma::vec Delta = Para.Delta;
  arma::mat Sigma = Para.Sigma;
  arma::vec Phi = Para.Phi;

  //Save raw samples
  int counter = 0;
  arma::colvec col(5 + 15 + 5 + 5 * M);
  for (int i = 0; i < 5; i++) col(i) = Alpha(i);
  for (int i = 0; i < 5; i++) col(i + 5) = Delta(i);
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j <= i; j++) {
      col(10 + counter) = Sigma(i, j);
      counter++;
    }
  }
  for (int i = 0; i < (5 * M); i++) col(i + 25) = Phi(i);
  return col;
}



//Function for storing raw MCMC samples to to an object in memory
arma::colvec StoreSamples_novar(datobj_novar DatObj, para_novar Para) {

  //Set data object
  int M = DatObj.M;

  //Set parameter objects
  double Alpha = Para.Alpha;
  arma::vec Delta = Para.Delta;
  arma::mat Sigma = Para.Sigma;
  arma::vec Phi = Para.Phi;

  //Save raw samples
  int counter = 0;
  arma::colvec col(4 + 10 + 1 + 4 * M);
  col(0) = Alpha;
  for (int i = 0; i < 4; i++) col(i + 1) = Delta(i);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j <= i; j++) {
      col(5 + counter) = Sigma(i, j);
      counter++;
    }
  }
  for (int i = 0; i < (4 * M); i++) col(i + 15) = Phi(i);
  return col;
}



// //Update burn-in progress bar----------------------------------------------------------------------------
// void UpdateBurnInBar(int s, mcmcobj McmcObj) {
//
//   //Set MCMC object
//   arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
//   int BarLength = McmcObj.BarLength;
//
//   //Add a new star
//   arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
//   arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
//   int NewStar = NewStarBooleanVec(0);
//   for (int i = 0; i < (BarLength + 1 - NewStar); i++) Rcpp::Rcout << std::fixed << "\b";
//   Rcpp::Rcout << std::fixed << "*";
//   for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
//   Rcpp::Rcout << std::fixed << "|";
//
// }




//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateBurnInBarInt(int s, mcmcobj McmcObj) {

  //Set MCMC object
  arma::vec WhichBurnInProgressInt = McmcObj.WhichBurnInProgressInt;
  arma::uvec NewStarBoolean = find(s == WhichBurnInProgressInt);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);

  //Add percentage to submited job mode
  Rcpp::Rcout.precision(0);
  if (NewStar == 0) Rcpp::Rcout << std::fixed << "10%..  ";
  if (NewStar == 1) Rcpp::Rcout << std::fixed << "20%..  ";
  if (NewStar == 2) Rcpp::Rcout << std::fixed << "30%..  ";
  if (NewStar == 3) Rcpp::Rcout << std::fixed << "40%..  ";
  if (NewStar == 4) Rcpp::Rcout << std::fixed << "50%..  ";
  if (NewStar == 5) Rcpp::Rcout << std::fixed << "60%..  ";
  if (NewStar == 6) Rcpp::Rcout << std::fixed << "70%..  ";
  if (NewStar == 7) Rcpp::Rcout << std::fixed << "80%..  ";
  if (NewStar == 8) Rcpp::Rcout << std::fixed << "90%..  ";
  if (NewStar == 9) Rcpp::Rcout << std::fixed << "100%!  ";

}



// //Update burn-in progress bar----------------------------------------------------------------------------
// void UpdateBurnInBar(int s, mcmcobj McmcObj) {
//
//   //Set MCMC object
//   arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
//   int BarLength = McmcObj.BarLength;
//
//   //Number of new star in interactive mode
//   arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
//   arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
//   int NewStar = NewStarBooleanVec(0);
//   for (int i = 0; i < (BarLength + 1 - NewStar); i++) Rcpp::Rcout << std::fixed << "\b";
//   Rcpp::Rcout << std::fixed << "*";
//   for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
//   Rcpp::Rcout << std::fixed << "|";
//
// }



//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateBurnInBar(int s, mcmcobj McmcObj) {

  //Set MCMC object
  arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
  int BarLength = McmcObj.BarLength;

  //Add a new star
  arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);
  Rcpp::Rcout << std::fixed << "\rBurn-in progress:  |";
  for (int i = 0; i < NewStar; i++) Rcpp::Rcout << std::fixed << "*";
  for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
  Rcpp::Rcout << std::fixed << "|";

}

