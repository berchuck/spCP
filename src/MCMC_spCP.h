//
//  functions used in spCP package
//

#ifndef __spCP__
#define __spCP__

//MCMC Sampler
Rcpp::List spCP_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                        Rcpp::List MetrObj_List, Rcpp::List Para_List,
                        Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                        arma::mat RawSamples, bool Interactive);
Rcpp::List CP_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                   Rcpp::List MetrObj_List, Rcpp::List Para_List,
                   Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                   arma::mat RawSamples, bool Interactive);
Rcpp::List spCP_lmc_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                         Rcpp::List MetrObj_List, Rcpp::List Para_List,
                         Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                         arma::mat RawSamples, bool Interactive);
Rcpp::List spCP_novar_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                           Rcpp::List MetrObj_List, Rcpp::List Para_List,
                           Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                           arma::mat RawSamples, bool Interactive);


//STRUCT DEFINITIONS
struct datobj {
  double Rho;
  double ScaleY;
  double ScaleDM;
  double tNu;
  int N;
  int M;
  int Nu;
  int FamilyInd;
  int WeightsInd;
  arma::vec YStar;
  arma::mat YStarWide;
  arma::vec DM;
  arma::Mat<int> W;
  arma::vec Time;
  arma::vec TimeVec;
  arma::vec OneM;
  arma::vec OneNu;
  arma::vec OneN;
  arma::mat EyeM;
  arma::mat EyeNu;
  arma::mat EyeN;
  arma::mat Eye5;
  arma::mat Eye5M;
  arma::mat ZDelta;
  arma::vec DMLong;
  arma::umat AdjacentEdgesBoolean;
  arma::umat PhiIndeces;
  arma::uvec XThetaInd;
};
struct hypara {
  double Kappa2;
  double Xi;
  arma::mat Psi;
  double AAlpha;
  double BAlpha;
};
struct metrobj {
  arma::vec MetropLambda0Vec;
  arma::vec AcceptanceLambda0Vec;
  arma::vec MetropLambda1Vec;
  arma::vec AcceptanceLambda1Vec;
  arma::vec MetropEtaVec;
  arma::vec AcceptanceEtaVec;
  double MetropAlpha;
  double AcceptanceAlpha;
};
struct para {
  arma::vec Beta;
  arma::vec Lambda;
  arma::vec Eta;
  arma::vec Delta;
  double Alpha;
  arma::mat Sigma;
  arma::vec Sigma2;
  arma::mat Omega;
  arma::mat OmegaInv;
  arma::mat WAlpha;
  arma::mat QInv;
  arma::mat Q;
  arma::mat SigmaInv;
  arma::vec Theta;
  arma::mat XTheta;
  arma::vec Mu;
  arma::vec Phi;
  arma::mat PhiPrec;
  arma::mat PhiCov;
  arma::vec PhiMean;
};
struct dataug {
  int NBelow;
  int NAbove;
  arma::uvec WhichAbove;
  arma::uvec WhichBelow;
};
struct mcmcobj {
  int NBurn;
  int NSims;
  int NThin;
  int NPilot;
  int NTotal;
  int NKeep;
  arma::vec WhichKeep;
  arma::vec WhichPilotAdapt;
  arma::vec WhichBurnInProgress;
  arma::vec WhichBurnInProgressInt;
  arma::vec WhichSamplerProgress;
  arma::vec BurnInProgress;
  int BarLength;
  int PilotAdaptDenominator;
};
struct datobj_lmc {
  double Rho;
  double ScaleY;
  double ScaleDM;
  double tNu;
  int N;
  int M;
  int Nu;
  int FamilyInd;
  int WeightsInd;
  arma::vec YStar;
  arma::mat YStarWide;
  arma::vec DM;
  arma::Mat<int> W;
  arma::vec Time;
  arma::vec TimeVec;
  arma::vec OneM;
  arma::vec OneNu;
  arma::vec OneN;
  arma::mat EyeM;
  arma::mat EyeNu;
  arma::mat EyeN;
  arma::mat Eye5;
  arma::mat Eye5M;
  arma::mat ZDelta;
  arma::vec DMLong;
  arma::uvec AdjacentEdgesBoolean;
  arma::umat PhiIndeces;
  arma::uvec XThetaInd;
  arma::umat GammaIndeces;
};
struct hypara_lmc {
  double Kappa2;
  double Xi;
  arma::mat Psi;
  arma::vec AAlpha;
  arma::vec BAlpha;
};
struct metrobj_lmc {
  arma::vec MetropBeta0Vec;
  arma::vec AcceptanceBeta0Vec;
  arma::vec MetropBeta1Vec;
  arma::vec AcceptanceBeta1Vec;
  arma::vec MetropLambda0Vec;
  arma::vec AcceptanceLambda0Vec;
  arma::vec MetropLambda1Vec;
  arma::vec AcceptanceLambda1Vec;
  arma::vec MetropEtaVec;
  arma::vec AcceptanceEtaVec;
  arma::vec MetropSigma;
  arma::vec AcceptanceSigma;
  arma::vec MetropAlpha;
  arma::vec AcceptanceAlpha;
};
struct para_lmc {
  arma::vec Beta0;
  arma::vec Beta1;
  arma::vec Lambda0;
  arma::vec Lambda1;
  arma::vec Eta;
  arma::vec Delta;
  arma::vec Alpha;
  arma::mat Sigma;
  arma::vec Sigma2;
  arma::mat Omega;
  arma::mat OmegaInv;
  arma::mat Gamma;
  arma::mat GammaInv;
  arma::mat SigmaInv;
  arma::mat A;
  arma::vec Theta;
  arma::mat XTheta;
  arma::vec Mu;
  arma::vec Phi;
  arma::mat PhiPrec;
  arma::mat PhiCov;
  arma::vec PhiMean;
};
// struct dataug_lmc {
//   int NBelow;
//   int NAbove;
//   arma::uvec WhichAbove;
//   arma::uvec WhichBelow;
// };
// struct mcmcobj_lmc {
//   int NBurn;
//   int NSims;
//   int NThin;
//   int NPilot;
//   int NTotal;
//   int NKeep;
//   arma::vec WhichKeep;
//   arma::vec WhichPilotAdapt;
//   arma::vec WhichBurnInProgress;
//   arma::vec WhichBurnInProgressInt;
//   arma::vec WhichSamplerProgress;
//   arma::vec BurnInProgress;
//   int BarLength;
//   int PilotAdaptDenominator;
// };
struct datobj_novar {
  double Rho;
  double ScaleY;
  double ScaleDM;
  double tNu;
  int N;
  int M;
  int Nu;
  int FamilyInd;
  int WeightsInd;
  arma::vec YStar;
  arma::mat YStarWide;
  arma::vec DM;
  arma::Mat<int> W;
  arma::vec Time;
  arma::vec TimeVec;
  arma::vec OneM;
  arma::vec OneNu;
  arma::vec OneN;
  arma::mat EyeM;
  arma::mat EyeNu;
  arma::mat EyeN;
  arma::mat Eye4;
  arma::mat Eye4M;
  arma::mat ZDelta;
  arma::vec DMLong;
  arma::uvec AdjacentEdgesBoolean;
  arma::umat PhiIndeces;
  arma::uvec XThetaInd;
};
struct hypara_novar {
  double Kappa2;
  double Xi;
  arma::mat Psi;
  double AAlpha;
  double BAlpha;
};
struct metrobj_novar {
  arma::vec MetropLambdaVec;
  arma::vec AcceptanceLambdaVec;
  arma::vec MetropEtaVec;
  arma::vec AcceptanceEtaVec;
  double MetropAlpha;
  double AcceptanceAlpha;
};
struct para_novar {
  arma::vec Beta;
  arma::vec Lambda;
  arma::vec Eta;
  arma::vec Delta;
  double Alpha;
  arma::mat Sigma;
  arma::vec Sigma2;
  arma::mat Omega;
  arma::mat OmegaInv;
  arma::mat WAlpha;
  arma::mat QInv;
  arma::mat Q;
  arma::mat SigmaInv;
  arma::vec Theta;
  arma::mat XTheta;
  arma::vec Mu;
  arma::vec Phi;
  arma::mat PhiPrec;
  arma::mat PhiCov;
  arma::vec PhiMean;
};
// struct dataug_novar {
//   int NBelow;
//   int NAbove;
//   arma::uvec WhichAbove;
//   arma::uvec WhichBelow;
// };
// struct mcmcobj_novar {
//   int NBurn;
//   int NSims;
//   int NThin;
//   int NPilot;
//   int NTotal;
//   int NKeep;
//   arma::vec WhichKeep;
//   arma::vec WhichPilotAdapt;
//   arma::vec WhichBurnInProgress;
//   arma::vec WhichBurnInProgressInt;
//   arma::vec WhichSamplerProgress;
//   arma::vec BurnInProgress;
//   int BarLength;
//   int PilotAdaptDenominator;
// };

//COVARIANCE FUNCTIONS
arma::mat GetRooti(arma::mat const& Cov, arma::mat const& Eye);
arma::mat WAlphaFnc(double Alpha, arma::colvec const& DMLong, arma::umat const& AdjacentEdgesBoolean, arma::Mat<int> const& W, int M, int WeightsInd);
arma::mat QFnc(arma::mat const& WAlpha, arma::mat const& EyeM, double Rho, int M);
arma::mat QInvFnc(arma::mat const& WAlpha, arma::mat const& EyeM, double Rho, int M);

//DISTRIBUTION FUNCTIONS
double lndMvn(arma::vec const& Y, arma::vec const& Mu, arma::mat const& Rooti);
double randuRcpp();
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);
arma::vec rnormRcpp(int n, double mean, double sd);
double rtnormRcpp(double mean, double sd, bool Above);
// double rtnormRcppMSM(double mean, double sd, double lower, double upper);
arma::vec rtnormRcppMSM(int N, arma::vec const& mean, arma::vec const& sd, double lower, double upper);
arma::mat rwishRcpp(double n, arma::mat const& V);

//MCMC CONVERSION FUNCTIONS
datobj ConvertDatObj(Rcpp::List DatObj_List);
hypara ConvertHyPara(Rcpp::List HyPara_List);
metrobj ConvertMetrObj(Rcpp::List MetrObj_List);
para ConvertPara(Rcpp::List Para_List);
mcmcobj ConvertMcmcObj(Rcpp::List McmcObj_List);
dataug ConvertDatAug(Rcpp::List DatAug_List);
datobj_lmc ConvertDatObj_lmc(Rcpp::List DatObj_List);
hypara_lmc ConvertHyPara_lmc(Rcpp::List HyPara_List);
metrobj_lmc ConvertMetrObj_lmc(Rcpp::List MetrObj_List);
para_lmc ConvertPara_lmc(Rcpp::List Para_List);
// mcmcobj_lmc ConvertMcmcObj_lmc(Rcpp::List McmcObj_List);
// dataug_lmc ConvertDatAug_lmc(Rcpp::List DatAug_List);
datobj_novar ConvertDatObj_novar(Rcpp::List DatObj_List);
hypara_novar ConvertHyPara_novar(Rcpp::List HyPara_List);
metrobj_novar ConvertMetrObj_novar(Rcpp::List MetrObj_List);
para_novar ConvertPara_novar(Rcpp::List Para_List);
// mcmcobj_novar ConvertMcmcObj_novar(Rcpp::List McmcObj_List);
// dataug_novar ConvertDatAug_novar(Rcpp::List DatAug_List);

//MCMC SAMPLER FUNCTIONS
para SampleDelta(datobj DatObj, para Para, hypara HyPara);
std::pair<para, metrobj> SampleAlpha(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj);
para SampleSigma(datobj DatObj, para Para, hypara HyPara);
para SampleBeta(datobj DatObj, para Para);
std::pair<para, metrobj> SampleLambda0(datobj DatObj, para Para, metrobj MetrObj);
std::pair<para, metrobj> SampleLambda1(datobj DatObj, para Para, metrobj MetrObj);
std::pair<para, metrobj> SampleEta(datobj DatObj, para Para, metrobj MetrObj);
arma::vec SampleProbit(datobj DatObj, para Para, dataug DatAug);
arma::vec SampleProbit_lmc(datobj DatObj, para_lmc Para, dataug DatAug);
arma::vec SampleTobit(datobj DatObj, para Para, dataug DatAug);
arma::vec SampleTobit_lmc(datobj DatObj, para_lmc Para, dataug DatAug);
datobj SampleY(datobj DatObj, para Para, dataug DatAug);
datobj SampleY_lmc(datobj DatObj, para_lmc Para, dataug DatAug);

//MCMC SAMPLER FUNCTIONS FOR LMC MODEL
para_lmc SampleDelta_lmc(datobj_lmc DatObj, para_lmc Para, hypara_lmc HyPara);
std::pair<para_lmc, metrobj_lmc> SampleBeta0_lmc(datobj_lmc DatObj, para_lmc Para, metrobj_lmc MetrObj);
std::pair<para_lmc, metrobj_lmc> SampleBeta1_lmc(datobj_lmc DatObj, para_lmc Para, metrobj_lmc MetrObj);
std::pair<para_lmc, metrobj_lmc> SampleLambda0_lmc(datobj_lmc DatObj, para_lmc Para, metrobj_lmc MetrObj);
std::pair<para_lmc, metrobj_lmc> SampleLambda1_lmc(datobj_lmc DatObj, para_lmc Para, metrobj_lmc MetrObj);
std::pair<para_lmc, metrobj_lmc> SampleEta_lmc(datobj_lmc DatObj, para_lmc Para, metrobj_lmc MetrObj);
std::pair<para_lmc, metrobj_lmc> SampleSigma_lmc(datobj_lmc DatObj, para_lmc Para, hypara_lmc HyPara, metrobj_lmc MetrObj);
std::pair<para_lmc, metrobj_lmc> SampleAlpha_lmc(datobj_lmc DatObj, para_lmc Para, hypara_lmc HyPara, metrobj_lmc MetrObj);
arma::vec SampleProbit_lmc(datobj_lmc DatObj, para_lmc Para, dataug DatAug);
arma::vec SampleTobit_lmc(datobj_lmc DatObj, para_lmc Para, dataug DatAug);
datobj_lmc SampleY_lmc(datobj_lmc DatObj, para_lmc Para, dataug DatAug);

//MCMC SAMPLER FUNCTIONS FOR NO VARIANCE MODEL
para_novar SampleDelta_novar(datobj_novar DatObj, para_novar Para, hypara_novar HyPara);
std::pair<para_novar, metrobj_novar> SampleAlpha_novar(datobj_novar DatObj, para_novar Para, hypara_novar HyPara, metrobj_novar MetrObj);
para_novar SampleSigma_novar(datobj_novar DatObj, para_novar Para, hypara_novar HyPara);
para_novar SampleBeta_novar(datobj_novar DatObj, para_novar Para);
std::pair<para_novar, metrobj_novar> SampleLambda_novar(datobj_novar DatObj, para_novar Para, metrobj_novar MetrObj);
std::pair<para_novar, metrobj_novar> SampleEta_novar(datobj_novar DatObj, para_novar Para, metrobj_novar MetrObj);
arma::vec SampleProbit_novar(datobj_novar DatObj, para_novar Para, dataug DatAug);
arma::vec SampleTobit_novar(datobj_novar DatObj, para_novar Para, dataug DatAug);
datobj_novar SampleY_novar(datobj_novar DatObj, para_novar Para, dataug DatAug);

//MCMC UTILITY FUNCTIONS
void BeginBurnInProgress(mcmcobj McmcObj, bool Interactive);
Rcpp::List OutputMetrObj(metrobj MetrObj);
Rcpp::List OutputMetrObj_lmc(metrobj_lmc MetrObj);
Rcpp::List OutputMetrObj_novar(metrobj_novar MetrObj);
metrobj PilotAdaptation(datobj DatObj, metrobj MetrObj, mcmcobj McmcObj);
metrobj_lmc PilotAdaptation_lmc(datobj_lmc DatObj, metrobj_lmc MetrObj, mcmcobj McmcObj);
metrobj_novar PilotAdaptation_novar(datobj_novar DatObj, metrobj_novar MetrObj, mcmcobj McmcObj);
void SamplerProgress(int s, mcmcobj McmcObj);
arma::colvec StoreSamples(datobj DatObj, para Para);
arma::colvec StoreSamples_lmc(datobj_lmc DatObj, para_lmc Para);
arma::colvec StoreSamples_novar(datobj_novar DatObj, para_novar Para);
void UpdateBurnInBar(int s, mcmcobj McmcObj);
void UpdateBurnInBarInt(int s, mcmcobj McmcObj);

//UTILITY FUNCTIONS
arma::mat GetXThetaLoc(double ThetaLoc, arma::vec const& Time, arma::vec const& OneNu, int Nu);
arma::colvec GetXThetaLoc_lmc(double ThetaLoc, arma::vec const& Time, arma::vec const& OneNu, int Nu);
arma::vec CreatePhi(arma::vec const& Beta, arma::vec const& Lambda, arma::vec const& Eta, int M);
arma::vec CreatePhi_lmc(arma::vec const& Beta0, arma::vec const& Beta1, arma::vec const& Lambda0, arma::vec const& Lambda1, arma::vec const& Eta, int M);
arma::vec CreatePhi_novar(arma::vec const& Beta, arma::vec const& Lambda, arma::vec const& Eta, int M);
arma::mat CholInv(arma::mat const& Cov);
arma::mat Inv2(arma::mat const& A);
arma::mat Inv3(arma::mat const& A);
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol);

#endif // __spCP__
