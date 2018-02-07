#include <RcppArmadillo.h>
#include "MCMC_spCP.h"

//Function to sample delta using a Gibbs sampler step---------------------------------------------------------------
para_lmc SampleDelta_lmc(datobj_lmc DatObj, para_lmc Para, hypara_lmc HyPara) {

  //Set data objects
  arma::mat Eye5 = DatObj.Eye5;
  arma::mat EyeM = DatObj.EyeM;
  arma::mat ZDelta = DatObj.ZDelta;
  arma::mat OneM = DatObj.OneM;

  //Set parameters
  arma::vec Phi = Para.Phi;
  arma::mat PhiPrec = Para.PhiPrec;
  arma::mat A = Para.A;
  arma::mat Gamma = Para.Gamma;

  //Set hyperparameter objects
  double Kappa2 = HyPara.Kappa2;

  //Sample delta
  arma::mat tZDelta = arma::trans(ZDelta);
  arma::mat CovDelta = CholInv(tZDelta * PhiPrec * ZDelta + Eye5 / Kappa2);
  arma::vec MeanDelta = CovDelta * (tZDelta * PhiPrec * Phi);
  arma::vec Delta = rmvnormRcpp(1, MeanDelta, CovDelta);

  //Update parameters dependent on delta
  arma::vec PhiMean = ZDelta * Delta;

  //Update parameters object
  Para.Delta = Delta;
  Para.PhiMean = PhiMean;
  return Para;
}



//Function to sample new value of beta0 using a Metropolis sampler step--------------------------------
std::pair<para_lmc, metrobj_lmc> SampleBeta0_lmc(datobj_lmc DatObj, para_lmc Para, metrobj_lmc MetrObj) {

  //Set data objects
  int M = DatObj.M;
  int Nu = DatObj.Nu;
  arma::umat PhiIndeces = DatObj.PhiIndeces;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  arma::vec OneNu = DatObj.OneNu;
  arma::mat YStarWide = arma::trans(DatObj.YStarWide);
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat Eye5M = DatObj.Eye5M;

  //Set Metropolis Tuning Objects
  arma::vec MetropBeta0Vec = MetrObj.MetropBeta0Vec;
  arma::vec AcceptanceBeta0Vec = MetrObj.AcceptanceBeta0Vec;

  //Set parameter objects
  arma::vec Phi = Para.Phi;
  arma::vec Beta0 = Para.Beta0;
  arma::vec Beta1 = Para.Beta1;
  arma::vec Eta = Para.Eta;
  arma::mat Sigma = Para.Sigma;
  arma::vec Delta = Para.Delta;
  arma::vec Sigma2 = Para.Sigma2;
  arma::mat MuMatrix = arma::trans(arma::reshape(Para.Mu, M, Nu));
  arma::vec Lambda0 = Para.Lambda0;
  arma::vec Lambda1 = Para.Lambda1;
  arma::mat XTheta = Para.XTheta;
  arma::mat PhiMean = Para.PhiMean;
  arma::mat PhiCov = Para.PhiCov;

  //Initialize objects
  arma::uvec LocInd(1), RowIndeces(Nu);
  arma::colvec YStarWideLoc(Nu), MeanLikelihoodLoc(Nu), MeanLikelihoodLocProposal(Nu);
  arma::colvec Beta0Proposal(M), PhiProposal(5 * M), XThetaLoc(Nu);;
  arma::mat OmegaLocChol(Nu, Nu), OmegaLocProposalChol(Nu, Nu);
  arma::vec RandU(1);
  arma::mat RootiLikelihood(Nu, Nu), RootiLikelihoodProposal(Nu, Nu), RootiPhi(5 * M, 5 * M);
  double Component1A, Component1B, Component1;
  double Component2A, Component2B, Component2;
  double TuningSD, LogR;
  double Beta0Loc, Beta0LocProposal, Beta1Loc;

  //Indeces
  arma::vec SeqNu(Nu);
  for (int i = 0; i < Nu; i++) SeqNu(i) = i;

  //Loop over visits
  for (int Loc = 0; Loc < M; Loc++) {

    //Indeces
    RowIndeces = arma::conv_to<arma::uvec>::from(SeqNu * M + Loc);

    //Visit specific objects
    LocInd(0) = Loc;
    YStarWideLoc = YStarWide.col(Loc);
    Beta0Loc = arma::as_scalar(Beta0.row(Loc));
    Beta1Loc = arma::as_scalar(Beta1.row(Loc));
    MeanLikelihoodLoc = MuMatrix.col(Loc);
    OmegaLocChol = arma::diagmat(sqrt(Sigma2(RowIndeces)));
    TuningSD = sqrt(MetropBeta0Vec(Loc));
    XThetaLoc = XTheta.submat(RowIndeces, LocInd);

    //Numerical fix for when a propopsal is rounded to -Inf or Inf
    Beta0LocProposal = arma::datum::inf;
    Beta0Proposal = Beta0;
    while ((!arma::is_finite(Beta0LocProposal))) {

      //Sample proposal
      Beta0LocProposal = arma::as_scalar(rnormRcpp(1, Beta0Loc, TuningSD));
      Beta0Proposal.row(Loc) = Beta0LocProposal;
      MeanLikelihoodLocProposal = Beta0LocProposal * OneNu + Beta1Loc * XThetaLoc;
      PhiProposal = CreatePhi_lmc(Beta0Proposal, Beta1, Lambda0, Lambda1, Eta, M);

    }

    //Likelihood Component
    RootiLikelihood = arma::solve(arma::trimatu(OmegaLocChol), EyeNu);
    Component1A = lndMvn(YStarWideLoc, MeanLikelihoodLocProposal, RootiLikelihood);
    Component1B = lndMvn(YStarWideLoc, MeanLikelihoodLoc, RootiLikelihood);
    Component1 = Component1A - Component1B;

    //Prior components
    RootiPhi = arma::solve(arma::trimatu(arma::chol(PhiCov)), Eye5M);
    Component2A = lndMvn(PhiProposal, PhiMean, RootiPhi);
    Component2B = lndMvn(Phi, PhiMean, RootiPhi);
    Component2 = Component2A - Component2B;

    //Log acceptance ratio
    LogR = Component1 + Component2;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep count of acceptances
      AcceptanceBeta0Vec(Loc)++;

      //Update parameters output
      Beta0 = Beta0Proposal;
      Phi = PhiProposal;

    }

    //End loop over locations
  }

  //Update Metropolis object
  MetrObj.AcceptanceBeta0Vec = AcceptanceBeta0Vec;

  //Update Para objects
  Para.Beta0 = Beta0;
  Para.Phi = Phi;
  Para.Mu = arma::kron(OneNu, Beta0) + XTheta * Beta1;

  //Return final object
  return std::pair<para_lmc, metrobj_lmc>(Para, MetrObj);

}



//Function to sample new value of beta1 using a Metropolis sampler step--------------------------------
std::pair<para_lmc, metrobj_lmc> SampleBeta1_lmc(datobj_lmc DatObj, para_lmc Para, metrobj_lmc MetrObj) {

  //Set data objects
  int M = DatObj.M;
  int Nu = DatObj.Nu;
  arma::umat PhiIndeces = DatObj.PhiIndeces;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  arma::vec OneNu = DatObj.OneNu;
  arma::mat YStarWide = arma::trans(DatObj.YStarWide);
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat Eye5M = DatObj.Eye5M;

  //Set Metropolis Tuning Objects
  arma::vec MetropBeta1Vec = MetrObj.MetropBeta1Vec;
  arma::vec AcceptanceBeta1Vec = MetrObj.AcceptanceBeta1Vec;

  //Set parameter objects
  arma::vec Phi = Para.Phi;
  arma::vec Beta0 = Para.Beta0;
  arma::vec Beta1 = Para.Beta1;
  arma::vec Eta = Para.Eta;
  arma::mat Sigma = Para.Sigma;
  arma::vec Delta = Para.Delta;
  arma::vec Sigma2 = Para.Sigma2;
  arma::mat MuMatrix = arma::trans(arma::reshape(Para.Mu, M, Nu));
  arma::vec Lambda0 = Para.Lambda0;
  arma::vec Lambda1 = Para.Lambda1;
  arma::mat XTheta = Para.XTheta;
  arma::mat PhiMean = Para.PhiMean;
  arma::mat PhiCov = Para.PhiCov;

  //Initialize objects
  arma::uvec LocInd(1), RowIndeces(Nu);
  arma::colvec YStarWideLoc(Nu), MeanLikelihoodLoc(Nu), MeanLikelihoodLocProposal(Nu);
  arma::colvec Beta1Proposal(M), PhiProposal(5 * M), XThetaLoc(Nu);;
  arma::mat OmegaLocChol(Nu, Nu), OmegaLocProposalChol(Nu, Nu);
  arma::vec RandU(1);
  arma::mat RootiLikelihood(Nu, Nu), RootiLikelihoodProposal(Nu, Nu), RootiPhi(5 * M, 5 * M);
  double Component1A, Component1B, Component1;
  double Component2A, Component2B, Component2;
  double TuningSD, LogR;
  double Beta0Loc, Beta1LocProposal, Beta1Loc;

  //Indeces
  arma::vec SeqNu(Nu);
  for (int i = 0; i < Nu; i++) SeqNu(i) = i;

  //Loop over visits
  for (int Loc = 0; Loc < M; Loc++) {

    //Indeces
    RowIndeces = arma::conv_to<arma::uvec>::from(SeqNu * M + Loc);

    //Visit specific objects
    LocInd(0) = Loc;
    YStarWideLoc = YStarWide.col(Loc);
    Beta0Loc = arma::as_scalar(Beta0.row(Loc));
    Beta1Loc = arma::as_scalar(Beta1.row(Loc));
    MeanLikelihoodLoc = MuMatrix.col(Loc);
    OmegaLocChol = arma::diagmat(sqrt(Sigma2(RowIndeces)));
    TuningSD = sqrt(MetropBeta1Vec(Loc));
    XThetaLoc = XTheta.submat(RowIndeces, LocInd);

    //Numerical fix for when a propopsal is rounded to -Inf or Inf
    Beta1LocProposal = arma::datum::inf;
    Beta1Proposal = Beta1;
    while ((!arma::is_finite(Beta1LocProposal))) {

      //Sample proposal
      Beta1LocProposal = arma::as_scalar(rnormRcpp(1, Beta1Loc, TuningSD));
      Beta1Proposal.row(Loc) = Beta1LocProposal;
      MeanLikelihoodLocProposal = Beta0Loc * OneNu + Beta1LocProposal * XThetaLoc;
      PhiProposal = CreatePhi_lmc(Beta0, Beta1Proposal, Lambda0, Lambda1, Eta, M);

    }

    //Likelihood Component
    RootiLikelihood = arma::solve(arma::trimatu(OmegaLocChol), EyeNu);
    Component1A = lndMvn(YStarWideLoc, MeanLikelihoodLocProposal, RootiLikelihood);
    Component1B = lndMvn(YStarWideLoc, MeanLikelihoodLoc, RootiLikelihood);
    Component1 = Component1A - Component1B;

    //Prior components
    RootiPhi = arma::solve(arma::trimatu(arma::chol(PhiCov)), Eye5M);
    Component2A = lndMvn(PhiProposal, PhiMean, RootiPhi);
    Component2B = lndMvn(Phi, PhiMean, RootiPhi);
    Component2 = Component2A - Component2B;

    //Log acceptance ratio
    LogR = Component1 + Component2;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep count of acceptances
      AcceptanceBeta1Vec(Loc)++;

      //Update parameters output
      Beta1 = Beta1Proposal;
      Phi = PhiProposal;

    }

    //End loop over locations
  }

  //Update Metropolis object
  MetrObj.AcceptanceBeta1Vec = AcceptanceBeta1Vec;

  //Update Para objects
  Para.Beta1 = Beta1;
  Para.Phi = Phi;
  Para.Mu = arma::kron(OneNu, Beta0) + XTheta * Beta1;

  //Return final object
  return std::pair<para_lmc, metrobj_lmc>(Para, MetrObj);

}



//Function to sample new value of lambda0 using a Metropolis sampler step--------------------------------
std::pair<para_lmc, metrobj_lmc> SampleLambda0_lmc(datobj_lmc DatObj, para_lmc Para, metrobj_lmc MetrObj) {

  //Set data objects
  int M = DatObj.M;
  int Nu = DatObj.Nu;
  arma::umat PhiIndeces = DatObj.PhiIndeces;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  arma::vec OneNu = DatObj.OneNu;
  arma::mat YStarWide = arma::trans(DatObj.YStarWide);
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat Eye5M = DatObj.Eye5M;

  //Set Metropolis Tuning Objects
  arma::vec MetropLambda0Vec = MetrObj.MetropLambda0Vec;
  arma::vec AcceptanceLambda0Vec = MetrObj.AcceptanceLambda0Vec;

  //Set parameter objects
  arma::vec Phi = Para.Phi;
  arma::vec Beta0 = Para.Beta0;
  arma::vec Beta1 = Para.Beta1;
  arma::vec Eta = Para.Eta;
  arma::mat Sigma = Para.Sigma;
  arma::vec Delta = Para.Delta;
  arma::vec Sigma2 = Para.Sigma2;
  arma::mat MuMatrix = arma::trans(arma::reshape(Para.Mu, M, Nu));
  arma::vec Lambda0 = Para.Lambda0;
  arma::vec Lambda1 = Para.Lambda1;
  arma::mat XTheta = Para.XTheta;
  arma::mat PhiMean = Para.PhiMean;
  arma::mat PhiCov = Para.PhiCov;

  //Initialize objects
  arma::uvec LocInd(1), RowIndeces(Nu);
  arma::colvec YStarWideLoc(Nu), MeanLikelihoodLoc(Nu), MeanLikelihoodLocProposal(Nu);
  arma::colvec Lambda0Proposal(M), PhiProposal(5 * M), XThetaLoc(Nu);;
  arma::mat OmegaLocChol(Nu, Nu), OmegaLocProposalChol(Nu, Nu);
  arma::vec RandU(1), Sigma2LocProposal(Nu);
  arma::mat RootiLikelihood(Nu, Nu), RootiLikelihoodProposal(Nu, Nu), RootiPhi(5 * M, 5 * M);
  double Component1A, Component1B, Component1;
  double Component2A, Component2B, Component2;
  double TuningSD, LogR;
  double Lambda0Loc, Lambda0LocProposal, Lambda1Loc;

  //Indeces
  arma::vec SeqNu(Nu);
  for (int i = 0; i < Nu; i++) SeqNu(i) = i;

  //Loop over visits
  for (int Loc = 0; Loc < M; Loc++) {

    //Indeces
    RowIndeces = arma::conv_to<arma::uvec>::from(SeqNu * M + Loc);

    //Visit specific objects
    LocInd(0) = Loc;
    YStarWideLoc = YStarWide.col(Loc);
    Lambda0Loc = arma::as_scalar(Lambda0.row(Loc));
    Lambda1Loc = arma::as_scalar(Lambda1.row(Loc));
    MeanLikelihoodLoc = MuMatrix.col(Loc);
    OmegaLocChol = arma::diagmat(sqrt(Sigma2(RowIndeces)));
    TuningSD = sqrt(MetropLambda0Vec(Loc));
    XThetaLoc = XTheta.submat(RowIndeces, LocInd);

    //Numerical fix for when a propopsal is rounded to -Inf or Inf
    Lambda0LocProposal = arma::datum::inf;
    Lambda0Proposal = Lambda0;
    while ((!arma::is_finite(Lambda0LocProposal))) {

      //Sample proposal
      Lambda0LocProposal = arma::as_scalar(rnormRcpp(1, Lambda0Loc, TuningSD));
      Lambda0Proposal.row(Loc) = Lambda0LocProposal;
      Sigma2LocProposal = arma::exp(2 * (Lambda0LocProposal * OneNu + Lambda1Loc * XThetaLoc));
      OmegaLocProposalChol = arma::diagmat(sqrt(Sigma2LocProposal));
      PhiProposal = CreatePhi_lmc(Beta0, Beta1, Lambda0Proposal, Lambda1, Eta, M);

    }

    //Likelihood Component
    RootiLikelihoodProposal = arma::solve(arma::trimatu(OmegaLocProposalChol), EyeNu);
    RootiLikelihood = arma::solve(arma::trimatu(OmegaLocChol), EyeNu);
    Component1A = lndMvn(YStarWideLoc, MeanLikelihoodLoc, RootiLikelihoodProposal);
    Component1B = lndMvn(YStarWideLoc, MeanLikelihoodLoc, RootiLikelihood);
    Component1 = Component1A - Component1B;

    //Prior components
    RootiPhi = arma::solve(arma::trimatu(arma::chol(PhiCov)), Eye5M);
    Component2A = lndMvn(PhiProposal, PhiMean, RootiPhi);
    Component2B = lndMvn(Phi, PhiMean, RootiPhi);
    Component2 = Component2A - Component2B;

    //Log acceptance ratio
    LogR = Component1 + Component2;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep count of acceptances
      AcceptanceLambda0Vec(Loc)++;

      //Update parameters output
      Lambda0 = Lambda0Proposal;
      Phi = PhiProposal;
      Sigma2.rows(RowIndeces) = Sigma2LocProposal;

    }

    //End loop over locations
  }

  //Update Metropolis object
  MetrObj.AcceptanceLambda0Vec = AcceptanceLambda0Vec;

  //Update Para objects
  Para.Lambda0 = Lambda0;
  Para.Sigma2 = Sigma2;
  Para.Phi = Phi;
  Para.Omega = arma::diagmat(Sigma2);
  Para.OmegaInv = arma::diagmat(1 / Sigma2);

  //Return final object
  return std::pair<para_lmc, metrobj_lmc>(Para, MetrObj);

}



//Function to sample new value of lambda1 using a Metropolis sampler step--------------------------------
std::pair<para_lmc, metrobj_lmc> SampleLambda1_lmc(datobj_lmc DatObj, para_lmc Para, metrobj_lmc MetrObj) {

  //Set data objects
  int M = DatObj.M;
  int Nu = DatObj.Nu;
  arma::umat PhiIndeces = DatObj.PhiIndeces;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  arma::vec OneNu = DatObj.OneNu;
  arma::mat YStarWide = arma::trans(DatObj.YStarWide);
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat Eye5M = DatObj.Eye5M;

  //Set Metropolis Tuning Objects
  arma::vec MetropLambda1Vec = MetrObj.MetropLambda1Vec;
  arma::vec AcceptanceLambda1Vec = MetrObj.AcceptanceLambda1Vec;

  //Set parameter objects
  arma::vec Phi = Para.Phi;
  arma::vec Beta0 = Para.Beta0;
  arma::vec Beta1 = Para.Beta1;
  arma::vec Eta = Para.Eta;
  arma::mat Sigma = Para.Sigma;
  arma::vec Delta = Para.Delta;
  arma::vec Sigma2 = Para.Sigma2;
  arma::mat MuMatrix = arma::trans(arma::reshape(Para.Mu, M, Nu));
  arma::vec Lambda0 = Para.Lambda0;
  arma::vec Lambda1 = Para.Lambda1;
  arma::mat XTheta = Para.XTheta;
  arma::mat PhiMean = Para.PhiMean;
  arma::mat PhiCov = Para.PhiCov;

  //Initialize objects
  arma::uvec LocInd(1), RowIndeces(Nu);
  arma::colvec YStarWideLoc(Nu), MeanLikelihoodLoc(Nu), MeanLikelihoodLocProposal(Nu);
  arma::colvec Lambda1Proposal(M), PhiProposal(5 * M), XThetaLoc(Nu);;
  arma::mat OmegaLocChol(Nu, Nu), OmegaLocProposalChol(Nu, Nu);
  arma::vec RandU(1), Sigma2LocProposal(Nu);
  arma::mat RootiLikelihood(Nu, Nu), RootiLikelihoodProposal(Nu, Nu), RootiPhi(5 * M, 5 * M);
  double Component1A, Component1B, Component1;
  double Component2A, Component2B, Component2;
  double TuningSD, LogR;
  double Lambda0Loc, Lambda1LocProposal, Lambda1Loc;

  //Indeces
  arma::vec SeqNu(Nu);
  for (int i = 0; i < Nu; i++) SeqNu(i) = i;

  //Loop over visits
  for (int Loc = 0; Loc < M; Loc++) {

    //Indeces
    RowIndeces = arma::conv_to<arma::uvec>::from(SeqNu * M + Loc);

    //Visit specific objects
    LocInd(0) = Loc;
    YStarWideLoc = YStarWide.col(Loc);
    Lambda0Loc = arma::as_scalar(Lambda0.row(Loc));
    Lambda1Loc = arma::as_scalar(Lambda1.row(Loc));
    MeanLikelihoodLoc = MuMatrix.col(Loc);
    OmegaLocChol = arma::diagmat(sqrt(Sigma2(RowIndeces)));
    TuningSD = sqrt(MetropLambda1Vec(Loc));
    XThetaLoc = XTheta.submat(RowIndeces, LocInd);

    //Numerical fix for when a propopsal is rounded to -Inf or Inf
    Lambda1LocProposal = arma::datum::inf;
    Lambda1Proposal = Lambda1;
    while ((!arma::is_finite(Lambda1LocProposal))) {

      //Sample proposal
      Lambda1LocProposal = arma::as_scalar(rnormRcpp(1, Lambda1Loc, TuningSD));
      Lambda1Proposal.row(Loc) = Lambda1LocProposal;
      Sigma2LocProposal = arma::exp(2 * (Lambda0Loc * OneNu + Lambda1LocProposal * XThetaLoc));
      OmegaLocProposalChol = arma::diagmat(sqrt(Sigma2LocProposal));
      PhiProposal = CreatePhi_lmc(Beta0, Beta1, Lambda0, Lambda1Proposal, Eta, M);

    }

    //Likelihood Component
    RootiLikelihoodProposal = arma::solve(arma::trimatu(OmegaLocProposalChol), EyeNu);
    RootiLikelihood = arma::solve(arma::trimatu(OmegaLocChol), EyeNu);
    Component1A = lndMvn(YStarWideLoc, MeanLikelihoodLoc, RootiLikelihoodProposal);
    Component1B = lndMvn(YStarWideLoc, MeanLikelihoodLoc, RootiLikelihood);
    Component1 = Component1A - Component1B;

    //Prior components
    RootiPhi = arma::solve(arma::trimatu(arma::chol(PhiCov)), Eye5M);
    Component2A = lndMvn(PhiProposal, PhiMean, RootiPhi);
    Component2B = lndMvn(Phi, PhiMean, RootiPhi);
    Component2 = Component2A - Component2B;

    //Log acceptance ratio
    LogR = Component1 + Component2;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep count of acceptances
      AcceptanceLambda1Vec(Loc)++;

      //Update parameters output
      Lambda1 = Lambda1Proposal;
      Phi = PhiProposal;
      Sigma2.rows(RowIndeces) = Sigma2LocProposal;

    }

    //End loop over locations
  }

  //Update Metropolis object
  MetrObj.AcceptanceLambda1Vec = AcceptanceLambda1Vec;

  //Update Para objects
  Para.Lambda1 = Lambda1;
  Para.Sigma2 = Sigma2;
  Para.Phi = Phi;
  Para.Omega = arma::diagmat(Sigma2);
  Para.OmegaInv = arma::diagmat(1 / Sigma2);

  //Return final object
  return std::pair<para_lmc, metrobj_lmc>(Para, MetrObj);

}



//Function to sample new value of eta using a Metropolis sampler step--------------------------------
std::pair<para_lmc, metrobj_lmc> SampleEta_lmc(datobj_lmc DatObj, para_lmc Para, metrobj_lmc MetrObj) {

  //Set data objects
  int M = DatObj.M;
  int Nu = DatObj.Nu;
  double tNu = DatObj.tNu;
  arma::umat PhiIndeces = DatObj.PhiIndeces;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  arma::vec OneNu = DatObj.OneNu;
  arma::mat YStarWide = arma::trans(DatObj.YStarWide);
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat Eye5M = DatObj.Eye5M;
  arma::vec Time = DatObj.Time;

  //Set Metropolis Tuning Objects
  arma::vec MetropEtaVec = MetrObj.MetropEtaVec;
  arma::vec AcceptanceEtaVec = MetrObj.AcceptanceEtaVec;

  //Set parameter objects
  arma::vec Phi = Para.Phi;
  arma::vec Beta0 = Para.Beta0;
  arma::vec Beta1 = Para.Beta1;
  arma::vec Eta = Para.Eta;
  arma::mat Sigma = Para.Sigma;
  arma::vec Delta = Para.Delta;
  arma::vec Sigma2 = Para.Sigma2;
  arma::mat MuMatrix = arma::trans(arma::reshape(Para.Mu, M, Nu));
  arma::vec Lambda0 = Para.Lambda0;
  arma::vec Lambda1 = Para.Lambda1;
  arma::mat XTheta = Para.XTheta;
  arma::mat PhiMean = Para.PhiMean;
  arma::mat PhiCov = Para.PhiCov;
  arma::vec Theta = Para.Theta;

  //Initialize objects
  arma::uvec LocInd(1), RowIndeces(Nu);
  arma::colvec YStarWideLoc(Nu), MeanLikelihoodLoc(Nu), MeanLikelihoodLocProposal(Nu);
  arma::colvec EtaProposal(M), PhiProposal(5 * M), XThetaLoc(Nu);;
  arma::mat OmegaLocChol(Nu, Nu), OmegaLocProposalChol(Nu, Nu);
  arma::vec RandU(1), Sigma2LocProposal(Nu);
  arma::mat RootiLikelihood(Nu, Nu), RootiLikelihoodProposal(Nu, Nu), RootiPhi(5 * M, 5 * M);
  double Component1A, Component1B, Component1;
  double Component2A, Component2B, Component2;
  double TuningSD, LogR;
  double EtaLocProposal, EtaLoc, ThetaLocProposal;
  double Beta0Loc, Beta1Loc, Lambda0Loc, Lambda1Loc;

  //Indeces
  arma::vec SeqNu(Nu);
  for (int i = 0; i < Nu; i++) SeqNu(i) = i;

  //Loop over visits
  for (int Loc = 0; Loc < M; Loc++) {

    //Indeces
    RowIndeces = arma::conv_to<arma::uvec>::from(SeqNu * M + Loc);

    //Visit specific objects
    LocInd(0) = Loc;
    YStarWideLoc = YStarWide.col(Loc);
    EtaLoc = arma::as_scalar(Eta.row(Loc));
    Beta0Loc = arma::as_scalar(Beta0.row(Loc));
    Beta1Loc = arma::as_scalar(Beta1.row(Loc));
    Lambda0Loc = arma::as_scalar(Lambda0.row(Loc));
    Lambda1Loc = arma::as_scalar(Lambda1.row(Loc));
    MeanLikelihoodLoc = MuMatrix.col(Loc);
    OmegaLocChol = arma::diagmat(sqrt(Sigma2(RowIndeces)));
    TuningSD = sqrt(MetropEtaVec(Loc));

    //Numerical fix for when a propopsal is rounded to -Inf or Inf
    EtaLocProposal = arma::datum::inf;
    EtaProposal = Eta;
    while ((!arma::is_finite(EtaLocProposal))) {

      //Sample proposal
      EtaLocProposal = arma::as_scalar(rnormRcpp(1, EtaLoc, TuningSD));
      EtaProposal.row(Loc) = EtaLocProposal;
      ThetaLocProposal = std::min(tNu, exp(EtaLocProposal));
      XThetaLoc = GetXThetaLoc_lmc(ThetaLocProposal, Time, OneNu, Nu);
      MeanLikelihoodLocProposal = OneNu * Beta0Loc + XThetaLoc * Beta1Loc;
      Sigma2LocProposal = arma::exp(2 * (Lambda0Loc * OneNu + Lambda1Loc * XThetaLoc));
      OmegaLocProposalChol = arma::diagmat(sqrt(Sigma2LocProposal));
      PhiProposal = CreatePhi_lmc(Beta0, Beta1, Lambda0, Lambda1, EtaProposal, M);

    }

    //Likelihood Component
    RootiLikelihoodProposal = arma::solve(arma::trimatu(OmegaLocProposalChol), EyeNu);
    RootiLikelihood = arma::solve(arma::trimatu(OmegaLocChol), EyeNu);
    Component1A = lndMvn(YStarWideLoc, MeanLikelihoodLocProposal, RootiLikelihoodProposal);
    Component1B = lndMvn(YStarWideLoc, MeanLikelihoodLoc, RootiLikelihood);
    Component1 = Component1A - Component1B;

    //Prior components
    RootiPhi = arma::solve(arma::trimatu(arma::chol(PhiCov)), Eye5M);
    Component2A = lndMvn(PhiProposal, PhiMean, RootiPhi);
    Component2B = lndMvn(Phi, PhiMean, RootiPhi);
    Component2 = Component2A - Component2B;

    //Log acceptance ratio
    LogR = Component1 + Component2;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep count of acceptances
      AcceptanceEtaVec(Loc)++;

      //Update parameters output
      Eta = EtaProposal;
      Phi = PhiProposal;
      Sigma2.rows(RowIndeces) = Sigma2LocProposal;
      Theta.row(Loc) = ThetaLocProposal;
      XTheta.submat(RowIndeces, LocInd) = XThetaLoc;

    }

    //End loop over locations
  }

  //Update Metropolis object
  MetrObj.AcceptanceEtaVec = AcceptanceEtaVec;

  //Update Para objects
  Para.Eta = Eta;
  Para.Theta = Theta;
  Para.XTheta = XTheta;
  Para.Sigma2 = Sigma2;
  Para.Phi = Phi;
  Para.Omega = arma::diagmat(Sigma2);
  Para.OmegaInv = arma::diagmat(1 / Sigma2);
  Para.Mu = arma::kron(OneNu, Beta0) + XTheta * Beta1;

  //Return final object
  return std::pair<para_lmc, metrobj_lmc>(Para, MetrObj);

}



//Function to sample sigma using a Metropolis sampler step-------------------------------------------------------------------
std::pair<para_lmc, metrobj_lmc> SampleSigma_lmc(datobj_lmc DatObj, para_lmc Para, hypara_lmc HyPara, metrobj_lmc MetrObj) {

  //Set data objects
  arma::mat EyeM = DatObj.EyeM;
  arma::mat Eye5M = DatObj.Eye5M;
  int M = DatObj.M;

  //Set parameters
  arma::mat A = Para.A;
  arma::vec Phi = Para.Phi;
  arma::vec PhiMean = Para.PhiMean;
  arma::mat PhiCov = Para.PhiCov;
  arma::mat GammaInv = Para.GammaInv;
  arma::mat Gamma = Para.Gamma;
  arma::mat Sigma = Para.Sigma;
  arma::mat SigmaInv = Para.SigmaInv;

  //Set hyperparameter objects
  double Xi = HyPara.Xi;
  arma::mat Psi = HyPara.Psi;

  //Set hyperparameter objects
  arma::vec MetropSigma = MetrObj.MetropSigma;
  arma::vec AcceptanceSigma = MetrObj.AcceptanceSigma;

  //Indeces
  arma::umat AIndeces(10, 2);
  arma::uvec VarIndeces(5), CovIndeces(10);
  VarIndeces << 0 << 5 << 9 << 12 << 14 << arma::endr;
  CovIndeces << 1 << 2 << 3 << 4 << 6 << 7 << 8 << 10 << 11 << 13 << arma::endr;
  AIndeces << 0 << 1 << arma::endr
           << 0 << 2 << arma::endr
           << 0 << 3 << arma::endr
           << 0 << 4 << arma::endr
           << 1 << 2 << arma::endr
           << 1 << 3 << arma::endr
           << 1 << 4 << arma::endr
           << 2 << 3 << arma::endr
           << 2 << 4 << arma::endr
           << 3 << 4 << arma::endr;

  //Declare objects
  arma::uvec CompInd(1), Col(1), Row(1);
  arma::umat Indeces(1, 2);
  double AComp, ACompProposal, Delta, DeltaProposal, TuningSD, LogR;
  double logDetSigmaProposal, signDetProposal, logDetSigma, signDet;
  double Component1, Component1A, Component1B;
  double Component2, Component2A, Component2B;
  double Component3, Component3A, Component3B;
  double Component4, Component4A, Component4B;
  arma::mat RootiPhiProposal(5 * M, 5 * M);
  arma::mat PhiCovProposal(5 * M, 5 * M), PhiCovCholProposal(5 * M, 5 * M);
  arma::mat AProposal(5, 5), SigmaProposal(5, 5), SigmaInvProposal(5, 5);
  arma::mat RootiPhi(5 * M, 5 * M);
  arma::vec AInsert(1), RandU(1), AdiagProp(5), Adiag(5);
  bool PhiBool = false, SigmaBool = false;

  //Loop over variance components
  for (int Comp = 0; Comp < 5; Comp++) {

    //Variance components
    CompInd(0) = Comp;
    AComp = arma::as_scalar(A.submat(CompInd, CompInd));
    Delta = log(AComp);
    TuningSD = sqrt(MetropSigma(VarIndeces(Comp)));

    //Sample a proposal
    while (!(PhiBool) | !(SigmaBool)) {
      DeltaProposal = arma::as_scalar(rnormRcpp(1, Delta, TuningSD));
      ACompProposal = exp(DeltaProposal);
      AInsert(0) = ACompProposal;
      AProposal = A;
      AProposal.submat(CompInd, CompInd) = AInsert;
      PhiCovProposal = arma::kron(AProposal, EyeM) * GammaInv * arma::kron(arma::trans(AProposal), EyeM);
      PhiBool = arma::chol(PhiCovCholProposal, PhiCovProposal);
      RootiPhiProposal = arma::solve(arma::trimatu(PhiCovCholProposal), Eye5M);
      SigmaProposal = arma::trans(AProposal) * AProposal;
      SigmaBool = arma::inv_sympd(SigmaInvProposal, SigmaProposal);
      // log_det(logDetSigmaProposal, signDetProposal, SigmaProposal);
    }

    //Spatial components
    RootiPhi = arma::solve(arma::trimatu(arma::chol(PhiCov)), Eye5M);
    Component1A = lndMvn(Phi, PhiMean, RootiPhiProposal);
    Component1B = lndMvn(Phi, PhiMean, RootiPhi);
    Component1 = Component1A - Component1B;

    //Prior components
    // log_det(logDetSigma, signDet, Sigma);
    Component2A = -0.5 * (Xi + 6) * log(det(SigmaProposal)) - 0.5 * arma::trace(Psi * SigmaInvProposal);
    Component2B = -0.5 * (Xi + 6) * log(det(Sigma)) - 0.5 * arma::trace(Psi * SigmaInv);
    Component2 = Component2A - Component2B;

    //Jacobian
    AdiagProp = AProposal.diag(0);
    Adiag = A.diag(0);
    Component3A = log(32) + 5 * log(arma::as_scalar(AdiagProp.at(0))) + 4 * log(arma::as_scalar(AdiagProp.at(1))) + 3 * log(arma::as_scalar(AdiagProp.at(2))) + 2 * log(arma::as_scalar(AdiagProp.at(3))) + log(arma::as_scalar(AdiagProp.at(4)));
    Component3B = log(32) + 5 * log(arma::as_scalar(Adiag.at(0))) + 4 * log(arma::as_scalar(Adiag.at(1))) + 3 * log(arma::as_scalar(Adiag.at(2))) + 2 * log(arma::as_scalar(Adiag.at(3))) + log(arma::as_scalar(Adiag.at(4)));
    Component3 = Component3A - Component3B;

    //Log normal transformation
    Component4A = ACompProposal;
    Component4B = AComp;
    Component4 = Component4A - Component4B;

    //Compute log acceptance ratio
    LogR = Component1 + Component2 + Component3 + Component4;

    //Reset booleans
    PhiBool = false;
    SigmaBool = false;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep count of acceptances
      AcceptanceSigma(VarIndeces(Comp))++;

      //Update parameters output
      Sigma = SigmaProposal;
      SigmaInv = SigmaInvProposal;
      A = AProposal;
      PhiCov = PhiCovProposal;

    }

  //End loop over variance components
  }

  //Loop over covariance components
  for (int Comp = 0; Comp < 10; Comp++) {

    //Variance components
    CompInd(0) = Comp;
    TuningSD = sqrt(MetropSigma(CovIndeces(Comp)));
    Indeces = AIndeces.rows(CompInd);
    Row = Indeces(0);
    Col = Indeces(1);
    AComp = arma::as_scalar(A.submat(Row, Col));

    //Sample a proposal
    while (!(PhiBool) | !(SigmaBool)) {
      ACompProposal = arma::as_scalar(rnormRcpp(1, AComp, TuningSD));
      AInsert(0) = ACompProposal;
      AProposal = A;
      AProposal.submat(Row, Col) = AInsert;
      PhiCovProposal = arma::kron(AProposal, EyeM) * GammaInv * arma::kron(arma::trans(AProposal), EyeM);
      PhiBool = arma::chol(PhiCovCholProposal, PhiCovProposal);
      RootiPhiProposal = arma::solve(arma::trimatu(PhiCovCholProposal), Eye5M);
      SigmaProposal = arma::trans(AProposal) * AProposal;
      SigmaBool = arma::inv_sympd(SigmaInvProposal, SigmaProposal);
      // log_det(logDetSigmaProposal, signDet, SigmaProposal);
    }

    //Spatial components
    RootiPhi = arma::solve(arma::trimatu(arma::chol(PhiCov)), Eye5M);
    Component1A = lndMvn(Phi, PhiMean, RootiPhiProposal);
    Component1B = lndMvn(Phi, PhiMean, RootiPhi);
    Component1 = Component1A - Component1B;

    //Prior components
    // log_det(logDetSigma, signDet, Sigma);
    Component2A = -0.5 * (Xi + 6) * log(det(SigmaProposal)) - 0.5 * arma::trace(Psi * SigmaInvProposal);
    Component2B = -0.5 * (Xi + 6) * log(det(Sigma)) - 0.5 * arma::trace(Psi * SigmaInv);
    Component2 = Component2A - Component2B;

    //Jacobian
    AdiagProp = AProposal.diag(0);
    Adiag = A.diag(0);
    Component3A = log(32) + 5 * log(arma::as_scalar(AdiagProp.at(0))) + 4 * log(arma::as_scalar(AdiagProp.at(1))) + 3 * log(arma::as_scalar(AdiagProp.at(2))) + 2 * log(arma::as_scalar(AdiagProp.at(3))) + log(arma::as_scalar(AdiagProp.at(4)));
    Component3B = log(32) + 5 * log(arma::as_scalar(Adiag.at(0))) + 4 * log(arma::as_scalar(Adiag.at(1))) + 3 * log(arma::as_scalar(Adiag.at(2))) + 2 * log(arma::as_scalar(Adiag.at(3))) + log(arma::as_scalar(Adiag.at(4)));
    Component3 = Component3A - Component3B;

    //Compute log acceptance ratio
    LogR = Component1 + Component2 + Component3;

    //Reset booleans
    PhiBool = false;
    SigmaBool = false;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep count of acceptances
      AcceptanceSigma(CovIndeces(Comp))++;

      //Update parameters output
      Sigma = SigmaProposal;
      SigmaInv = SigmaInvProposal;
      A = AProposal;
      PhiCov = PhiCovProposal;

    }

  //End loop over variance components
  }

  //Update metropolis tuning
  MetrObj.AcceptanceSigma = AcceptanceSigma;

  //Update parameters object
  Para.Sigma = Sigma;
  Para.SigmaInv = SigmaInv;
  Para.A = A;
  Para.PhiPrec = arma::kron(arma::inv(arma::trans(A)), EyeM) * Gamma * arma::kron(arma::inv(A), EyeM);
  Para.PhiCov = PhiCov;
  return std::pair<para_lmc, metrobj_lmc>(Para, MetrObj);
}



//Function to sample new value of alpha using a Metropolis sampler step-----------------------------------------------
std::pair<para_lmc, metrobj_lmc> SampleAlpha_lmc(datobj_lmc DatObj, para_lmc Para, hypara_lmc HyPara, metrobj_lmc MetrObj) {

  //Set data objects
  int M = DatObj.M;
  int WeightsInd = DatObj.WeightsInd;
  arma::vec DMLong = DatObj.DMLong;
  arma::umat AdjacentEdgesBoolean = DatObj.AdjacentEdgesBoolean;
  arma::Mat<int> W = DatObj.W;
  arma::mat EyeM = DatObj.EyeM;
  arma::mat Eye5M = DatObj.Eye5M;
  double Rho = DatObj.Rho;
  arma::umat GammaIndeces = DatObj.GammaIndeces;

  //Set hyperparameter objects
  arma::vec AAlpha = HyPara.AAlpha;
  arma::vec BAlpha = HyPara.BAlpha;

  //Set parameter objects
  arma::vec Alpha = Para.Alpha;
  arma::mat Phi = Para.Phi;
  arma::vec PhiMean = Para.PhiMean;
  arma::mat PhiCov = Para.PhiCov;
  arma::mat GammaInv = Para.GammaInv;
  arma::mat Gamma = Para.Gamma;
  arma::mat A = Para.A;

  //Set metropolis objects
  arma::vec MetropAlpha = sqrt(MetrObj.MetropAlpha);
  arma::vec AcceptanceAlpha = MetrObj.AcceptanceAlpha;

  //Initialize objects
  double Component1, Component1A, Component1B;
  double Component2, Component2A, Component2B, Component3;
  double BigDelta, BigDeltaProposal, AlphaProposal, LogR;
  arma::vec AlphaProposalVec(1), AAlphaVec(1), BAlphaVec(1), RandU(1);
  double TOL = 0.000001;
  arma::mat GammaInvProposal(5 * M, 5 * M), WAlphaProposal_p(M, M), QInv_Proposal(M, M);
  arma::mat PhiCovProposal(5 * M, 5 * M), RootiPhi(5 * M, 5 * M), RootiPhiProposal(5 * M, 5 * M);
  arma::uvec Indeces(M * M);

  //Loop over alphas
  for (int p = 0; p < 5; p++) {

    //Transform current state to real line
    BigDelta = log((Alpha(p) - AAlpha(p)) / (BAlpha(p) - Alpha(p)));

    //Sample a new Proposal
    BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropAlpha(p)));

    //Compute Phi Proposal
    AlphaProposal = (BAlpha(p) * exp(BigDeltaProposal) + AAlpha(p)) / (1 + exp(BigDeltaProposal));

    //Fix numerical issue where AlphaProposal can equal AAlpha or BAlpha
    AlphaProposalVec = AlphaProposal;
    AAlphaVec = AAlpha(p);
    BAlphaVec = BAlpha(p);
    if ((rows_equal(AlphaProposalVec, AAlphaVec, TOL)) || (rows_equal(AlphaProposalVec, BAlphaVec, TOL))) {
      if (rows_equal(AlphaProposalVec, AAlphaVec, TOL)) AlphaProposal *= 1.1;
      if (rows_equal(AlphaProposalVec, BAlphaVec, TOL)) AlphaProposal *= 0.99;
      BigDeltaProposal = log((AlphaProposal - AAlpha(p)) / (BAlpha(p) - AlphaProposal));
    }

    //Proposal alpha objects
    GammaInvProposal = GammaInv;
    WAlphaProposal_p = WAlphaFnc(AlphaProposal, DMLong, AdjacentEdgesBoolean, W, M, WeightsInd);
    QInv_Proposal = QInvFnc(WAlphaProposal_p, EyeM, Rho, M);
    Indeces = GammaIndeces.col(p);
    GammaInvProposal(Indeces) = arma::vectorise(QInv_Proposal);
    PhiCovProposal = arma::kron(A, EyeM) * GammaInvProposal * arma::kron(arma::trans(A), EyeM);

    //Phi components
    RootiPhi = arma::solve(arma::trimatu(arma::chol(PhiCov)), Eye5M);
    RootiPhiProposal = arma::solve(arma::trimatu(arma::chol(PhiCovProposal)), Eye5M);
    Component1A = lndMvn(Phi, PhiMean, RootiPhiProposal);
    Component1B = lndMvn(Phi, PhiMean, RootiPhi);
    Component1 = Component1A - Component1B;

    //Prior components
    Component2A = BigDeltaProposal;
    Component2B = BigDelta;
    Component2 = Component2A - Component2B;

    //Jacobian
    Component3 = 2 * log((1 + exp(BigDelta)) / (1 + exp(BigDeltaProposal)));

    //Compute log acceptance ratio
    LogR = Component1 + Component2 + Component3;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep Count of Acceptances
      AcceptanceAlpha(p)++;

      //Update parameters object
      Alpha(p) = AlphaProposal;
      GammaInv = GammaInvProposal;
      Gamma(Indeces) = arma::vectorise(CholInv(QInv_Proposal));
      PhiCov = PhiCovProposal;

    }

  //End loop over alphas
  }

  //Update metropolis object
  MetrObj.AcceptanceAlpha = AcceptanceAlpha;

  //Update parameter object
  Para.Alpha = Alpha;
  Para.GammaInv = GammaInv;
  Para.Gamma = Gamma;
  Para.PhiPrec = arma::kron(arma::inv(arma::trans(A)), EyeM) * Gamma * arma::kron(arma::inv(A), EyeM);
  Para.PhiCov = PhiCov;

  //Return output object
  return std::pair<para_lmc, metrobj_lmc>(Para, MetrObj);

}



//Function to sample latent probit process using Gibbs sampling step------------------------------------------------
arma::vec SampleProbit_lmc(datobj_lmc DatObj, para_lmc Para, dataug DatAug) {

  //Set data objects
  arma::vec YStar = DatObj.YStar;

  //Set parameters
  arma::vec Mu = Para.Mu;
  arma::vec Sigma2 = Para.Sigma2;

  //Set data augmentation objects
  arma::uvec WhichBelow = DatAug.WhichBelow;
  arma::uvec WhichAbove = DatAug.WhichAbove;
  int NBelow = DatAug.NBelow;
  int NAbove = DatAug.NAbove;

  //Sample latent Variable from full conditional
  YStar(WhichBelow) = rtnormRcppMSM(NBelow, Mu(WhichBelow), arma::sqrt(Sigma2(WhichBelow)), -arma::datum::inf, 0);
  YStar(WhichAbove) = rtnormRcppMSM(NAbove, Mu(WhichAbove), arma::sqrt(Sigma2(WhichAbove)), 0, arma::datum::inf);
  return YStar;

}



//Function to sample latent tobit process using Gibbs sampling step------------------------------------------------
arma::vec SampleTobit_lmc(datobj_lmc DatObj, para_lmc Para, dataug DatAug) {

  //Set data objects
  arma::vec YStar = DatObj.YStar;

  //Set parameters
  arma::vec Mu = Para.Mu;
  arma::vec Sigma2 = Para.Sigma2;

  //Set data augmentation objects
  int NBelow = DatAug.NBelow;
  arma::uvec WhichBelow = DatAug.WhichBelow;

  //Moments
  arma::vec Mean = Mu(WhichBelow);
  arma::vec SD = arma::sqrt(Sigma2(WhichBelow));

  //Sample latent Variable from full conditional
  for (int i = 0; i < NBelow; i++) YStar(WhichBelow(i)) = rtnormRcpp(Mean(i), SD(i), true);
  return YStar;

}



//Function to sample latent process from its full conditional------------------------------------------------------
datobj_lmc SampleY_lmc(datobj_lmc DatObj, para_lmc Para, dataug DatAug) {

  //Set data objects
  int FamilyInd = DatObj.FamilyInd;
  int N = DatObj.N;
  int M = DatObj.M;
  int Nu = DatObj.Nu;

  //Sample latent process
  arma::vec YStar(N);
  if (FamilyInd == 1) YStar = SampleProbit_lmc(DatObj, Para, DatAug);
  if (FamilyInd == 2) YStar = SampleTobit_lmc(DatObj, Para, DatAug);

  //Save output
  DatObj.YStar = YStar;
  DatObj.YStarWide = arma::reshape(YStar, M, Nu);
  return DatObj;

}
