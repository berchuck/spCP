###Function for reading in sampler inputs and creating a list object that contains all relavent data objects--------------------
CreateDatObj_lmc <- function(Y, DM, W, Time, Rho, ScaleY, ScaleDM, Family, Weights, Distance) {

  ###Data objects
  YObserved <- Y / ScaleY #scale observed data
  N <- length(YObserved)  #total observations
  M <- dim(W)[1] #number of spatial locations
  Nu <- N / M #number of visits

  ###Dynamic Objects (updated with data augmentation)
  YStar <- matrix(YObserved, ncol = 1)
  YStarWide <- matrix(YStar, nrow = M, ncol = Nu)

  ###Matrix Objects
  OneM <- matrix(rep(1, M), nrow = M , ncol = 1)
  OneNu <- matrix(rep(1, Nu), nrow = Nu, ncol = 1)
  OneN <- matrix(rep(1, N), nrow = N, ncol = 1)
  One5 <- matrix(rep(1, 5), nrow = 5, ncol = 1)
  EyeM <- diag(M)
  EyeNu <- diag(Nu)
  EyeN <- diag(N)
  Eye5 <- diag(5)
  Eye5M <- diag(5 * M)

  ###Create Z
  if (Distance == "euclidean") Dist <- function(x, y) abs(x - y)
  if (Distance == "circumference") Dist <- function(x, y) pmin(abs(x - y), (360 - pmax(x, y) + pmin(x, y))) #circumference of optic nerve
  DM_Grid <- expand.grid(DM, DM)
  DM_Vector <- Dist(DM_Grid[ , 1], DM_Grid[ , 2]) / ScaleDM
  DM_Matrix <- matrix(DM_Vector, nrow = dim(W)[1], ncol = dim(W)[1], byrow = TRUE)
  AdjacentEdgesBoolean <- (W == 1) & (!lower.tri(W))
  DMLong <- matrix(DM_Matrix[AdjacentEdgesBoolean], ncol = 1)
  AdjacentEdgesBoolean <- matrix(which(AdjacentEdgesBoolean) - 1, ncol = 1)

  ###Weights indicator
  if (Weights == "continuous") WeightsInd <- 0
  if (Weights == "binary") WeightsInd <- 1

  ###Family indicator
  if (Family == "normal") FamilyInd <- 0
  if (Family == "probit") FamilyInd <- 1
  if (Family == "tobit") FamilyInd <- 2

  ###Object for XTheta and XThetaLoc
  XThetaLogical <- do.call(rbind, replicate(Nu, EyeM, simplify = FALSE))
  XThetaInd <- cbind(which(XThetaLogical == 1, arr.ind = TRUE), which(XThetaLogical == 1))
  XThetaInd <- XThetaInd[order(XThetaInd[, 1]), ]
  XThetaInd <- XThetaInd[, 3] - 1

  ###Delta matrix
  ZDelta <- kronecker(Eye5, OneM)

  ###Phi indeces
  PhiIndeces <- matrix(1:(5 * M), ncol = 5)

  ###Gamma indeces
  GammaIndeces <- matrix(nrow = M * M, ncol = 5)
  Gamma <- Eye5M
  for (p in 1:5) {
    Gamma[as.matrix(expand.grid(1:M + (p - 1) * M, 1:M + (p - 1) * M))] <- p
  }
  for (p in 1:5) {
    GammaIndeces[, p] <- which(Gamma == p) - 1
  }

  ###Make parameters global
  DatObj <- list()
  DatObj$YObserved <- YObserved
  DatObj$ScaleY <- ScaleY
  DatObj$ScaleDM <- ScaleDM
  DatObj$YStar <- YStar
  DatObj$YStarWide <- YStarWide
  DatObj$DM <- DM
  DatObj$W <- W
  DatObj$N <- N
  DatObj$M <- M
  DatObj$Nu <- Nu
  DatObj$OneM <- OneM
  DatObj$OneNu <- OneNu
  DatObj$OneN <- OneN
  DatObj$EyeM <- EyeM
  DatObj$EyeNu <- EyeNu
  DatObj$EyeN <- EyeN
  DatObj$Eye5 <- Eye5
  DatObj$Eye5M <- Eye5M
  DatObj$DMLong <- DMLong
  DatObj$AdjacentEdgesBoolean <- AdjacentEdgesBoolean
  DatObj$FamilyInd <- FamilyInd
  DatObj$WeightsInd <- WeightsInd
  DatObj$Time <- Time
  DatObj$TimeVec <- kronecker(Time, OneM)
  DatObj$tNu <- max(Time)
  DatObj$Rho <- Rho
  DatObj$XThetaInd <- XThetaInd
  DatObj$ZDelta <- ZDelta
  DatObj$PhiIndeces <- PhiIndeces
  DatObj$GammaIndeces <- GammaIndeces
  return(DatObj)

}




###Function to create Hyperparameter Object------------------------------------------------------------------------------------
CreateHyPara_lmc <- function(Hypers, DatObj) {

  ###Set data objects
  DMLong <- DatObj$DMLong

  ###Which parameters are user defined?
  UserHypers <- names(Hypers)

  ###Set hyperparameters for Delta
  if ("Delta" %in% UserHypers) {
    Kappa2 <- Hypers$Delta$Kappa2
  }
  if (!("Delta" %in% UserHypers)) {
    Kappa2 <- 1000
  }

  ###Set hyperparameters for Sigma
  if ("Sigma" %in% UserHypers) {
    Xi <- Hypers$Sigma$Xi
    Psi <- Hypers$Sigma$Psi
  }
  if (!("Sigma" %in% UserHypers)) {
    Xi <- 5 + 1
    Psi <- diag(5)
  }

  ###Set hyperparameters for Alpha
  if ("Alpha" %in% UserHypers) {
    AAlpha <- Hypers$Alpha$AAlpha
    BAlpha <- Hypers$Alpha$BAlpha
  }
  if (!("Alpha" %in% UserHypers)) {
    BAlpha <- rep(-log(0.5) / min(DMLong[DMLong > 0]), 5)
    AAlpha <- rep(-log(0.5) / max(DMLong), 5)
  }

  ###Create object for hyperparameters
  HyPara <- list()
  HyPara$Kappa2 <- Kappa2
  HyPara$Xi <- Xi
  HyPara$Psi <- Psi
  HyPara$AAlpha <- AAlpha
  HyPara$BAlpha <- BAlpha
  return(HyPara)

}



###Function for creating an object containing relevant Metropolis information---------------------------------------------------
CreateMetrObj_lmc <- function(Tuning, DatObj) {

  ###Set Data Objects
  M <- DatObj$M

  ###Which parameters are user defined?
  UserTuners <- names(Tuning)

  ###Set tuning parameters for Beta0Vec
  if ("Beta0Vec" %in% UserTuners) MetropBeta0Vec <- Tuning$Beta0Vec
  if (!("Beta0Vec" %in% UserTuners)) MetropBeta0Vec <- rep(1, M)

  ###Set tuning parameters for Beta1Vec
  if ("Beta1Vec" %in% UserTuners) MetropBeta1Vec <- Tuning$Beta1Vec
  if (!("Beta1Vec" %in% UserTuners)) MetropBeta1Vec <- rep(1, M)

  ###Set tuning parameters for Lambda0Vec
  if ("Lambda0Vec" %in% UserTuners) MetropLambda0Vec <- Tuning$Lambda0Vec
  if (!("Lambda0Vec" %in% UserTuners)) MetropLambda0Vec <- rep(1, M)

  ###Set tuning parameters for Lambda1Vec
  if ("Lambda1Vec" %in% UserTuners) MetropLambda1Vec <- Tuning$Lambda1Vec
  if (!("Lambda1Vec" %in% UserTuners)) MetropLambda1Vec <- rep(1, M)

  ###Set tuning parameters for EtaVec
  if ("EtaVec" %in% UserTuners) MetropEtaVec <- Tuning$EtaVec
  if (!("EtaVec" %in% UserTuners)) MetropEtaVec <- rep(1, M)

  ###Set tuning parameters for Sigma
  if ("Sigma" %in% UserTuners) MetropSigma <- Tuning$Sigma
  if (!("Sigma" %in% UserTuners)) MetropSigma <- rep(1, 15)

  ###Set tuning parameters for Alpha
  if ("Alpha" %in% UserTuners) MetropAlpha <- Tuning$Alpha
  if (!("Alpha" %in% UserTuners)) MetropAlpha <- rep(1, 5)

  ###Set acceptance rate counters
  AcceptanceBeta0Vec <- AcceptanceBeta1Vec <- AcceptanceLambda0Vec <- AcceptanceLambda1Vec <- AcceptanceEtaVec <- rep(0, M)
  AcceptanceAlpha <- rep(0, 5)
  AcceptanceSigma <- rep(0, 15)

  ###Return metropolis object
  MetrObj <- list()
  MetrObj$MetropBeta0Vec <- MetropBeta0Vec
  MetrObj$AcceptanceBeta0Vec <- AcceptanceBeta0Vec
  MetrObj$MetropBeta1Vec <- MetropBeta1Vec
  MetrObj$AcceptanceBeta1Vec <- AcceptanceBeta1Vec
  MetrObj$MetropLambda0Vec <- MetropLambda0Vec
  MetrObj$AcceptanceLambda0Vec <- AcceptanceLambda0Vec
  MetrObj$MetropLambda1Vec <- MetropLambda1Vec
  MetrObj$AcceptanceLambda1Vec <- AcceptanceLambda1Vec
  MetrObj$MetropEtaVec <- MetropEtaVec
  MetrObj$AcceptanceEtaVec <- AcceptanceEtaVec
  MetrObj$MetropSigma <- MetropSigma
  MetrObj$AcceptanceSigma <- AcceptanceSigma
  MetrObj$MetropAlpha <- MetropAlpha
  MetrObj$AcceptanceAlpha <- AcceptanceAlpha
  MetrObj$OriginalTuners <- c(MetropBeta0Vec, MetropBeta1Vec, MetropLambda0Vec, MetropLambda1Vec, MetropEtaVec, MetropSigma, MetropAlpha)
  return(MetrObj)

}



###Function for creating inital parameter object-------------------------------------------------------------------------------
CreatePara_lmc <- function(Starting, DatObj, HyPara) {

  ###Set data objects
  M <- DatObj$M
  N <- DatObj$N
  tNu <- DatObj$tNu
  OneNu <- DatObj$OneNu
  OneN <- DatObj$OneN
  OneM <- DatObj$OneM
  EyeM <- DatObj$EyeM
  TimeVec <- DatObj$TimeVec
  W <- DatObj$W
  DMLong <- DatObj$DMLong
  AdjacentEdgesBoolean <- DatObj$AdjacentEdgesBoolean
  Rho <- DatObj$Rho
  WeightsInd <- DatObj$WeightsInd
  XThetaInd <- DatObj$XThetaInd
  ZDelta <- DatObj$ZDelta

  ###Set hyperparameter objects
  AAlpha <- HyPara$AAlpha
  BAlpha <- HyPara$BAlpha

  ###Which parameters are user defined?
  UserStarters <- names(Starting)

  ###Set initial value of Sigma
  if ("Sigma" %in% UserStarters) Sigma <- matrix(Starting$Sigma, nrow = 5, ncol = 5)
  if (!("Sigma" %in% UserStarters)) Sigma <- diag(5)

  ###Set initial values of Alpha
  if ("Alpha" %in% UserStarters) Alpha <- Starting$Alpha
  if ((!"Alpha" %in% UserStarters)) Alpha <- AAlpha + BAlpha / 2

  ###Set initial values of Delta
  if ("Delta" %in% UserStarters) Delta <- matrix(Starting$Delta, ncol = 1)
  if ((!"Delta" %in% UserStarters)) Delta <- matrix(c(0, 0, 0, 0, 0), ncol = 1)

  ###Set inital value of random effects
  Beta0 <- matrix(0, nrow = M, ncol = 1)
  Beta1 <- matrix(0, nrow = M, ncol = 1)
  Lambda0 <- matrix(0, nrow = M, ncol = 1)
  Lambda1 <- matrix(0, nrow = M, ncol = 1)
  Eta <- matrix(0, nrow = M, ncol = 1)
  Phi <- CreatePhi_lmc(Beta0, Beta1, Lambda0, Lambda1, Eta, M)

  ###Initialize likelihood covariance objects
  Gamma <- GammaInv <- matrix(0, nrow = M * 5, ncol = M * 5)
  for (p in 1:5) {
    Indeces <- as.matrix(expand.grid(1:M + (p - 1) * M, 1:M + (p - 1) * M))
    WAlpha_p <- WAlphaFnc(Alpha[p], DMLong, AdjacentEdgesBoolean, W, M, WeightsInd)
    GammaInv[Indeces] <- QInv_p <- QInvFnc(WAlpha_p, EyeM, Rho, M)
    Gamma[Indeces] <- CholInv(QInv_p)
  }
  SigmaInv <- CholInv(Sigma)
  A <- chol(Sigma)
  PhiCov <- kronecker(A, EyeM) %*% GammaInv %*% kronecker(t(A), EyeM)
  PhiPrec <- kronecker(solve(t(A)), EyeM) %*% Gamma %*% kronecker(solve(A), EyeM)
  PhiMean <- ZDelta %*% Delta

  ###Get Theta objects
  Theta <- exp(Eta)
  XTheta <- GetXTheta_lmc(Theta, XThetaInd, TimeVec, OneNu, OneN, tNu, N, M)

  ###Joint likelihood objects
  Mu <- kronecker(OneNu, Beta0) + XTheta %*% Beta1
  Sigma2 <- exp(2 * (kronecker(OneNu, Lambda0) + XTheta %*% Lambda1))
  Omega <- diag(as.numeric(Sigma2))
  OmegaInv <- diag(as.numeric(1 / Sigma2))

  ###Save parameter objects
  Para <- list()
  Para$Beta0 <- Beta0
  Para$Beta1 <- Beta1
  Para$Lambda0 <- Lambda0
  Para$Lambda1 <- Lambda1
  Para$Eta <- Eta
  Para$Delta <- Delta
  Para$Alpha <- Alpha
  Para$Sigma <- Sigma
  Para$Sigma2 <- Sigma2
  Para$Omega <- Omega
  Para$OmegaInv <- OmegaInv
  Para$GammaInv <- GammaInv
  Para$Gamma <- Gamma
  Para$SigmaInv <- SigmaInv
  Para$A <- A
  Para$Theta <- Theta
  Para$XTheta <- XTheta
  Para$Mu <- Mu
  Para$Phi <- Phi
  Para$PhiMean <- PhiMean
  Para$PhiPrec <- PhiPrec
  Para$PhiCov <- PhiCov
  return(Para)

}



###Function that creates the data augmentation (i.e. Tobit) booleans------------------------------------------------------------
CreateDatAug_lmc <- function(DatObj) {

  ###Set data object
  YObserved <- DatObj$YObserved
  FamilyInd <- DatObj$FamilyInd
  YStarWide <- DatObj$YStarWide
  M <- DatObj$M
  Nu <- DatObj$Nu

  ###Initialize Data Augmentation Object
  DatAug <- NULL

  ###Normal objects
  if (FamilyInd == 0) {
    DatAug$NBelow <- 0
    DatAug$NAbove <- 0
    DatAug$WhichBelow <- 0
    DatAug$WhichAbove <- 0
    # DatAug$TobitIndeces <- matrix(c(0,0,0,0),2,2)
    # DatAug$ProbitIndeces <- matrix(c(0,0,0,0),2,2)
  }

  ###Probit objects
  if (FamilyInd == 1) {
    TobitBoolean <- YObserved <= 0
    WhichBelow <- which(TobitBoolean)
    NBelow <- length(WhichBelow)
    TobitBooleanMat <- matrix(TobitBoolean, nrow = M, ncol = Nu)
    YStarBelow <- list()
    for (i in 1 : Nu) YStarBelow[[i]] <- YStarWide[!TobitBooleanMat[,i], i]
    NBelowList <- unlist(lapply(YStarBelow, f<-function(x) M - length(x)))
    TobitIndeces <- which(TobitBooleanMat, arr.ind = TRUE)
    TobitIndeces <- TobitIndeces - 1
    ProbitBoolean <- YObserved > 0
    WhichAbove <- which(ProbitBoolean)
    NAbove <- length(WhichAbove)
    ProbitBooleanMat <- matrix(ProbitBoolean, nrow = M, ncol = Nu)
    YStarAbove <- list()
    for (i in 1 : Nu) YStarAbove[[i]] <- YStarWide[!ProbitBooleanMat[,i], i]
    NAboveList <- unlist(lapply(YStarAbove, f<-function(x) M - length(x)))
    ProbitIndeces <- which(ProbitBooleanMat, arr.ind = TRUE)
    ProbitIndeces <- ProbitIndeces - 1

    ###Save objects
    DatAug <- list()
    DatAug$WhichBelow <- WhichBelow - 1
    DatAug$WhichAbove <- WhichAbove - 1
    DatAug$NBelow <- NBelow
    DatAug$NAbove <- NAbove
    DatAug$TobitIndeces <- TobitIndeces
    DatAug$ProbitIndeces <- ProbitIndeces
  }

  ###Tobit objects
  if (FamilyInd == 2) {
    TobitBoolean <- YObserved <= 0
    WhichBelow <- which(TobitBoolean)
    NBelow <- length(WhichBelow)
    TobitBooleanMat <- matrix(TobitBoolean, nrow = M, ncol = Nu)
    YStarNonZero <- list()
    for (i in 1 : Nu) YStarNonZero[[i]] <- YStarWide[!TobitBooleanMat[,i], i]
    NBelowCount <- unlist(lapply(YStarNonZero, f<-function(x) M - length(x)))
    TobitIndeces <- which(TobitBooleanMat, arr.ind = TRUE) - 1
    # TobitIndeces <- TobitIndeces - 1
    # ZDatAug <- model.matrix(~-1 + as.factor(TobitIndeces[,2]))
    # attributes(ZDatAug) <- NULL
    # ZDatAug <- structure(ZDatAug, class = "matrix", dim = c(NBelow, Nu))
    # WDatAug <- array(FALSE, dim = c(M, M, Nu))
    # for (i in 1:NBelow) {
    #   Visit <- TobitIndeces[i, 2] + 1
    #   Location <- TobitIndeces[i, 1] + 1
    #   WDatAug[ , Location, Visit] <- rep(TRUE, M)
    # }
    # WDatAug <- matrix(which(WDatAug) - 1, ncol = 1)

    ###Save objects
    DatAug <- list()
    DatAug$WhichBelow <- WhichBelow - 1
    DatAug$NBelow <- NBelow
    DatAug$TobitBooleanMat <- TobitBooleanMat
    DatAug$YStarNonZero <- YStarNonZero
    DatAug$NBelowCount <- NBelowCount
    DatAug$TobitIndeces <- TobitIndeces
    DatAug$ProbitIndeces <- matrix(c(0,0,0,0),2,2)
    DatAug$NAbove <- 0
    DatAug$WhichAbove <- 0
    # DatAug$ZDatAug <- ZDatAug
    # DatAug$WDatAug <- WDatAug
  }
  return(DatAug)

}



###Function that creates inputs for MCMC sampler--------------------------------------------------------------------------------
CreateMcmc_lmc <- function(MCMC, DatObj) {

  ###Which parameters are user defined?
  UserMCMC <- names(MCMC)

  ###Set MCMC objects
  if ("NBurn" %in% UserMCMC) NBurn <- MCMC$NBurn
  if (!("NBurn" %in% UserMCMC)) NBurn <- 10000
  if ("NSims" %in% UserMCMC) NSims <- MCMC$NSims
  if (!("NSims" %in% UserMCMC)) NSims <- 100000
  if ("NThin" %in% UserMCMC) NThin <- MCMC$NThin
  if (!("NThin" %in% UserMCMC)) NThin <- 10
  if ("NPilot" %in% UserMCMC) NPilot <- MCMC$NPilot
  if (!("NPilot" %in% UserMCMC)) NPilot <- 20

  ###One last check of MCMC user inputs
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!(is.wholenumber(NSims / NThin))) stop('MCMC: "NThin" must be a factor of "NSims"')
  if (!(is.wholenumber(NBurn / NPilot))) stop('MCMC: "NPilot" must be a factor of "NBurn"')

  ###Create MCMC objects
  NTotal <- NBurn + NSims
  WhichKeep <- NBurn + (1:(NSims / NThin)) * NThin
  NKeep <- length(WhichKeep)

  ###Pilot adaptation objects
  WhichPilotAdapt <- (1:NPilot) * NBurn / NPilot
  PilotAdaptDenominator <- WhichPilotAdapt[1]

  ###Burn-in progres bar
  BarLength <- 50 #Burn-in bar length (arbitrary)
  BurnInProgress <- seq(1 / BarLength, 1, 1 / BarLength)
  WhichBurnInProgress <- sapply(BurnInProgress, function(x) tail(which(1 : NBurn <= x * NBurn), 1))

  ###Progress output objects
  SamplerProgress <- seq(0.1, 1.0, 0.1) #Intervals of progress update (arbitrary)
  WhichSamplerProgress <- sapply(SamplerProgress, function(x) tail(which(1:NSims <= x * NSims), 1)) + NBurn
  WhichBurnInProgressInt <- sapply(SamplerProgress, function(x) tail(which(1:NBurn <= x * NBurn), 1))

  ###Save objects
  MCMC <- list()
  MCMC$NBurn <- NBurn
  MCMC$NSims <- NSims
  MCMC$NThin <- NThin
  MCMC$NPilot <- NPilot
  MCMC$NTotal <- NTotal
  MCMC$WhichKeep <- WhichKeep
  MCMC$NKeep <- NKeep
  MCMC$WhichPilotAdapt <- WhichPilotAdapt
  MCMC$PilotAdaptDenominator <- PilotAdaptDenominator
  MCMC$BurnInProgress <- BurnInProgress
  MCMC$WhichBurnInProgress <- WhichBurnInProgress
  MCMC$WhichBurnInProgressInt <- WhichBurnInProgressInt
  MCMC$BarLength <- BarLength
  MCMC$WhichSamplerProgress <- WhichSamplerProgress
  return(MCMC)

}



###Function that creates a storage object for raw samples-----------------------------------------------------------------------
CreateStorage_lmc <- function(DatObj, McmcObj) {

  ###Set data objects
  M <- DatObj$M

  ###Set MCMC objects
  NKeep <- McmcObj$NKeep

  ###Create storage object
  Out <- matrix(nrow = (5 + 15 + 5 + 5 * M), ncol = NKeep)
  return(Out)

}
