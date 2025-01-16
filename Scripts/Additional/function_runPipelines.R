library(minfi)
library(ENmix)
library(wateRmelon)

########################################################################
########################################################################

# The minfi preprocessIllumina function does not work on EPICv2 due to a small incompatibility
# I have corrected the incompatiility in the following functions
# functions adapted from: https://github.com/hansenlab/minfi

.isRGOrStop <- function(object) {
    if (!is(object, "RGChannelSet")) {
        stop("object is of class '", class(object), "', but needs to be of ",
             "class 'RGChannelSet' or 'RGChannelSetExtended'")
    }
}

.is27k <- function(object) {
    annotation(object)["array"] == "IlluminaHumanMethylation27k"
}

.is450k <- function(object) {
    annotation(object)["array"] == "IlluminaHumanMethylation450k"
}

.isEPIC <- function(object) {
    annotation(object)["array"] == "IlluminaHumanMethylationEPIC"
}

.isEPICv2 <- function(object) {
    annotation(object)["array"] == "IlluminaHumanMethylationEPICv2"
}

#get .getManifestString() function from utils.R
.getManifestString <- function(annotation) {
    if (length(annotation) == 1) {
        if (annotation == "Unknown") {
            stop("Cannot get Manifest object for an 'Unknown' array")
        }
        return(paste0(annotation, "manifest"))
    }
    if ("array" %in% names(annotation)) {
        if (annotation["array"] == "Unknown") {
            stop("Cannot get Manifest object for an 'Unknown' array")
        }
        return(paste0(annotation["array"], "manifest"))
    }
    stop("unable to get the manifest string for this object")
}

#define pmax2() and pmin2() functions
pmax2 <- function(x, y) {
    pmax(x, y)
}
pmin2 <- function(x, y) {
    pmin(x, y)
}

normalize.illumina.control <- function(rgSet, reference = 1) {
    # This function returns an rgset, not a methylset
    # code duplication
    Green <- getGreen(rgSet)
    Red <- getRed(rgSet)

    if (.is450k(rgSet) || .isEPIC(rgSet) || .isEPICv2(rgSet)) {
        AT.controls <- getControlAddress(
            object = rgSet,
            controlType = c("NORM_A", "NORM_T"))
        CG.controls <- getControlAddress(
            object = rgSet,
            controlType = c("NORM_C", "NORM_G"))
    }
    if (.is27k(rgSet)) {
        AT.controls <- getControlAddress(
            object = rgSet,
            controlType = "Normalization-Red")
        CG.controls <- getControlAddress(
            object = rgSet,
            controlType = "Normalization-Green")
    }

    Green.avg <- colMeans2(Green, rows = match(CG.controls, rownames(Green)))
    Red.avg <- colMeans2(Red, rows = match(AT.controls, rownames(Red)))
    ref <- (Green.avg + Red.avg)[reference] / 2
    if (is.na(ref)) {
        stop("perhaps 'reference' refer to an array that is not present.")
    }
    Green.factor <- ref / Green.avg
    Red.factor <- ref / Red.avg
    Green <- sweep(Green, 2, FUN = "*", Green.factor)
    Red <- sweep(Red, 2, FUN = "*", Red.factor)
    assay(rgSet, "Green") <- Green
    assay(rgSet, "Red") <- Red
    rgSet
}

bgcorrect.illumina <- function(rgSet) {
    .isRGOrStop(rgSet)
    Green <- getGreen(rgSet)
    Red <- getRed(rgSet)
    if (.is450k(rgSet) || .isEPIC(rgSet) || .isEPICv2(rgSet)) {
        NegControls <- getControlAddress(rgSet, controlType = "NEGATIVE")
    }
    if (.is27k(rgSet)) {
        NegControls <- getControlAddress(rgSet, controlType = "Negative")
    }
    Green.bg <- apply(Green[NegControls, , drop = FALSE], 2, function(xx) {
        sort(as.vector(xx))[31]
    })
    Red.bg <- apply(Red[NegControls, , drop = FALSE], 2, function(xx) {
        sort(as.vector(xx))[31]
    })
    Green <- pmax2(sweep(Green, 2, Green.bg), 0)
    Red <- pmax2(sweep(Red, 2, Red.bg), 0)
    assay(rgSet, "Green") <- Green
    assay(rgSet, "Red") <- Red
    rgSet
}

# Exported functions -----------------------------------------------------------

# TODO: Document: This does not realize the result for a DelayedMatrix-backed
#       RGChannelSet{Extended}
preprocessIllumina <- function(rgSet, bg.correct = TRUE,
                               normalize = c("controls", "no"), reference = 1) {
    .isRGOrStop(rgSet)
    normalize <- match.arg(normalize)

    if (normalize == "controls") {
        rgSet <- normalize.illumina.control(rgSet, reference = reference)
    }
    if (bg.correct) {
        rgSet <- bgcorrect.illumina(rgSet)
    }
    out <- preprocessRaw(rgSet)
    preprocess <- sprintf(
        "Illumina, bg.correct = %s, normalize = %s, reference = %d",
        bg.correct, normalize, reference)

    # TODO: The manifest package version is currently not updated since
    #       `packageVersion(getManifest(rgSet))` fails. `packageVersion()`
    #       expects a string
    out@preprocessMethod <- c(
        rg.norm = preprocess,
        minfi = as.character(packageVersion("minfi")),
        manifest = as.character(
            packageVersion(.getManifestString(rgSet@annotation))))
    #packageVersion(getManifest(rgSet)))
    out
}

########################################################################
########################################################################

# This function runs different DNAm data processing pipelines from the wateRmelon, minfi, en ENmix package
# Adapted from: Ori, A. P. S., Lu, A. T., Horvath, S. & Ophoff, R. A. Significant variation in the performance of DNA methylation predictors across data preprocessing and normalization strategies. Genome Biol 23, 225 (2022).
# DOI: 10.1186/s13059-022-02793-w
# GitHub: https://github.com/anilpsori/DNAm_pipelines_and_biomarkers/tree/main

runPipelines = function(RGset, method = c(1:101), returnAll = FALSE, selectCpG = NULL, removeCpG = NULL){
  ##RGset: RGChannelSetExtended object of DNAm data
  ##this function can call different DNAm data processing pipelines from the wateRmelon, minfi, en ENmix package
  ##method variable: a numeric single variable or vector that says which pipeline to run
  ##selectCpG: CpGs you would like to subset and write to .csv as output (e.g. those in datMiniAnnotation3.csv)
  ##returnALL: if processed values for ALL probes should be returned and no subsetting be done
  
  ##define methylation data
  M = RGset #this should be an RGChannelSet(Extended) object
  
  #########################
  #### ENmix methods
  #########################
  #apply normalization and background correction
  #bgParaEst: oob, est, neg
  #dyeCorr: mean, RELIC, none
  #normalization: q1,q2,q3, none
  #RCP (dyebias): TRUE or FALSE
  
  if(method == 1){
    pipeline = "enmix_oob_mean_nonorm_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=5)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 2){
    pipeline = "enmix_est_mean_nonorm_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=5)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 3){
    pipeline = "enmix_neg_mean_nonorm_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=5)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 4){
    pipeline = "enmix_oob_relic_nonorm_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=5)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 5){
    pipeline = "enmix_est_relic_nonorm_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=5)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 6){
    pipeline = "enmix_neg_relic_nonorm_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=5)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 7){
    pipeline = "enmix_oob_nodye_nonorm_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=5)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 8){
    pipeline = "enmix_est_nodye_nonorm_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=5)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 9){
    pipeline = "enmix_neg_nodye_nonorm_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=5)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 10){
    pipeline = "enmix_oob_mean_q1_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 11){
    pipeline = "enmix_est_mean_q1_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 12){
    pipeline = "enmix_neg_mean_q1_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 13){
    pipeline = "enmix_oob_relic_q1_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 14){
    pipeline = "enmix_est_relic_q1_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 15){
    pipeline = "enmix_neg_relic_q1_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 16){
    pipeline = "enmix_oob_nodye_q1_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 17){
    pipeline = "enmix_est_nodye_q1_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 18){
    pipeline = "enmix_neg_nodye_q1_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 19){
    pipeline = "enmix_oob_mean_q2_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 20){
    pipeline = "enmix_est_mean_q2_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 21){
    pipeline = "enmix_neg_mean_q2_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 22){
    pipeline = "enmix_oob_relic_q2_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 23){
    pipeline = "enmix_est_relic_q2_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 24){
    pipeline = "enmix_neg_relic_q2_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 25){
    pipeline = "enmix_oob_nodye_q2_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 26){
    pipeline = "enmix_est_nodye_q2_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 27){
    pipeline = "enmix_neg_nodye_q2_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 28){
    pipeline = "enmix_oob_mean_q3_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 29){
    pipeline = "enmix_est_mean_q3_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 30){
    pipeline = "enmix_neg_mean_q3_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 31){
    pipeline = "enmix_oob_relic_q3_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 32){
    pipeline = "enmix_est_relic_q3_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 33){
    pipeline = "enmix_neg_relic_q3_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 34){
    pipeline = "enmix_oob_nodye_q3_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 35){
    pipeline = "enmix_est_nodye_q3_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 36){
    pipeline = "enmix_neg_nodye_q3_norcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getB(normalized)
  }

  else if (method == 37){
    pipeline = "enmix_oob_mean_nonorm_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=5)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 38){
    pipeline = "enmix_est_mean_nonorm_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=5)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 39){
    pipeline = "enmix_neg_mean_nonorm_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=5)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 40){
    pipeline = "enmix_oob_relic_nonorm_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=5)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 41){
    pipeline = "enmix_est_relic_nonorm_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=5)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 42){
    pipeline = "enmix_neg_relic_nonorm_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=5)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 43){
    pipeline = "enmix_oob_nodye_nonorm_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=5)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 44){
    pipeline = "enmix_est_nodye_nonorm_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=5)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 45){
    pipeline = "enmix_neg_nodye_nonorm_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=5)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 46){
    pipeline = "enmix_oob_mean_q1_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 47){
    pipeline = "enmix_est_mean_q1_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 48){
    pipeline = "enmix_neg_mean_q1_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 49){
    pipeline = "enmix_oob_relic_q1_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 50){
    pipeline = "enmix_est_relic_q1_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 51){
    pipeline = "enmix_neg_relic_q1_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 52){
    pipeline = "enmix_oob_nodye_q1_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 53){
    pipeline = "enmix_est_nodye_q1_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 54){
    pipeline = "enmix_neg_nodye_q1_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=5)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 55){
    pipeline = "enmix_oob_mean_q2_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 56){
    pipeline = "enmix_est_mean_q2_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 57){
    pipeline = "enmix_neg_mean_q2_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 58){
    pipeline = "enmix_oob_relic_q2_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 59){
    pipeline = "enmix_est_relic_q2_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 60){
    pipeline = "enmix_neg_relic_q2_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 61){
    pipeline = "enmix_oob_nodye_q2_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 62){
    pipeline = "enmix_est_nodye_q2_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 63){
    pipeline = "enmix_neg_nodye_q2_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=5)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 64){
    pipeline = "enmix_oob_mean_q3_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 65){
    pipeline = "enmix_est_mean_q3_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 66){
    pipeline = "enmix_neg_mean_q3_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 67){
    pipeline = "enmix_oob_relic_q3_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 68){
    pipeline = "enmix_est_relic_q3_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 69){
    pipeline = "enmix_neg_relic_q3_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 70){
    pipeline = "enmix_oob_nodye_q3_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 71){
    pipeline = "enmix_est_nodye_q3_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }

  else if (method == 72){
    pipeline = "enmix_neg_nodye_q3_rcp"
    print("Now running background:  ");cat(method, pipeline, "\n")
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=5)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  #########################
  #### WateRmelon Methods
  #########################
  else if (method == 73){
    pipeline = "watermelon_dasen"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #dasen: same as nasen but type I and type II backgrounds are normalized first. This is our recommended   method
    probecorrected = dasen(M)
    probecorrected = getB(probecorrected) # I ADDED this line, because it doesnt give betas otherwise, maybe needed in the other wateRmelon methods as well
  }

  else if (method == 74){
    pipeline = "watermelon_naten"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #naten: quantile normalizes methylated and unmethylated intensities separately, then calculates betas
    probecorrected = naten(M)
    probecorrected = getB(probecorrected)
  }

  else if (method == 75){
    pipeline = "watermelon_nanet"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #nanet: quantile normalizes methylated and unmethylated intensities together, then calculates betas. This should equalize dye bias.
    probecorrected = nanet(M)
    probecorrected = getB(probecorrected)
  }

  else if (method == 76){
    pipeline = "watermelon_nanes"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #nanes: quantile normalizes methylated and unmethylated intensities separately, except for type II probes where methylated and unmethylated are normalized together
    #This should equalize dye bias without affecting type I probes which are not susceptible
    probecorrected = nanes(M)
    probecorrected = getB(probecorrected)
  }

  else if (method == 77){
    pipeline = "watermelon_danes"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #danes: same as nanes, except typeI and type II background are equalised first
    probecorrected = danes(M)
    probecorrected = getB(probecorrected)
  }

  else if (method == 78){
    pipeline = "watermelon_danet"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #danet: same as nanet, except typeI and type II background are equalised first
    probecorrected = danet(M)
    probecorrected = getB(probecorrected)
  }

  else if (method == 79){
    pipeline = "watermelon_danen"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #danen: background equalisation only, no normalization
    probecorrected = danen(M)
    probecorrected = getB(probecorrected)
  }

  else if (method == 80){
    pipeline = "watermelon_daten1"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #daten1: same as naten, except typeI and type II background are equalised first (smoothed only for methylated)
    probecorrected = daten1(M)
    probecorrected = getB(probecorrected)
  }

  else if (method == 81){
    pipeline = "watermelon_daten2"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #daten2: same as naten, except typeI and type II background are equalised first (smoothed for methylated an unmethylated)
    probecorrected = daten2(M)
    probecorrected = getB(probecorrected)
  }

  else if (method == 82){
    pipeline = "watermelon_nasen"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #nasen: same as naten but typeI and typeII intensities quantile normalized separately
    temp = nasen(M)
    probecorrected = getB(temp)
    #probecorrected = temp$beta
  }

  else if (method == 83){
    pipeline = "watermelon_raw_BMIQ"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #BMIQ: is an intra-sample normalisation procedure, correcting the bias of type-2 probe values
    raw = preprocessRaw(M)
    probecorrected = BMIQ(raw)
  }

  else if (method == 84){
    pipeline = "watermelon_raw_PBC"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #BMIQ: is an intra-sample normalisation procedure, correcting the bias of type-2 probe values
    raw = preprocessRaw(M)
    probecorrected = fuks(raw)
  }
  
  #########################
  #### Minfi methods
  #########################
  else if (method == 85){
    pipeline = "minfi_raw"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #preprocessRaw: Converts the Red/Green channel for an Illumina methylation array into methylation signal, without using any normalization.
    temp = preprocessRaw(M)
    probecorrected = getB(temp)
  }

  else if (method == 86){
    pipeline = "minfi_illumina_bg_nonorm"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #bgcorrect.illumina: implements Illumina GenomeStudio background correction
    temp = preprocessIllumina(M, bg.correct=T,normalize="no")
    probecorrected = getB(temp)
  }

  else if (method == 87){
    pipeline = "minfi_illumina_bg_normcontrol"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #bgcorrect.illumina: implements Illumina GenomeStudio background correction and normalization
    temp = preprocessIllumina(M, bg.correct=T,normalize="controls",reference=1)
    probecorrected = getB(temp)
  }

  else if (method == 88){
    pipeline = "minfi_illumina_nobg_normcontrol"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #bgcorrect.illumina: implements Illumina GenomeStudio background correction and normalization
    temp = preprocessIllumina(M, bg.correct=F,normalize="controls",reference=1)
    probecorrected = getB(temp)
  }

  else if (method == 89){
    pipeline = "minfi_noob_dyecorr"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #preprocessNoob: Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays
    temp = preprocessNoob(M,dyeCorr=T,dyeMethod="single")
    probecorrected = getB(temp)
  } else if (method == 90){
    pipeline = "minfi_noob_nodyecorr"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #preprocessNoob: Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays
    temp = preprocessNoob(M,dyeCorr=F)
    probecorrected = getB(temp)
  } 
  else if (method == 91){
    pipeline = "minfi_funnorm_nobg_nodyecorr"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    ##preprocessFunnorm: a between-array normalization method for the Illumina Infinium HumanMethylation450 platform
    ##It removes unwanted variation by regressing out variability explained by the control probes present on the array.
    temp = preprocessFunnorm(M, bgCorr=F,dyeCorr=F)
    probecorrected = getBeta(temp)
  }

  else if (method == 92){
    pipeline = "minfi_funnorm_bg_dyecorr"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #preprocessNoob: Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays
    ##preprocessFunnorm: a between-array normalization method for the Illumina Infinium HumanMethylation450 platform. It removes unwanted variation by regressing out variability explained by the control probes present on the array.
    temp = preprocessFunnorm(M, bgCorr=T,dyeCorr=T)
    probecorrected = getBeta(temp)
  }

  else if (method == 93){
    pipeline = "minfi_funnorm_bg_nodyecorr"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #preprocessNoob: Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays
    ##preprocessFunnorm: a between-array normalization method for the Illumina Infinium HumanMethylation450 platform. It removes unwanted variation by regressing out variability explained by the control probes present on the array.
    temp = preprocessFunnorm(M, bgCorr=T,dyeCorr=F)
    probecorrected = getBeta(temp)
  }

  else if (method == 94){
    pipeline = "minfi_raw_quantile_strat"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    ##PreprocessQuantile:This function implements stratified quantile normalization preprocessing for Illumina methylation microarrays.
    ##Probes are stratified by region (CpG island, shore, etc.)
    temp = preprocessRaw(M)
    temp2 = preprocessQuantile(temp, quantileNormalize = T,stratified = T)
    probecorrected = getBeta(temp2)
  }

  else if (method == 95){
    pipeline = "minfi_raw_quantile_nostrat"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    ##PreprocessQuantile:This function implements stratified quantile normalization preprocessing for Illumina methylation microarrays.
    ##Probes are stratified by region (CpG island, shore, etc.)
    temp = preprocessRaw(M)
    temp2 = preprocessQuantile(temp, quantileNormalize = T,stratified = F)
    probecorrected = getBeta(temp2)
  }

  else if (method == 96){
    pipeline = "minfi_illumina_bg_quantile_strat"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    ##PreprocessQuantile:This function implements stratified quantile normalization preprocessing for Illumina methylation microarrays. Probes are stratified by region (CpG island, shore, etc.)
    temp = bgcorrect.illumina(M)
    temp2 = preprocessQuantile(temp, quantileNormalize = T,stratified = T)
    probecorrected = getBeta(temp2)
  }

  else if (method == 97){
    pipeline = "minfi_illumina_bg_quantile_nostrat"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    ##preprocessQuantile: This function implements stratified quantile normalization preprocessing for Illumina methylation microarrays. Probes are stratified by region (CpG island, shore, etc.)
    temp = bgcorrect.illumina(M)
    temp2 = preprocessQuantile(temp, quantileNormalize = T,stratified = F)
    probecorrected = getBeta(temp2)
  }

  else if (method == 98){
    pipeline = "minfi_raw_SWAN"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    ##preprocessSWAN: Subset-quantile Within Array Normalisation (SWAN) is a within array normalisation method 
    temp = preprocessSWAN(M)
    probecorrected = getB(temp)
  }

  else if (method == 99){
    pipeline = "minfi_illumina_bg_SWAN"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    ##preprocessSWAN: Subset-quantile Within Array Normalisation (SWAN) is a within array normalisation method 
    temp = bgcorrect.illumina(M)
    temp2 = preprocessSWAN(temp)
    probecorrected = getB(temp2)
  }

  else if (method == 100){
    pipeline = "cross_noob_dyecorr_BMIQ"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #preprocessNoob: Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays
    #BMIQ: is an intra-sample normalisation procedure, correcting the bias of type-2 probe values
    temp = preprocessNoob(M,dyeCorr=T,dyeMethod="single")
    probecorrected = BMIQ(temp)
  }

  else if (method == 101){
    pipeline = "cross_noob_nodyecorr_BMIQ"
    print("Now running method:  ");cat(method, pipeline, "\n")
    
    #preprocessNoob: Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays
    #BMIQ: is an intra-sample normalisation procedure, correcting the bias of type-2 probe values
    temp = preprocessNoob(M,dyeCorr=F,dyeMethod="single")
    probecorrected = BMIQ(temp)
  }

  else {
    stop("Method does not exist")
  }
  
  ##get betas and write 30K probes from datProbes_30491.csv
  dat0 = data.frame(probecorrected)
  probes = data.frame(ProbeID=rownames(dat0))
  dat1=cbind(probes,dat0)
  
  if(returnAll){
    output = list(pipeline=pipeline,probes = dat1)
    return(output)
  }
  else{
    if(!is.null(selectCpG)) {
      dat2=dat1[dat1$ProbeID %in% selectCpG,]
      output = list(pipeline=pipeline,probes = dat2)
      return(output)  
    } else if (!is.null(removeCpG)) {
      dat2=dat1[!(dat1$ProbeID %in% removeCpG),]
      output = list(pipeline=pipeline,probes = dat2)
      return(output)  
    } else {
      warning("No CpGs have been excluded as the selectCpG and removeCpG arguments are both NULL")
      output = list(pipeline=pipeline,probes = dat1)
      return(output)
    }
  }
}

# This function will translate the method number to the pipeline name
method_to_pipeline = function(method) {
    switch(method,
           "enmix_oob_mean_nonorm_norcp",
           "enmix_est_mean_nonorm_norcp",
           "enmix_neg_mean_nonorm_norcp",
           "enmix_oob_relic_nonorm_norcp",
           "enmix_est_relic_nonorm_norcp",
           "enmix_neg_relic_nonorm_norcp",
           "enmix_oob_nodye_nonorm_norcp",
           "enmix_est_nodye_nonorm_norcp",
           "enmix_neg_nodye_nonorm_norcp",
           "enmix_oob_mean_q1_norcp",
           "enmix_est_mean_q1_norcp",
           "enmix_neg_mean_q1_norcp",
           "enmix_oob_relic_q1_norcp",
           "enmix_est_relic_q1_norcp",
           "enmix_neg_relic_q1_norcp",
           "enmix_oob_nodye_q1_norcp",
           "enmix_est_nodye_q1_norcp",
           "enmix_neg_nodye_q1_norcp",
           "enmix_oob_mean_q2_norcp",
           "enmix_est_mean_q2_norcp",
           "enmix_neg_mean_q2_norcp",
           "enmix_oob_relic_q2_norcp",
           "enmix_est_relic_q2_norcp",
           "enmix_neg_relic_q2_norcp",
           "enmix_oob_nodye_q2_norcp",
           "enmix_est_nodye_q2_norcp",
           "enmix_neg_nodye_q2_norcp",
           "enmix_oob_mean_q3_norcp",
           "enmix_est_mean_q3_norcp",
           "enmix_neg_mean_q3_norcp",
           "enmix_oob_relic_q3_norcp",
           "enmix_est_relic_q3_norcp",
           "enmix_neg_relic_q3_norcp",
           "enmix_oob_nodye_q3_norcp",
           "enmix_est_nodye_q3_norcp",
           "enmix_neg_nodye_q3_norcp",
           "enmix_oob_mean_nonorm_rcp",
           "enmix_est_mean_nonorm_rcp",
           "enmix_neg_mean_nonorm_rcp",
           "enmix_oob_relic_nonorm_rcp",
           "enmix_est_relic_nonorm_rcp",
           "enmix_neg_relic_nonorm_rcp",
           "enmix_oob_nodye_nonorm_rcp",
           "enmix_est_nodye_nonorm_rcp",
           "enmix_neg_nodye_nonorm_rcp",
           "enmix_oob_mean_q1_rcp",
           "enmix_est_mean_q1_rcp",
           "enmix_neg_mean_q1_rcp",
           "enmix_oob_relic_q1_rcp",
           "enmix_est_relic_q1_rcp",
           "enmix_neg_relic_q1_rcp",
           "enmix_oob_nodye_q1_rcp",
           "enmix_est_nodye_q1_rcp",
           "enmix_neg_nodye_q1_rcp",
           "enmix_oob_mean_q2_rcp",
           "enmix_est_mean_q2_rcp",
           "enmix_neg_mean_q2_rcp",
           "enmix_oob_relic_q2_rcp",
           "enmix_est_relic_q2_rcp",
           "enmix_neg_relic_q2_rcp",
           "enmix_oob_nodye_q2_rcp",
           "enmix_est_nodye_q2_rcp",
           "enmix_neg_nodye_q2_rcp",
           "enmix_oob_mean_q3_rcp",
           "enmix_est_mean_q3_rcp",
           "enmix_neg_mean_q3_rcp",
           "enmix_oob_relic_q3_rcp",
           "enmix_est_relic_q3_rcp",
           "enmix_neg_relic_q3_rcp",
           "enmix_oob_nodye_q3_rcp",
           "enmix_est_nodye_q3_rcp",
           "enmix_neg_nodye_q3_rcp",
           "watermelon_dasen",
           "watermelon_naten",
           "watermelon_nanet",
           "watermelon_nanes",
           "watermelon_danes",
           "watermelon_danet",
           "watermelon_danen",
           "watermelon_daten1",
           "watermelon_daten2",
           "watermelon_nasen",
           "watermelon_raw_BMIQ",
           "watermelon_raw_PBC",
           "minfi_raw",
           "minfi_illumina_bg_nonorm",
           "minfi_illumina_bg_normcontrol",
           "minfi_illumina_nobg_normcontrol",
           "minfi_noob_dyecorr",
           "minfi_noob_nodyecorr",
           "minfi_funnorm_nobg_nodyecorr",
           "minfi_funnorm_bg_dyecorr",
           "minfi_funnorm_bg_nodyecorr",
           "minfi_raw_quantile_strat",
           "minfi_raw_quantile_nostrat",
           "minfi_illumina_bg_quantile_strat",
           "minfi_illumina_bg_quantile_nostrat",
           "minfi_raw_SWAN",
           "minfi_illumina_bg_SWAN",
           "cross_noob_dyecorr_BMIQ",
           "cross_noob_nodyecorr_BMIQ",
           stop("Invalid method")
    )
}

# This function will translate the pipeline name to the method number
pipeline_to_method <- function(pipeline) {
    method_list <- c(
        "enmix_oob_mean_nonorm_norcp" = 1,
        "enmix_est_mean_nonorm_norcp" = 2,
        "enmix_neg_mean_nonorm_norcp" = 3,
        "enmix_oob_relic_nonorm_norcp" = 4,
        "enmix_est_relic_nonorm_norcp" = 5,
        "enmix_neg_relic_nonorm_norcp" = 6,
        "enmix_oob_nodye_nonorm_norcp" = 7,
        "enmix_est_nodye_nonorm_norcp" = 8,
        "enmix_neg_nodye_nonorm_norcp" = 9,
        "enmix_oob_mean_q1_norcp" = 10,
        "enmix_est_mean_q1_norcp" = 11,
        "enmix_neg_mean_q1_norcp" = 12,
        "enmix_oob_relic_q1_norcp" = 13,
        "enmix_est_relic_q1_norcp" = 14,
        "enmix_neg_relic_q1_norcp" = 15,
        "enmix_oob_nodye_q1_norcp" = 16,
        "enmix_est_nodye_q1_norcp" = 17,
        "enmix_neg_nodye_q1_norcp" = 18,
        "enmix_oob_mean_q2_norcp" = 19,
        "enmix_est_mean_q2_norcp" = 20,
        "enmix_neg_mean_q2_norcp" = 21,
        "enmix_oob_relic_q2_norcp" = 22,
        "enmix_est_relic_q2_norcp" = 23,
        "enmix_neg_relic_q2_norcp" = 24,
        "enmix_oob_nodye_q2_norcp" = 25,
        "enmix_est_nodye_q2_norcp" = 26,
        "enmix_neg_nodye_q2_norcp" = 27,
        "enmix_oob_mean_q3_norcp" = 28,
        "enmix_est_mean_q3_norcp" = 29,
        "enmix_neg_mean_q3_norcp" = 30,
        "enmix_oob_relic_q3_norcp" = 31,
        "enmix_est_relic_q3_norcp" = 32,
        "enmix_neg_relic_q3_norcp" = 33,
        "enmix_oob_nodye_q3_norcp" = 34,
        "enmix_est_nodye_q3_norcp" = 35,
        "enmix_neg_nodye_q3_norcp" = 36,
        "enmix_oob_mean_nonorm_rcp" = 37,
        "enmix_est_mean_nonorm_rcp" = 38,
        "enmix_neg_mean_nonorm_rcp" = 39,
        "enmix_oob_relic_nonorm_rcp" = 40,
        "enmix_est_relic_nonorm_rcp" = 41,
        "enmix_neg_relic_nonorm_rcp" = 42,
        "enmix_oob_nodye_nonorm_rcp" = 43,
        "enmix_est_nodye_nonorm_rcp" = 44,
        "enmix_neg_nodye_nonorm_rcp" = 45,
        "enmix_oob_mean_q1_rcp" = 46,
        "enmix_est_mean_q1_rcp" = 47,
        "enmix_neg_mean_q1_rcp" = 48,
        "enmix_oob_relic_q1_rcp" = 49,
        "enmix_est_relic_q1_rcp" = 50,
        "enmix_neg_relic_q1_rcp" = 51,
        "enmix_oob_nodye_q1_rcp" = 52,
        "enmix_est_nodye_q1_rcp" = 53,
        "enmix_neg_nodye_q1_rcp" = 54,
        "enmix_oob_mean_q2_rcp" = 55,
        "enmix_est_mean_q2_rcp" = 56,
        "enmix_neg_mean_q2_rcp" = 57,
        "enmix_oob_relic_q2_rcp" = 58,
        "enmix_est_relic_q2_rcp" = 59,
        "enmix_neg_relic_q2_rcp" = 60,
        "enmix_oob_nodye_q2_rcp" = 61,
        "enmix_est_nodye_q2_rcp" = 62,
        "enmix_neg_nodye_q2_rcp" = 63,
        "enmix_oob_mean_q3_rcp" = 64,
        "enmix_est_mean_q3_rcp" = 65,
        "enmix_neg_mean_q3_rcp" = 66,
        "enmix_oob_relic_q3_rcp" = 67,
        "enmix_est_relic_q3_rcp" = 68,
        "enmix_neg_relic_q3_rcp" = 69,
        "enmix_oob_nodye_q3_rcp" = 70,
        "enmix_est_nodye_q3_rcp" = 71,
        "enmix_neg_nodye_q3_rcp" = 72,
        "watermelon_dasen" = 73,
        "watermelon_naten" = 74,
        "watermelon_nanet" = 75,
        "watermelon_nanes" = 76,
        "watermelon_danes" = 77,
        "watermelon_danet" = 78,
        "watermelon_danen" = 79,
        "watermelon_daten1" = 80,
        "watermelon_daten2" = 81,
        "watermelon_nasen" = 82,
        "watermelon_raw_BMIQ" = 83,
        "watermelon_raw_PBC" = 84,
        "minfi_raw" = 85,
        "minfi_illumina_bg_nonorm" = 86,
        "minfi_illumina_bg_normcontrol" = 87,
        "minfi_illumina_nobg_normcontrol" = 88,
        "minfi_noob_dyecorr" = 89,
        "minfi_noob_nodyecorr" = 90,
        "minfi_funnorm_nobg_nodyecorr" = 91,
        "minfi_funnorm_bg_dyecorr" = 92,
        "minfi_funnorm_bg_nodyecorr" = 93,
        "minfi_raw_quantile_strat" = 94,
        "minfi_raw_quantile_nostrat" = 95,
        "minfi_illumina_bg_quantile_strat" = 96,
        "minfi_illumina_bg_quantile_nostrat" = 97,
        "minfi_raw_SWAN" = 98,
        "minfi_illumina_bg_SWAN" = 99,
        "cross_noob_dyecorr_BMIQ" = 100,
        "cross_noob_nodyecorr_BMIQ" = 101,
    )
    
    if (pipeline %in% names(method_list)) {
        return(method_list[[pipeline]])
    } else {
        stop("Invalid pipeline")
    }
}
