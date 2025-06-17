
#*######################################################
#*######################################################
#*######################################################
#*##########      BIOCLOCK: Preprocessing      #########
#*######################################################
#*######################################################
#*######################################################

###################################
####           Set-up           ###
###################################

library(minfi)
library(ENmix)
library(wateRmelon)
library(EpiDISH)
library(dplyr)
library(tidyr)
library(arrow)
library(readxl)
library(writexl)

setwd(paste0("/", file.path("data", "user_homes", "mennovd", "BIOKLOK")))

# In order to calculate the PC clocks, the PC-Clocks GitHub repository is required.
# Clone it from: https://github.com/MorganLevineLab/PC-Clocks/tree/main
# Put the folder in Scripts/Additional and name the folder "PC_Clocks".

# install EPICv2 manifest and annotation from GitHub using Bioconductor (not available on conda)
BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest", release = "0.99.1", force = FALSE, update = FALSE)
BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38", release = "0.99.0", force = FALSE, update = FALSE)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

dir.create("Data/Objects", showWarnings = FALSE, recursive = TRUE) # stores saved R objects
dir.create("Data/Pheno", showWarnings = FALSE, recursive = TRUE) # stores phenotypic data (sample data)
dir.create("Data/Infinium", showWarnings = FALSE, recursive = TRUE) # stores Infinium EPICv2 array idat files
dir.create("Results/PrePro", showWarnings = FALSE, recursive = TRUE) # stores results of preprocessing pipelines
dir.create("Images/QC", showWarnings = FALSE, recursive = TRUE) # stores general quality control images
dir.create("Images/PrePro", showWarnings = FALSE, recursive = TRUE) # stores quality control images of preprocessing pipelines

###################################
####          Functions         ###
###################################

# PC clock script is retrieved from: https://github.com/MorganLevineLab/PC-Clocks
source("Scripts/Additional/run_calcPCClocks.R")

# Preprocessing pipeline script is retrieved from: https://github.com/anilpsori/_pipelines_and_biomarkers/tree/main
source("Scripts/Additional/function_runPipelines.R")

# Function to summarize duplicated CpGs for EPICv2 (retrieved from ENmix: https://github.com/xuz1/ENmix)
rm.cgsuffix <- function(datMeth) {
    cgid <- sapply(strsplit(rownames(datMeth), split = "_"), unlist)[1, ]
    dupcg <- unique(cgid[duplicated(cgid)])
    datMeth2 <- datMeth[cgid %in% dupcg, ]
    cid <- sapply(strsplit(rownames(datMeth2), split = "_"), unlist)[1, ]
    datMeth2 <- aggregate(datMeth2, by = list(cid), FUN = function(x) mean(x, na.rm = TRUE))
    rownames(datMeth2) <- datMeth2[, 1]
    datMeth2 <- as.matrix(datMeth2[, -1])
    datMeth <- datMeth[!(cgid %in% dupcg), ]
    rownames(datMeth) <- sapply(strsplit(rownames(datMeth), split = "_"), unlist)[1, ]
    rbind(datMeth, datMeth2)
}

###################################
####     Data Loading & QC      ###
###################################

# Load methylation data
RGset <- read.metharray.exp(base = file.path("Data", "Infinium"), recursive = TRUE, extended = TRUE, verbose = TRUE) # need extended for QC
RGset@annotation <- c(array = "IlluminaHumanMethylationEPIC2", annotation = "20a1.hg38")
RGset <- RGset[sort(rownames(RGset), index.return = TRUE)$ix, sort(colnames(RGset), index.return = TRUE)$ix]
save(RGset, file = "Data/Objects/RGset_raw.Rdata")

# ENmix Quality control
plotCtrl(RGset)
QC <- QCinfo(RGset)
save(QC, file = "Data/Objects/QC.Rdata")
# 0  samples with percentage of low quality CpG value greater than  0.05  or bisulfite intensity less than  7049.815 
# 5613  CpGs with percentage of low quality value greater than  0.05
# 0  samples are outliers based on averaged total intensity value 
# 0  samples are outliers in beta value distribution 
# 0  outlier samples were added into badsample list

QC$badsample
QC$badCpG
QC$outlier_sample

# Check betas density
Mset_raw <- getmeth(RGset)
save(Mset_raw, file = "Data/Objects/Mset_raw.Rdata")

betas_raw <- getB(Mset_raw)
save(betas_raw, file = "Data/Objects/betas_raw.Rdata")

anno <- getProbeType(Mset_raw)
betas_raw_1 <- betas_raw[anno == "I",]
betas_raw_2 <- betas_raw[anno == "II",]
jpeg("Images/QC/freqpolygon_beta_raw_probe_types.jpg", height = 900, width = 600)
par(mfrow = c(3, 1))
multifreqpoly(betas_raw, main = "Multifreqpoly: All betas raw", xlab = "Beta value")
multifreqpoly(betas_raw_1, main = "Multifreqpoly: Infinium I", xlab = "Beta value")
multifreqpoly(betas_raw_2, main = "Multifreqpoly: Infinium II", xlab = "Beta value")
dev.off()

# Load phenotypic sample data
pheno <- read_xlsx("Data/Pheno/Sample Data.xlsx", sheet = "PhenotypicData")
pheno <- as.data.frame(pheno)
rownames(pheno) <- pheno$array_id
pheno <- pheno[colnames(RGset), ]
save(pheno, file = "Data/Objects/pheno.Rdata")

###################################
####   Preprocessing Pipelines  ###
###################################

# Run all preprocessing pipelines from Ori et al. (DOI: 10.1186/s13059-022-02793-w)

for (method in 1:101) {

    pipeline <- method_to_pipeline(method)
    pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")
    pipeline_image_dir <- paste0("Images/PrePro/", pipeline, "/")
    dir.create(pipeline_image_dir, recursive = TRUE, showWarnings = FALSE)

    if (method %in% c(89, 92, 100)) next # not compatible with EPIC2
    
    if(!(file.exists(paste0(pipeline_data_dir, "betas_no_suffix.feather")))) {
        output <- runPipelines(RGset, method = method, returnAll = FALSE, removeCpG = QC$badCpG)
        betas <- output$probes
        dir.create(pipeline_data_dir)
        arrow::write_feather(betas, paste0(pipeline_data_dir, "betas.feather"))

        # Remove probe ID column
        CpGs <- betas$ProbeID
        betas$ProbeID <- NULL
        betas <- as.matrix(betas)
        rownames(betas) <- CpGs

        # Check betas density plot
        load("Data/Objects/Mset_raw.Rdata")

        anno <- getAnnotation(Mset)
        anno <- data.frame(ProbeID = rownames(anno), Type = anno$Type)
        anno <- anno[match(anno$ProbeID, rownames(betas)),]

        betas_1 <- betas[anno$Type == "I",]
        betas_2 <- betas[anno$Type == "II",]
        jpeg(paste0(pipeline_image_dir, pipeline,".jpg"), height = 900, width = 600)
        par(mfrow = c(3, 1))
        multifreqpoly(betas, main = paste("Multifreqpoly: All betas", pipeline), xlab = "Beta value", legend = NULL)
        multifreqpoly(betas_1, main = "Multifreqpoly: Infinium I", xlab = "Beta value", legend = NULL)
        multifreqpoly(betas_2, main = "Multifreqpoly: Infinium II", xlab = "Beta value", legend = NULL)
        dev.off()

        # Summarize replicate probes (remove probe suffixes) and save betas
        betas <- rm.cgsuffix(betas)
        betas <- as.data.frame(betas)
        betas$id <- rownames(betas)
        arrow::write_feather(betas, paste0(pipeline_data_dir, "betas_no_suffix.feather"))
    }
}

###################################
####      EpiAge Prediction     ###
###################################

# Estimate Horvath & GrimAge epigenetic ages using the Biolearn & PC-Clocks packages

# Specify path to mamba/conda python environment and the Additional/Biolearn_EpiAge.py script
# ! EDIT THIS PATH TO YOUR OWN PATH !
python <-  "/data/user_homes/mennovd/.conda/envs/Bioclock_py/bin/python"
script <- "/data/user_homes/mennovd/BIOKLOK/Scripts/Additional/Biolearn_EpiAge.py"
# see python script: Additional/Biolearn_EpiAge.py

for (method in 1:101) {

    if (method %in% c(89, 92, 100)) next

    pipeline <- method_to_pipeline(method)
    cat("Calculating Epigenetic Clocks:", method, pipeline, "\n")

    pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")

    if (!file.exists(paste0(pipeline_data_dir, "mAge.feather"))) {

        # Biolearn clocks
        load("Data/pheno.Rdata")
        mAge <- pheno[,c("array_id", "subject_id", "exercise_timepoint", "sex", "age", "has_replicate")]
        mAge$female <- as.numeric(mAge$sex == "female")
        mAge$sex <- as.numeric(mAge$sex == "male") + 1
        arrow::write_feather(mAge, paste0(pipeline_data_dir, "mAge.feather"))

        # see python script: Additional/Biolearn_EpiAge.py
        system2(python, args = c(script, pipeline), wait = TRUE)

        biolearn <- arrow::read_feather(paste0(pipeline_data_dir, "biolearn_mAge.feather"))
        biolearn <- biolearn[,c("Horvathv1", "Horvathv2", "DNAmGrimAge_V1", "DNAmGrimAge_V2")]
        colnames(biolearn) <- c("Horvath1", "Horvath2", "GrimAge1", "GrimAge2")
        mAge <- cbind(mAge, biolearn)

        # PC clocks
        betas <- arrow::read_feather(paste0(pipeline_data_dir, "betas_no_suffix.feather"))
        CpGs <- betas$id
        betas$id <- NULL
        betas <- as.matrix(betas)
        rownames(betas) <- CpGs

        PC_mAge <- calcPCClocks(path_to_PCClocks_directory = "Scripts/Additional/PC_Clocks/", datMeth = t(betas), datPheno = mAge[, c("age", "female")])
        PC_mAge <- PC_mAge[,c("PCHorvath1", "PCHorvath2", "PCGrimAge")]
        save(PC_mAge, file = paste0(pipeline_data_dir, "PC_mAge.Rdata"))
        
        mAge <- cbind(mAge, PC_mAge)
        arrow::write_feather(mAge, paste0(pipeline_data_dir, "mAge.feather"))
    }
}

###################################
####    Leukocyte Prediction    ###
###################################

# Estimate blood leukocyte fractions using the EpiDISH RCP method

# Load the reference matrix for the 12 leukocyte fractions on the EPIC array
data(cent12CT.m)

for (method in 1:101) {

    if (method %in% c(89, 92, 100)) next

    pipeline <- method_to_pipeline(method)
    cat("Calculating Blood Cell Composition:", method, pipeline, "\n")

    pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")

    if (!file.exists(paste0(pipeline_data_dir, "bloodcells.feather"))) {

        betas <- arrow::read_feather(paste0(pipeline_data_dir, "betas_no_suffix.feather"))
        CpGs <- betas$id
        betas$id <- NULL
        betas <- as.matrix(betas)
        rownames(betas) <- CpGs

        # Intersect betas and reference matrix CpGs
        cpgs <- intersect(rownames(betas), rownames(cent12CT.m))
        betas_subset <- betas[cpgs,]
        ref_subset <- as.matrix(cent12CT.m)[cpgs,]

        # Estimate leukocyte proportions
        bloodFractions <- epidish(beta.m = betas_subset, ref.m = ref_subset, method = "RPC")$estF

        # Ensure that the proportions sum to one and are percentage
        bloodFractions <- (bloodFractions / rowSums(bloodFractions))*100

        mAge <- arrow::read_feather(paste0(pipeline_data_dir, "mAge.feather"))
        mAge <- cbind(mAge, bloodFractions)
        arrow::write_feather(mAge, paste0(pipeline_data_dir, "mAge.feather"))
        save(bloodFractions, file = paste0(pipeline_data_dir, "bloodcells.Rdata"))
    }
}

########################################
####    Preprocessing optimization   ###
########################################

bloodcells <- c(
    "CD4Tnv",
    "CD4Tmem",
    "CD8Tmem",
    "CD8Tnv",
    "Treg",
    "Bmem",
    "Bnv",
    "Baso",
    "Eos",
    "NK",
    "Neu",
    "Mono"
)

clocks <- c(
    "Horvath1", 
    "Horvath2", 
    "PCHorvath1", 
    "PCHorvath2", 
    "GrimAge1", 
    "GrimAge2",
    "PCGrimAge"
)

predictors <- c(clocks, bloodcells)

# Calculate the biological and technical variance for each predictor and each preprocessing pipeline
TechVar <- data.frame(row.names = predictors)
BioVar <- data.frame(row.names = predictors)

for (method in 1:101) {

    if (method %in% c(89, 92, 100)) next

    # Load the predictor data
    pipeline <- method_to_pipeline(method)
    cat("Calculating measurement noise:", method, pipeline, "\n")
    pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")
    mAge <- arrow::read_feather(paste0(pipeline_data_dir, "mAge.feather"))

    # Calculate the technical variance
    techvar <- mAge[mAge$has_replicate, c("subject_id", "exercise_timepoint", predictors)] %>%
            group_by(subject_id, exercise_timepoint) %>%
            summarize_at(predictors, function(x) var(x)) %>%
            ungroup() %>%
            summarize_at(predictors, function(x) mean(x)) %>%
            as.matrix() %>%
            as.vector()

    TechVar[,pipeline] <- techvar

    # Calculate the biological variance
    biovar <- mAge[, c("subject_id", "exercise_timepoint", predictors)] %>%
            group_by(subject_id, exercise_timepoint) %>%
            summarize_at(predictors, function(x) mean(x)) %>%
            ungroup() %>%
            summarize_at(predictors, function(x) var(x)) %>%
            as.matrix() %>%
            as.vector()

    BioVar[,pipeline] <- biovar
}

save(TechVar, file = "Data/Objects/TechVar.Rdata")
save(BioVar, file = "Data/Objects/BioVar.Rdata")

# Calculate the ratio of biological variance over technical variance (for each predictor and preprocessing pipeline)
BT <- BioVar/(TechVar + 1e-10)
BT_max <- rowMaxs(as.matrix(BT), na.rm = TRUE)
names(BT_max) <- rownames(BT)
optimal_pipelines <- sapply(predictors, function(predictor) colnames(BT)[BT[predictor,] == BT_max[predictor]])
BT_pipelines <- data.frame(predictors = predictors, pipeline = optimal_pipelines)
BT_pipelines[predictors, "Bio/Tech"] <- BT_max[predictors]
BT_pipelines[predictors,]
#            predictors                     pipeline     Bio/Tech
# Horvath1     Horvath1             watermelon_nanet 2.850053e+01
# Horvath2     Horvath2       enmix_est_nodye_q3_rcp 1.056976e+02
# PCHorvath1 PCHorvath1    enmix_est_mean_nonorm_rcp 3.900147e+02
# PCHorvath2 PCHorvath2       enmix_oob_nodye_q2_rcp 3.707282e+02
# GrimAge1     GrimAge1   minfi_raw_quantile_nostrat 4.713229e+02
# GrimAge2     GrimAge2 minfi_funnorm_nobg_nodyecorr 9.156775e+01
# PCGrimAge   PCGrimAge     enmix_neg_nodye_q2_norcp 2.329356e+03
# CD4Tnv         CD4Tnv          watermelon_raw_BMIQ 4.369612e+01
# CD4Tmem       CD4Tmem         minfi_noob_nodyecorr 3.052212e+01
# CD8Tmem       CD8Tmem       enmix_oob_relic_q1_rcp 2.331840e+02
# CD8Tnv         CD8Tnv       enmix_neg_relic_q3_rcp 2.389515e+01
# Treg             Treg    enmix_neg_mean_nonorm_rcp 2.290357e+01
# Bmem             Bmem   enmix_est_nodye_nonorm_rcp 2.065470e+02
# Bnv               Bnv       enmix_est_nodye_q3_rcp 1.324842e+01
# Baso             Baso       enmix_est_nodye_q3_rcp 4.933858e+00
# Eos               Eos enmix_est_relic_nonorm_norcp 3.960289e+05
# NK                 NK             watermelon_danet 7.473183e+01
# Neu               Neu             watermelon_danes 8.381748e+01
# Mono             Mono             watermelon_nasen 5.656306e+01

# For the leukocyte fractions, calculate the mean z-scaled ratio of biological variance over technical variance
BT_BC <- BT[bloodcells,]
BT_BC_scaled <- t(scale(t(BT_BC), center = TRUE, scale = TRUE))
BT_BC_mean <- colMeans(as.matrix(BT_BC_scaled), na.rm = TRUE)
pipeline_bloodcells <- names(which.max(BT_BC_mean))
pipeline_bloodcells
BT_BC_mean[pipeline_bloodcells]
# enmix_est_nodye_q3_rcp 
#               1.002638

# Optimal pipelines:
# PCHorvath1 clock with the enmix_est_mean_nonorm_rcp pipeline
# PCGrimAge clock with the enmix_neg_nodye_q2_norcp pipeline
# enmix_est_nodye_q3_rcp is optimal for the leukocyte deconvolution

# Retrieve the optimal estimates for each predictor
load("Data/Objects/pheno.Rdata")

pipeline <- "enmix_est_mean_nonorm_rcp"
pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")
mAge <- arrow::read_feather(paste0(pipeline_data_dir, "mAge.feather"))

pheno$epigenetic_age_Horvath <- mAge$PCHorvath1

pipeline <- "enmix_neg_nodye_q2_norcp"
pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")
mAge <- arrow::read_feather(paste0(pipeline_data_dir, "mAge.feather"))

pheno$epigenetic_age_GrimAge <- mAge$PCGrimAge

pipeline <- "enmix_est_nodye_q3_rcp"
pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")
mAge <- arrow::read_feather(paste0(pipeline_data_dir, "mAge.feather"))

pheno$B_cells_naive <- mAge$Bnv
pheno$B_cells_memory <- mAge$Bmem
pheno$CD4T_cells_naive <- mAge$CD4Tnv
pheno$CD4T_cells_memory <- mAge$CD4Tmem
pheno$CD8T_cells_naive <- mAge$CD8Tnv
pheno$CD8T_cells_memory <- mAge$CD8Tmem
pheno$T_regulatory_cells <- mAge$Treg
pheno$natural_killer_cells <- mAge$NK
pheno$eosinophils <- mAge$Eos
pheno$basophils <- mAge$Baso
pheno$monocytes <- mAge$Mono
pheno$neutrophils <- mAge$Neu

# Save optimal estimates
save(pheno, file = "Data/Objects/pheno.Rdata")

###################################
########    Statistics    #########
###################################

# Calculate standard error of the mean on the measurement noise in months

pheno[pheno$has_replicate == TRUE, c("subject_id", "exercise_timepoint", "epigenetic_age_Horvath", "epigenetic_age_GrimAge")] %>%
    group_by(subject_id, exercise_timepoint) %>%
    summarize_at(c("epigenetic_age_Horvath", "epigenetic_age_GrimAge"), function(x) abs(diff(x))) %>%
    ungroup() %>%
    summarize(
        mean_noise_Horvath = mean(epigenetic_age_Horvath * 12),
        SEM_noise_Horvath = sd(epigenetic_age_Horvath * 12) / sqrt(n()),
        mean_noise_GrimAge = mean(epigenetic_age_GrimAge * 12),
        SEM_noise_GrimAge = sd(epigenetic_age_GrimAge * 12) / sqrt(n())
    )

#   mean_noise_Horvath SEM_noise_Horvath mean_noise_GrimAge SEM_noise_GrimAge
#                 5.05              1.36               1.91             0.500
