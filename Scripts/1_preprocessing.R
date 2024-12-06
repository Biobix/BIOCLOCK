
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
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(dplyr)
library(tidyr)
library(feather)
library(arrow)
library(readxl)
library(writexl)

setwd(paste0("/", file.path("data", "user_homes", "mennovd", "BIOKLOK")))

dir.create("Data/Clocks", showWarnings = FALSE, recursive = TRUE)
dir.create("Data/Objects", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/PrePro", showWarnings = FALSE, recursive = TRUE)
dir.create("Images/QC", showWarnings = FALSE, recursive = TRUE)
dir.create("Images/PrePro", showWarnings = FALSE, recursive = TRUE)

###################################
####          Functions         ###
###################################

# PC clock script is available at https://github.com/MorganLevineLab/PC-Clocks
source("Scripts/Additional/run_calcPCClocks.R")

# Preprocessing pipeline script is available at: https://github.com/anilpsori/_pipelines_and_biomarkers/tree/main
source("Scripts/Additional/function_runPipelines.R")

# Function to summarize duplicated CpGs for EPIC2 (retrieved from ENmix: https://github.com/xuz1/ENmix)
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
pheno <- read_xlsx("Data/Pheno/Sample_Data.xlsx")
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

    if (method %in% c(86, 87, 88, 89, 92, 96, 97, 99, 100)) next # not compatible with EPIC2
    
    if(!(file.exists(paste0(pipeline_data_dir, "betas_no_suffix.feather")))) {
        output <- runPipelines(RGset, method = method, returnAll = FALSE, removeCpG = QC$badCpG)
        betas <- output$probes
        dir.create(pipeline_data_dir)
        write_feather(betas, paste0(pipeline_data_dir, "betas.feather"))

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

        # Remove probe suffixes and save betas
        betas <- rm.cgsuffix(betas)
        betas <- as.data.frame(betas)
        betas$id <- rownames(betas)
        write_feather(betas, paste0(pipeline_data_dir, "betas_no_suffix.feather"))
    }
}

###################################
####      EpiAge Prediction     ###
###################################

# Estimate Horvath & GrimAge epigenetic ages using the Biolearn & PC-Clocks packages

for (method in 1:102) {

    if (method %in% c(86, 87, 88, 89, 92, 96, 97, 99, 100)) next

    pipeline <- method_to_pipeline(method)
    cat("Calculating Epigenetic Clocks:", method, pipeline, "\n")

    pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")

    if (!file.exists(paste0(pipeline_data_dir, "mAge.feather"))) {

        # Biolearn clocks
        load("Data/pheno.Rdata")
        mAge <- pheno[,c("array_id", "sex", "age")]
        colnames(mAge) <- c("id", "sex", "age")
        mAge$sex <- (mAge$sex == "male") + 1
        write_feather(mAge, paste0(pipeline_data_dir, "mAge.feather"))

        # see python script: Biolearn_EpiAge.py
        python <-  "/data/user_homes/mennovd/BIOKLOK/.venv/bin/python3.10"
        script <- "/data/user_homes/mennovd/BIOKLOK/Scripts/Additional/Biolearn_EpiAge.py"
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
        
        mAge$female <- as.numeric(mAge$sex == "female")

        PCmAge <- calcPCClocks(path_to_PCClocks_directory = "Clocks/Clock_Metadata/", datMeth = t(betas), datPheno = mAge[, c("age", "female")])
        PCmAge <- PCmAge[,c("PCHorvath1", "PCHorvath2", "PCGrimAge")]
        write_feather(PCmAge, paste0(pipeline_data_dir, "PCmAge.feather"))
        
        mAge <- cbind(mAge, PCmAge)
        write_feather(mAge, paste0(pipeline_data_dir, "mAge.feather"))
    }
}

###################################
####    Leukocyte Prediction    ###
###################################

# Estimate blood leukocyte fractions using the ENmix (estimateCellProp) method
# uses the Houseman et al. (2012) Quadratic Programming method

for (method in 1:102) {

    if (method %in% c(86, 87, 88, 89, 92, 96, 97, 99, 100)) next

    pipeline <- method_to_pipeline(method)
    cat("Calculating Blood Cell Composition:", method, pipeline, "\n")

    pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")

    if (!file.exists(paste0(pipeline_data_dir, "bloodcells.feather"))) {

        mAge <- arrow::read_feather(paste0(pipeline_data_dir, "mAge.feather"))
        betas <- arrow::read_feather(paste0(pipeline_data_dir, "betas_no_suffix.feather"))

        CpGs <- betas$id
        betas$id <- NULL
        betas <- as.matrix(betas)
        rownames(betas) <- CpGs

        # Estimate leukocyte proportions
        cellProp <- estimateCellProp(betas, refdata = "FlowSorted.Blood.EPIC", nonnegative = TRUE, nProbes = 100, normalize = FALSE)
        cellProp <- cellProp[, -1]

        # Ensure that the proportions sum to one and are percentage
        cellProp <- (cellProp / rowSums(cellProp))*100

        mAge <- cbind(mAge, cellProp)
        arrow::write_feather(mAge, paste0(pipeline_data_dir, "mAge.feather"))
        arrow::write_feather(cellProp, paste0(pipeline_data_dir, "bloodcells.feather"))
    }
}

###################################
####    Preprocessing Results   ###
###################################

bloodcells <- c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")
clocks <- c("Horvath1", "Horvath2", "GrimAge1", "GrimAge2", "PCHorvath1", "PCHorvath2", "PCGrimAge")
predictors <- c(clocks, bloodcells)

# Estimate the technical measurement noise using the replicate samples for each predictor and pipeline
Noise <- data.frame(row.names = predictors)
for (method in 1:102) {

    if (method %in% c(86, 87, 88, 89, 92, 96, 97, 99, 100)) next

    pipeline <- method_to_pipeline(method)
    cat("Calculating measurement noise:", method, pipeline, "\n")

    pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")

    mAge <- arrow::read_feather(paste0(pipeline_data_dir, "mAge.feather"))
    mAge <- as.data.frame(mAge)
    rownames(mAge) <- mAge$id

    load("Data/Objects/pheno.Rdata")

    mAge <- cbind(pheno, mAge[rownames(pheno), predictors])

    noise_reps <- mAge[mAge$has_replicate == TRUE, c("subject_id", "exercise_timepoint", predictors)] %>%
        group_by(subject_id, exercise_timepoint) %>%
        summarize_at(predictors, function(x) abs(diff(x))) %>%
        ungroup() %>%
        summarize_at(predictors, function(x) mean(x)) %>%
        as.matrix() %>%
        as.vector()
    
    Noise[,pipeline] <- noise_reps
}

save(Noise, file = "Data/Objects/Noise.Rdata")

# Retrieve pipeline with the lowest measurement noise for each predictor
OptimalPipelines <- data.frame(
    Predictor = c(clocks, bloodcells),
    Pipeline = rep("", length(clocks) + length(bloodcells)),
    Noise = rep(NA, length(clocks) + length(bloodcells))
)
for (clock in clocks) {
    ix <- which(OptimalPipelines$Predictor == clock)
    OptimalPipelines[ix, "Pipeline"] <- names(which.min(Noise[clock,]))
    OptimalPipelines[ix, "Noise"] <- min(Noise[clock,])
}
# Take the pipeline that has the lowest mean z-scaled (per cell type) measurement noise, over all cell types
Tech_bloodcells <- t(scale(t(Noise[bloodcells,]), center = TRUE, scale = TRUE))
Tech_bloodcells <- colMeans(Tech_bloodcells)
pipeline_bloodcells <- names(which.min(Tech_bloodcells))
for (bc in bloodcells) {
    ix <- which(OptimalPipelines$Predictor == bc)
    OptimalPipelines[ix, "Pipeline"] <- pipeline_bloodcells
    OptimalPipelines[ix, "Noise"] <- Noise[bc, pipeline_bloodcells]
}
save(OptimalPipelines, file = "Data/Objects/OptimalPipelines.Rdata")

# Horvath and GrimAge clocks with the lowest mean measurement error:
# PCHorvath1 clock with the enmix_est_nodye_nonorm_rcp has a mean measurement error of 0.3735464 years
# PCGrimAge clock with the enmix_oob_nodye_q2_norcp has a mean measurement error of 0.1583529 years
# watermelon_dasen is optimal for the leukocyte deconvolution

# Retrieve the optimal estimates for each predictor
load("Data/Objects/pheno.Rdata")

pipeline <- "enmix_est_nodye_nonorm_rcp"
pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")
mAge <- arrow::read_feather(paste0(pipeline_data_dir, "mAge.feather"))

pheno$epigenetic_age_Horvath <- mAge$PCHorvath1

pipeline <- "enmix_oob_nodye_q2_norcp"
pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")
mAge <- arrow::read_feather(paste0(pipeline_data_dir, "mAge.feather"))

pheno$epigenetic_age_GrimAge <- mAge$PCGrimAge

pipeline <- "watermelon_dasen"
pipeline_data_dir <- paste0("Results/PrePro/", pipeline, "/")
mAge <- arrow::read_feather(paste0(pipeline_data_dir, "mAge.feather"))

pheno$estimated_B_cells <- mAge$Bcell
pheno$estimated_CD4T_cells <- mAge$CD4T
pheno$estimated_CD8T_cells <- mAge$CD8T
pheno$estimated_natural_killer_cells <- mAge$NK
pheno$estimated_monocytes <- mAge$Mono
pheno$estimated_neutrophils <- mAge$Neu

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

# mean_noise_Horvath SEM_noise_Horvath mean_noise_GrimAge SEM_noise_GrimAge
#                4.48              1.84               1.90             0.63
