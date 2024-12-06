#*######################################################
#*######################################################
#*######################################################
#*#########    BIOCLOCK: Statistical Analysis   ########
#*######################################################
#*######################################################
#*######################################################

###################################
####           Set-up           ###
###################################

library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(stringr)
library(feather)
library(arrow)
library(lme4)
library(lmerTest)
library(MuMIn)

setwd(paste0("/", file.path("data", "user_homes", "mennovd", "BIOKLOK")))

dir.create("Results/Analysis", recursive = TRUE, showWarnings = FALSE)

###################################
####        Data Loading        ###
###################################

load("Data/Objects/pheno.Rdata")

clocks <- c(
    "epigenetic_age_GrimAge",
    "epigenetic_age_Horvath"
)

bloodcells <- c(
    "estimated_B_cells",
    "estimated_CD4T_cells",
    "estimated_CD8T_cells",
    "estimated_natural_killer_cells",
    "estimated_monocytes",
    "estimated_neutrophils"
)

predictors <- c(clocks, bloodcells)

biovars <- c(
    "vo2max_per_kg",
    "peak_power_output_per_kg",
    "hand_grip_strength",
    "body_mass_index",
    "fat_percent",
    "lean_mass",
    "bone_density",
    "nightly_heart_rate",
    "nightly_heart_rate_variability",
    "systolic_blood_pressure",  
    "diastolic_blood_pressure",
    "pulse_wave_velocity",
    "sleep_score",
    "sleep_duration",
    "diet_quality_score"
)

###################################
####        Data Summary        ###
###################################

# Calculate mean epigenetic ages & blood cell fractions for replicates
bio <- pheno[, c(
    "subject_id", 
    "exercise_timepoint", 
    predictors
    )] %>%
    group_by(subject_id, exercise_timepoint) %>%
    summarize_at(predictors, mean) %>%
    arrange(subject_id, exercise_timepoint) %>% 
    ungroup()

# Select relevant columns
pheno_subset <- pheno[, c(
        "subject_id",
        "exercise_timepoint", 
        "sex",
        "age",
        "vo2max_per_kg", 
        "peak_power_output_per_kg", 
        "hand_grip_strength",
        "weight",
        "height",
        "body_mass_index",
        "fat_percent",
        "lean_mass",
        "bone_density",
        "nightly_heart_rate", 
        "nightly_heart_rate_variability",
        "systolic_blood_pressure",
        "diastolic_blood_pressure",
        "pulse_wave_velocity",
        "sleep_duration",
        "sleep_score",
        "diet_quality_score",
        "planned_training",
        "executed_training",
        "training_adherence",
        "dropout",
        "excluded"
        )
    ] %>% unique()

# Merge pheno with bio
bio <- merge(bio, pheno_subset, by = c("subject_id", "exercise_timepoint"))

# Make blood cell fractions sum to 1 for replicates
bio[,bloodcells] <- (bio[,bloodcells] / rowSums(bio[,bloodcells]))*100

# Calculate epigenetic age acceleration
GrimAge_model <- lm(epigenetic_age_GrimAge ~ age, bio)
bio$epigenetic_age_acceleration_GrimAge <-  residuals(GrimAge_model)

Horvath_model <- lm(epigenetic_age_Horvath ~ age, bio)
bio$epigenetic_age_acceleration_Horvath <- residuals(Horvath_model)

save(bio, file = "Data/Objects/bio.Rdata")

# Longitudinal changes in variables
d_bio <- bio %>%
    arrange(subject_id, exercise_timepoint) %>%
    group_by(subject_id) %>%
    reframe(across(where(is.numeric), ~ -diff(.))) %>%
    ungroup()

save(d_bio, file = "Data/Objects/d_bio.Rdata")

###################################
#####  Changes in variables   #####
###################################

table_data <- bio %>% 
    filter(!excluded, !dropout) %>% 
    select(all_of(c("subject_id", "exercise_timepoint", bloodcells, biovars)))

results_table <- data.frame(variable = character(), mean_pre = numeric(), sd_pre = numeric(), 
                            mean_post = numeric(), sd_post = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

for (variable in colnames(table_data)) {
    if(is.numeric(table_data[[variable]])) {
    m_pre <- mean(table_data %>% filter(exercise_timepoint == "pre") %>% pull(variable), na.rm = TRUE)
    sd_pre <- sd(table_data %>% filter(exercise_timepoint == "pre") %>% pull(variable), na.rm = TRUE)
    m_post <- mean(table_data %>% filter(exercise_timepoint == "post") %>% pull(variable), na.rm = TRUE)
    sd_post <- sd(table_data %>% filter(exercise_timepoint == "post") %>% pull(variable), na.rm = TRUE)
    formula <- paste(variable, "~ (1|subject_id) + factor(exercise_timepoint)")
    model <- lmer(formula, data = table_data)
    p <- summary(model)$coefficients["factor(exercise_timepoint)pre", "Pr(>|t|)"]
    new_row <- data.frame(variable = variable, mean_pre = m_pre, sd_pre = sd_pre, 
                            mean_post = m_post, sd_post = sd_post, p_value = p, stringsAsFactors = FALSE)
    results_table <- rbind(results_table, new_row)
    }
}
write_xlsx(results_table, "Results/Analysis/Table_Variable_Changes.xlsx")

###################################
#######   clock statistics  #######
###################################

### Difference in EAA in males vs females, pre-EET

## GrimAge clock
t <- t.test(
    bio %>% filter(sex == "male", exercise_timepoint == "pre") %>% pull(epigenetic_age_acceleration_GrimAge), 
    bio %>% filter(sex == "female", exercise_timepoint == "pre") %>% pull(epigenetic_age_acceleration_GrimAge),
    var.equal = FALSE
)
print(t$p.value)
# 0.001499333
print(t$estimate[1] - t$estimate[2])
# 2.80849

## Horvath clock
t <- t.test(
    bio %>% filter(sex == "male", exercise_timepoint == "pre") %>% pull(epigenetic_age_acceleration_Horvath), 
    bio %>% filter(sex == "female", exercise_timepoint == "pre") %>% pull(epigenetic_age_acceleration_Horvath),
    var.equal = FALSE
)
print(t$p.value)
# 0.09842557
print(t$estimate[1] - t$estimate[2])
# 2.011997

### Calculate partial marginal variance of EA explained by CA (including replicates)

##  GrimAge clock

# Build the full and reduced mixed models
model <- lmer(epigenetic_age_GrimAge ~ (1|subject_id) + factor(exercise_timepoint) + factor(sex) + age, data = pheno)
reduced_model <- lmer(epigenetic_age_GrimAge ~ (1|subject_id) + factor(exercise_timepoint) + factor(sex), data = pheno)

# Perform a likelihood ratio test to compare the models (p-val for variance explained by age)
anova_result <- anova(model, reduced_model)
p_value <- anova_result$`Pr(>Chisq)`[2]
print(p_value)
# 2.099036e-20

# Calculate partial marginal R² for age (MuMin package)
r2_full <- r.squaredGLMM(model)
print(r2_full)
r2_reduced <- r.squaredGLMM(reduced_model)
print(r2_reduced)
partial_r2_Age <- r2_full[1] - r2_reduced[1]
print(partial_r2_Age)
# 0.8383908

## Horvath clock

# Build the full and reduced mixed models
model <- lmer(epigenetic_age_Horvath ~ (1|subject_id) + factor(exercise_timepoint) + factor(sex) + age, data = pheno)
reduced_model <- lmer(epigenetic_age_Horvath ~ (1|subject_id) + factor(exercise_timepoint) + factor(sex), data = pheno)

# Perform a likelihood ratio test to compare the models (p-val for variance explained by age)
anova_result <- anova(model, reduced_model)
p_value <- anova_result$`Pr(>Chisq)`[2]
print(p_value)
# 1.333332e-13

# Calculate partial marginal R² for age (MuMin package)
r2_full <- r.squaredGLMM(model)
print(r2_full)
r2_reduced <- r.squaredGLMM(reduced_model)
print(r2_reduced)
partial_r2_Age <- r2_full[1] - r2_reduced[1]
print(partial_r2_Age)
# 0.710164

### Difference in EAA pre- and post-EET (in months)

## unadjusted for blood cells

# GrimAge clock
model <- lmer(epigenetic_age_acceleration_GrimAge*12 ~ (1|subject_id) + factor(exercise_timepoint), data = bio %>% filter(!dropout, !excluded))
print(summary(model)$coefficients)
#                                Estimate Std. Error       df    t value   Pr(>|t|)
# (Intercept)                   -3.532929   6.313744 35.90903 -0.5595617 0.57925229
# factor(exercise_timepoint)pre  7.596151   3.034977 32.00000  2.5028691 0.01761841

# Horvath clock
model <- lmer(epigenetic_age_acceleration_Horvath*12 ~ (1|subject_id) + factor(exercise_timepoint), data = bio %>% filter(!dropout, !excluded))
print(summary(model)$coefficients)
#                                Estimate Std. Error       df    t value   Pr(>|t|)
# (Intercept)                   -7.444666   8.211871 33.32073 -0.9065737 0.37114185
# factor(exercise_timepoint)pre  6.099630   2.335849 32.00000  2.6113115 0.01361605

## Adjusted for blood cells

# GrimAge clock
model <- lmer(epigenetic_age_acceleration_GrimAge*12 ~ (1|subject_id) + factor(exercise_timepoint) + estimated_neutrophils, data = bio %>% filter(!dropout, !excluded))
print(summary(model)$coefficients)
#                                  Estimate Std. Error       df   t value     Pr(>|t|)
# (Intercept)                   -113.117896  13.226149 48.62557 -8.552595 2.939615e-11
# factor(exercise_timepoint)pre    7.131758   1.646289 31.13104  4.332021 1.430306e-04
# estimated_neutrophils          201.344035  21.957010 35.57866  9.169920 6.668941e-11

# Horvath clock
model <- lmer(epigenetic_age_acceleration_Horvath*12 ~ (1|subject_id) + factor(exercise_timepoint) + estimated_CD4T_cells, data = bio %>% filter(!dropout, !excluded))
print(summary(model)$coefficients)
#                                  Estimate Std. Error       df   t value     Pr(>|t|)
# (Intercept)                     34.583881  10.424628 57.98665  3.317517 1.572209e-03
# factor(exercise_timepoint)pre    5.653304   1.510155 30.95979  3.743526 7.421053e-04
# estimated_CD4T_cells          -314.676683  47.128469 35.16725 -6.676998 9.788410e-08

### Correlation of d_EAA between both clocks
model <- lm(epigenetic_age_acceleration_Horvath ~ epigenetic_age_acceleration_GrimAge, d_bio)
model <- lm(epigenetic_age_acceleration_GrimAge ~ epigenetic_age_acceleration_Horvath, d_bio)
print(summary(model))
# Residual standard error: 1.195 on 36 degrees of freedom
# Multiple R-squared:  0.2978,    Adjusted R-squared:  0.2783 
# F-statistic: 15.27 on 1 and 36 DF,  p-value: 0.0003949

###########################################
####  Associations d_EAA ~ d_variable  ####
###########################################

### Associations with blood cells

AssociationsTest <- function(data, outcome, predictors) {
    # Scale numeric columns
    numeric_columns <- sapply(data, is.numeric)
    data_numeric <- scale(data[, numeric_columns])

    # Helper function to test a single predictor
    ModelTest <- function(data, outcome, predictor) {
        formula <- as.formula(paste0(outcome, " ~ ", predictor))
        model <- lm(formula, data = as.data.frame(data), na.action = na.omit)
        summary_model <- summary(model)
        
        # Extract and return results
        list(
            Predictor = predictor,
            Coefficient = summary_model$coefficients[2, 1],
            P_value = summary_model$coefficients[2, 4],
            R_squared = summary_model$r.squared
        )
    }

    # Iterate over predictors and apply ModelTest
    results <- lapply(predictors, function(predictor) {
        ModelTest(data_numeric, outcome, predictor)
    })

    # Combine results into a data frame
    results_df <- do.call(rbind, lapply(results, as.data.frame))
    rownames(results_df) <- NULL
    return(results_df)
}

GrimAge_BC_associations <- AssociationsTest(
    data = d_bio,
    outcome = "epigenetic_age_acceleration_GrimAge", 
    predictors = bloodcells
)
write_xlsx(GrimAge_BC_associations, "Results/Analysis/BloodCell_Associations_GrimAge.xlsx")

Horvath_BC_associations <- AssociationsTest(
    data = d_bio,
    outcome = "epigenetic_age_acceleration_Horvath", 
    predictors = bloodcells
)
write_xlsx(Horvath_BC_associations, "Results/Analysis/BloodCell_Associations_Horvath.xlsx")

### Associations with other variables, unadjusted and adjusted for blood cells

adjustedAssociationsTest <- function(data, outcome, predictors, adjust) {
    # Scale numeric columns
    numeric_columns <- sapply(data, is.numeric)
    data_numeric <- scale(data[, numeric_columns])

    # Helper function to test a single predictor
    adjustedModelTest <- function(data, outcome, predictor, adjust) {
        # Unadjusted model
        unadjusted_formula <- as.formula(paste0(outcome, " ~ ", predictor))
        unadjusted_model <- lm(unadjusted_formula, data = as.data.frame(data), na.action = na.omit)
        unadjusted_summary <- summary(unadjusted_model)

        # Adjusted model
        adjusted_formula <- as.formula(paste0(outcome, " ~ ", predictor, " + ", adjust))
        adjusted_model <- lm(adjusted_formula, data = as.data.frame(data), na.action = na.omit)
        adjusted_summary <- summary(adjusted_model)

        # Return results
        list(
            Predictor = predictor,
            Coefficient_Unadjusted = unadjusted_summary$coefficients[2, 1],
            P_value_Unadjusted = unadjusted_summary$coefficients[2, 4],
            R_squared_Unadjusted = unadjusted_summary$r.squared,
            Coefficient_Adjusted = adjusted_summary$coefficients[2, 1],
            P_value_Adjusted = adjusted_summary$coefficients[2, 4],
            R_squared_Adjusted = adjusted_summary$r.squared
        )
    }

    # Iterate over predictors and apply adjustedModelTest
    results <- lapply(predictors, function(predictor) {
        adjustedModelTest(data_numeric, outcome, predictor, adjust)
    })

    # Combine results into a data frame
    results_df <- do.call(rbind, lapply(results, as.data.frame))
    rownames(results_df) <- NULL
    return(results_df)
}

GrimAge_associations <- adjustedAssociationsTest(
    data = d_bio,
    outcome = "epigenetic_age_acceleration_GrimAge", 
    predictors = c(biovars, "executed_training", "training_adherence"),
    adjust = "estimated_neutrophils"
)
write_xlsx(GrimAge_associations, "Results/Analysis/Associations_GrimAge.xlsx")

Horvath_associations <- adjustedAssociationsTest(
    data = d_bio,
    outcome = "epigenetic_age_acceleration_Horvath", 
    predictors = c(biovars, "executed_training", "training_adherence"),
    adjust = "estimated_CD4T_cells"
)
write_xlsx(Horvath_associations, "Results/Analysis/Associations_Horvath.xlsx")
