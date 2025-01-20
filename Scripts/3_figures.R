#*######################################################
#*######################################################
#*######################################################
#*###############    BIOCLOCK: Figures   ###############
#*######################################################
#*######################################################
#*######################################################

###################################
####           Set-up           ###
###################################

library(ggplot2)
library(ggnewscale)
library(patchwork)
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(stringr)
library(lme4)
library(lmerTest)

setwd(paste0("/", file.path("data", "user_homes", "mennovd", "BIOKLOK")))

dir.create("Images/Paper", recursive = TRUE, showWarnings = FALSE)

###################################
########   Data Loading    ########
###################################

load("Data/Objects/pheno.Rdata")
load("Data/Objects/bio.Rdata")
load("Data/Objects/d_bio.Rdata")

###################################
######  Plotting Functions   ######
###################################

source("Scripts/Additional/plottingFunctions.R")

##############################
########    Figures    #######
##############################

######## EA ifo CA

# make sure the "pre" timepoint comes first
pheno_plot <- pheno %>%
    mutate(exercise_timepoint = recode(exercise_timepoint, "pre" = "Pre", "post" = "Post")) %>%
    mutate(exercise_timepoint = factor(exercise_timepoint, levels = c("Pre", "Post")))
pos <- min(pheno[, "age"])

## GrimAge

clock <- "epigenetic_age_GrimAge"

# label the samples with replicates with a star
label_positions <- pheno_plot %>%
    dplyr::filter(
        has_replicate,
        exercise_timepoint == "Post"
    ) %>% # only for individuals with replicates
    group_by(subject_id) %>%
    summarize(
        age = mean(age), 
        clock_value = mean(!!sym(clock)),
        label = "*"
    )

plot_grimage <- ggplot(pheno_plot, aes(x = age, y = !!sym(clock))) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgray", linetype = "dashed", size = 0.35, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, color = "darkgray", linetype = "dotted", size = 0.35, show.legend = FALSE) +
    geom_path(aes(group = subject_id), size = 0.2, color = "black", show.legend = FALSE) +
    geom_text(
        data = label_positions,
        aes(x = age, y = clock_value, label = label),
        nudge_y = 0,
        nudge_x = 0.3,
        size = 2,
        vjust = 1,
        color = "black",
        show.legend = FALSE
    ) +
    geom_point(aes(color = exercise_timepoint), size = 0.6, , show.legend = TRUE) +
    scale_color_manual(
        name = "Timepoint",
        values = c("Pre" = colors[2], "Post" = colors[7])
    ) +
    mytheme +
    annotate( # trend line equation
        "text",
        x = max(pheno[, "age"]), y = min(pheno[, "age"]),
        label = paste("EA == ", format(round(coef(lm(pheno[, clock] ~ pheno[, "age"]))[2], 2), nsmall = 2), "%*% CA +", format(round(coef(lm(pheno[, clock] ~ pheno[, "age"]))[1], 2), nsmall = 2)),
        color = "black", hjust = 1, size = 1.5, parse = TRUE, show.legend = FALSE
    ) +
    xlab("CA (Years)") +
    labs(color = "Timepoint") +
    theme(
        axis.text.x = element_blank(), # remove tick text
        axis.text.y = element_blank(), # remove tick text
        axis.ticks = element_blank(), # remove ticks
    ) +
    annotate(
        "text",
        x = min(pheno[, "age"]), y = max(pheno[, clock]),
        label = "R^2 == 0.86 * ',' ~ P == 1.217 %*% 10^{-21}", # mixed model anova & partial R2
        hjust = 0, vjust = 1, color = "black", size = 1.5, parse = TRUE, show.legend = FALSE
    ) +
    labs(
        x = bquote(.(expression(CA ~ (Years)))),
        y = bquote(.(expression(GrimAge ~ EA ~ (Years))))
    ) +
    geom_point(data = data.frame(), aes(x = pos, y = pos), color = "darkgray", size = 0.6, show.legend = FALSE) +
    geom_point(data = data.frame(), aes(x = pos + 0.5, y = pos + 0.5), color = "darkgray", size = 0.6, show.legend = FALSE) +
    geom_line(data = data.frame(), aes(x = c(pos, pos + 0.5), y = c(pos, pos + 0.5)), color = "darkgray", size = 0.3, show.legend = FALSE)

## Horvath

clock <- "epigenetic_age_Horvath"

label_positions <- pheno_plot %>%
    dplyr::filter(
        has_replicate,
        exercise_timepoint == "Post"
    ) %>% # only for individuals with replicates
    group_by(subject_id) %>%
    summarize(
        age = mean(age), 
        clock_value = mean(!!sym(clock)),
        label = "*"
    )

plot_horvath <- ggplot(pheno_plot, aes(x = age, y = !!sym(clock))) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgray", linetype = "dashed", size = 0.35, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, color = "darkgray", linetype = "dotted", size = 0.35, show.legend = FALSE) +
    geom_path(aes(group = subject_id), size = 0.2, color = "black", show.legend = FALSE) +
    geom_text(
        data = label_positions,
        aes(x = age, y = clock_value, label = label),
        nudge_y = 0,
        nudge_x = 0.3,
        size = 2,
        vjust = 1,
        color = "black",
        show.legend = FALSE
    ) +
    geom_point(aes(color = exercise_timepoint), size = 0.6, show.legend = TRUE) +
    scale_color_manual(
        name = "Timepoint",
        values = c("Pre" = colors[2], "Post" = colors[7])
    ) +
    mytheme +
    annotate( # trend line equation
        "text",
        x = max(pheno[, "age"]), y = min(pheno[, "age"]),
        label = paste("EA == ", format(round(coef(lm(pheno[, clock] ~ pheno[, "age"]))[2], 2), nsmall = 2), "%*% CA +", format(round(coef(lm(pheno[, clock] ~ pheno[, "age"]))[1], 2), nsmall = 2)),
        color = "black", hjust = 1, size = 1.5, parse = TRUE,
        show.legend = FALSE
    ) +
    xlab("CA (Years)") +
    labs(color = "Timepoint") +
    theme(
        axis.text.x = element_blank(), # remove tick text
        axis.text.y = element_blank(), # remove tick text
        axis.ticks = element_blank(), # remove ticks
    ) +
    annotate(
        "text",
        x = min(pheno[, "age"]), y = max(pheno[, clock]),
        label = "R^2 == 0.71 * ',' ~ P == 1.764 %*% 10^{-13}",
        hjust = 0, vjust = 1, color = "black", size = 1.5, parse = TRUE,
        show.legend = FALSE
    ) +
    labs(
        x = bquote(.(expression(CA ~ (Years)))),
        y = bquote(.(expression(Horvath ~ EA ~ (Years))))
    ) +
    geom_point(data = data.frame(), aes(x = pos, y = pos), color = "darkgray", size = 0.6, show.legend = FALSE) +
    geom_point(data = data.frame(), aes(x = pos + 0.5, y = pos + 0.5), color = "darkgray", size = 0.6, show.legend = FALSE) +
    geom_line(data = data.frame(), aes(x = c(pos, pos + 0.5), y = c(pos, pos + 0.5)), color = "darkgray", size = 0.3, show.legend = FALSE)

dplot_horvath_acc <- deltaPlot(
    data = bio,
    vars = c("epigenetic_age_acceleration_Horvath"),
    ncol = 1, 
    dropouts = bio %>% filter(dropout | excluded) %>% pull(subject_id), # remove dropouts
    plotsegments = "all",
    ytitles = expression(Horvath ~ EAA ~ (Years)),
    xtitles = expression(Timepoint),
    display = c("effect", "pval", "mean"),
    collectaxis = FALSE,
    collectguides = TRUE,
    legend = FALSE,
    dim = 2000,
    width = 1.1
)

dplot_grimage_acc <- deltaPlot(
    data = bio,
    vars = c("epigenetic_age_acceleration_GrimAge"),
    ncol = 1, 
    dropouts = bio %>% filter(dropout | excluded) %>% pull(subject_id), # remove dropouts
    plotsegments = "all",
    ytitles = expression(GrimAge ~ EAA ~ (Years)),
    xtitles = expression(Timepoint),
    display = c("effect", "pval", "mean"),
    collectaxis = FALSE,
    collectguides = TRUE,
    legend = FALSE,
    dim = 2000,
    width = 1.1
)

combined_plot_horvath <- (plot_horvath | dplot_horvath_acc) +
    plot_layout(widths = c(2, 1), guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 9))
ggsave("Images/Paper/Horvath_EA.jpg", combined_plot_horvath, dpi = 1200, width = 16, height = 6, units = "cm", bg = "white", device = "jpg")

combined_plot_grimage <- (plot_grimage | dplot_grimage_acc) +
    plot_layout(widths = c(2, 1), guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 9))
ggsave("Images/Paper/GrimAge_EA.jpg", combined_plot_grimage, dpi = 1200, width = 16, height = 6, units = "cm", bg = "white", device = "jpg")


######### GrimAge correlation with VO2max & Neu

VO2max_plot <- correlationPlot(
    data = d_bio, 
    yvars = "epigenetic_age_acceleration_GrimAge", 
    xvars = "vo2max_per_kg",
    ncol = 1, 
    xtitles = expression(Delta * VO[2] ~ max ~ (ml/kg/min)), 
    ytitles = expression(GrimAge ~ Delta * EAA ~ (Years)),
    labelxposition = "max",
    labelyposition = "max",
    color = colors[5]
)

Neu_plot <- correlationPlot(
    data = d_bio, 
    yvars = "epigenetic_age_acceleration_GrimAge", 
    xvars = "neutrophils", 
    ncol = 1, 
    xtitles = expression(Delta * Neutrophils ~ ("%")),
    ytitles = expression(GrimAge ~ Delta * EAA ~ (Years)),
    labelxposition = "min",
    labelyposition = "max",
    color = colors[3]
)

combined_plot_corr <- (VO2max_plot | Neu_plot) +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 9))
ggsave("Images/Paper/GrimAge_Correlations.jpg", combined_plot_corr, dpi = 1200, width = 16, height = 6, units = "cm", bg = "white", device = "jpg")


######### (D) EAA correlations

plot_data <- bio %>%
    mutate(exercise_timepoint = factor(recode(exercise_timepoint, "pre" = "Pre", "post" = "Post"), levels = c("Pre", "Post"))) # set order right

EAA_plot <- ggplot(plot_data, aes(x = epigenetic_age_acceleration_Horvath, y = epigenetic_age_acceleration_GrimAge, color = exercise_timepoint)) +
    mytheme +
    labs(
        x = expression(Horvath ~ EAA ~ (Years)),
        y = expression(GrimAge ~ EAA ~ (Years))
    ) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgray", linetype = "dashed", size = 0.35) +
    geom_path(aes(group = subject_id), color = "black", size = 0.25) +
    geom_point(size = 0.6) +
    scale_color_manual(values = c("Pre" = colors[2], "Post" = colors[7])) +
    labs(color = "Timepoint")

D_EAA_plot <- correlationPlot(
    data = d_bio, 
    yvars = "epigenetic_age_acceleration_GrimAge", 
    xvars = "epigenetic_age_acceleration_Horvath", 
    ncol = 1, 
    xtitles = expression(Horvath ~ Delta * EAA ~ (Years)),
    ytitles = expression(GrimAge ~ Delta * EAA ~ (Years)),
    labelxposition = "min",
    labelyposition = "max",
    color = colors[1]
)

combined_plot <- (EAA_plot | D_EAA_plot) +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 9))
ggsave("Images/Paper/Clocks_Association.jpg", combined_plot, dpi = 1200, width = 16, height = 6, units = "cm", bg = "white", device = "jpg")

########## Delta Epigenetic age acceleration & bloodcells correlation

bloodcells <- c(
    "B_cells_naive",
    "B_cells_memory",
    "CD4T_cells_naive",
    "CD4T_cells_memory",
    "CD8T_cells_naive",
    "CD8T_cells_memory",
    "T_regulatory_cells",
    "natural_killer_cells",
    "eosinophils",
    "basophils",
    "monocytes",
    "neutrophils"
)

bloodcell_names <- c(
    expression(B ~ Lymphocytes ~ Naive),
    expression(B ~ Lymphocytes ~ Memory),
    expression(CD4^{"+"} ~ T ~ Lymphocytes ~ Naive),
    expression(CD4^{"+"} ~ T ~ Lymphocytes ~ Memory),
    expression(CD8^{"+"} ~ T ~ Lymphocytes ~ Naive),
    expression(CD8^{"+"} ~ T ~ Lymphocytes ~ Memory),
    expression(Regulatory ~ T ~ Lymphocytes),
    expression(Natural ~ Killer ~ Lymphocytes),
    expression(Eosinophils),
    expression(Basophils),
    expression(Monocytes),
    expression(Neutrophils)
)

names(bloodcell_names) <- bloodcells

color_vector <- c(
  eosinophils = "#B15928",
  B_cells_memory = "#FFFF99",
  CD8T_cells_naive = "#6A3D9A",
  basophils = "#CAB2D6",
  T_regulatory_cells = "#FF7F00",
  B_cells_naive = "#FDBF6F",
  natural_killer_cells = "#E31A1C",
  CD8T_cells_memory = "#FB9A99",
  monocytes = "#33A02C",
  CD4T_cells_naive = "#B2DF8A",
  CD4T_cells_memory = "#1F78B4",
  neutrophils = "#A6CEE3"
)

d_titles <- as.vector(lapply(bloodcell_names, function(name) {
    bquote(Delta * .(name) ~ ("%"))
}))

BC_grimage_plot <- correlationPlot(
    data = d_bio,
    yvars = "epigenetic_age_acceleration_GrimAge",
    xvars = bloodcells, 
    ncol = 3, 
    xtitles = d_titles,
    ytitles = rep(expression(GrimAge ~ Delta * EAA ~ (Years)), length(bloodcell_names)),
    color = color_vector[bloodcells],
    sortplotsbyX = FALSE
)
combined_plot <- BC_grimage_plot +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 9)
)
ggsave("Images/Paper/GrimAge_Bloodcells_Associations.jpg", combined_plot, dpi = 1200, width = 16, height = 18, units = "cm", bg = "white", device = "jpg")

BC_horvath_plot <- correlationPlot(
    data = d_bio,
    yvars = "epigenetic_age_acceleration_Horvath",
    xvars = bloodcells, 
    ncol = 3, 
    xtitles = d_titles,
    ytitles = rep(expression(Horvath ~ Delta * EAA ~ (Years)), length(bloodcell_names)),
    color = color_vector[bloodcells],
    sortplotsbyX = FALSE
)
combined_plot <- BC_horvath_plot +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 9)
)
ggsave("Images/Paper/Horvath_Bloodcells_Associations.jpg", combined_plot, dpi = 1200, width = 16, height = 18, units = "cm", bg = "white", device = "jpg")

######### Mean leukocyte composition bar chart

leukocytes <- bio %>%
    filter(!dropout, !excluded) %>%
    dplyr::select(exercise_timepoint, all_of(bloodcells)) %>%
    filter(exercise_timepoint %in% c("pre", "post")) %>%
    mutate(exercise_timepoint = factor(recode(exercise_timepoint, "pre" = "Pre", "post" = "Post"), levels = c("Post", "Pre"))) %>%
    pivot_longer(cols = -exercise_timepoint, names_to = "leukocyte", values_to = "concentration") %>%
    group_by(exercise_timepoint, leukocyte) %>%
    summarize(across(everything(), mean, na.rm = TRUE)) %>%
    arrange(concentration) %>%
    mutate(
        leukocyte = factor(leukocyte, levels = unique(leukocyte)),
        concentration = round(concentration, 2)
    )

bar_chart <- ggplot(leukocytes, aes(x = exercise_timepoint, y = concentration, fill = leukocyte)) +
    geom_bar(stat = "identity", position = "fill") +  # position = "fill" for proportional bars
    coord_flip() +  # Flip for horizontal bars
    scale_fill_manual(
        values = color_vector,
        labels = bloodcell_names
    ) +
    mytheme +
    theme(
                legend.title = element_blank()
    ) +
    geom_text(
        aes(label = ifelse(concentration > 0.1, concentration, "")),
        position = position_fill(vjust = 0.5),
        size = 0.8, angle = 90
    ) +
    geom_text(
        aes(label = ifelse(concentration < 0.1, concentration, "")),
        position = position_fill(vjust = 0.5),
        size = 0.8, angle = 90, vjust = 1,
    ) +
    scale_y_continuous(
        labels = function(x) x * 100,  # Horizontal axis in percentages
        breaks = seq(0, 1, 0.2)  # Define the tick marks
    ) +
    labs(
        title = NULL,
        x = NULL,
        y = "Proportion (%)"
    ) +
    guides(fill = guide_legend(reverse = TRUE))  # Reverse the order of the legend

# Deltaplots

BCs <- deltaPlot(
    data = bio, 
    vars = bloodcells, 
    ncol = 3, 
    dropouts = bio %>% filter(dropout | excluded) %>% pull(subject_id),
    plotsegments = "all",
    display =  c("mean", "effect", "pval"),
    ytitles = bloodcell_names,
    legend = TRUE,
    collectaxis = FALSE,
    width = 1.1
)

combined_leukocytes <- bar_chart / BCs + 
    plot_layout(heights = c(1, 4)) +
    plot_annotation(tag_levels = "a")  & 
    theme(plot.tag = element_text(size = 9))
ggsave("Images/Paper/Leukocyte_Composition.jpg", combined_leukocytes, dpi = 1200, width = 16, height = 27, units = "cm", bg = "white", device = "jpg")


######### Training time

trainingtime <- read_xlsx("Data/Pheno/Sample Data.xlsx", sheet = "TrainingTime")

trainingtime <- trainingtime %>%
    mutate(
        week = week - min(week),
        adherence = (executed / planned) * 100,
        # if planned is 0, set adherence to NA
        adherence = if_else(planned == 0, 100, adherence)
    ) %>%
    filter(
        !(subject_id %in% (bio %>% filter(dropout | excluded) %>% pull(subject_id)))
    )

# calculate mean & SE per week
trainingtime_means <- trainingtime %>%
    gather(variable, value, planned:adherence) %>%
    group_by(week, variable) %>%
    summarize(
        mean_value = mean(value, na.rm = TRUE)
    )
trainingtime_errors <- trainingtime %>%
    gather(variable, value, planned:adherence) %>%
    group_by(week, variable) %>%
    summarize(
        mean_value = mean(value, na.rm = TRUE),
        error = sd(value, na.rm = TRUE) / sqrt(n())
    )

plot <- ggplot() +
    geom_errorbar(
        data = filter(trainingtime_errors, variable %in% c("executed", "planned")),
        aes(x = week, ymin = mean_value - error, ymax = mean_value + error, color = variable),
        width = 0.2,
        size = 0.25
    ) +
    # adherence is scaled by 12 to match the hours scale
    geom_errorbar( 
        data = filter(trainingtime_errors, variable == "adherence"),
        aes(x = week, ymin = (mean_value - error)/12, ymax = (mean_value + error)/12, color = variable),
        width = 0.2,
        size = 0.25
    ) +
    geom_point(
        data = filter(trainingtime_means, variable %in% c("executed", "planned")),
        aes(x = week, y = mean_value, color = variable), 
        size = 0.5
    ) +
    geom_line(
        data = filter(trainingtime_means, variable %in% c("executed", "planned")),
        aes(x = week, y = mean_value, color = variable),
        size = 0.25
    ) +
    geom_point(
        data = filter(trainingtime_means, variable == "adherence"),
        aes(x = week, y = mean_value / 12, color = variable), 
        size = 0.25
    ) +
    geom_line(
        data = filter(trainingtime_means, variable == "adherence"),
        aes(x = week, y = mean_value / 12, color = variable),
        size = 0.25
    ) +
    geom_hline(yintercept = 100 / 12, linetype = "dotted", color = colors[7], size = 0.25) +
    scale_y_continuous(breaks = seq(2, 20, 2), sec.axis = sec_axis(~ . * 12, name = bquote(.(expression(Adherence ~ ("%")))), breaks = seq(0, 200, 25))) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
    scale_color_manual(
        values = c("executed" = colors[5], "planned" = colors[3], "adherence" = colors[7]),
        breaks = c("planned", "executed", "adherence"),
        labels = c("Executed", "Planned", "Adherence")   # Capitalized labels
    ) +
    labs(
        y = bquote(.(expression(Planned ~ "&" ~ Executed ~ Training ~ (Hours)))),
        x = bquote(.(expression(Week)))
    ) +
    mytheme +
    theme(
        legend.position = "right",
        legend.title = element_blank()
    ) +
    guides(color = guide_legend(title = NULL))
ggsave("Images/Paper/Training_Time.jpg", plot, dpi = 1200, width = 16, height = 6, units = "cm", bg = "white", device = "jpg")
