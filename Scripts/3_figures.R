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

pvalue_to_stars <- function(pvalue) {
    if (is.na(pvalue)) {
        return("")
    } else if (pvalue < 0.001) {
        return("***")
    } else if (pvalue < 0.01) {
        return("**")
    } else if (pvalue < 0.05) {
        return("*")
    } else {
        return("")
    }
}

mytheme <- theme_minimal() +
    theme(
        axis.title.x = element_text(margin = margin(t = 1, b = 0), size = 6),
        axis.title.y = element_text(margin = margin(r = 1, l = 0), size = 6),
        axis.title.y.right = element_text(margin = margin(l = 1, r = 0), size = 6),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_text(margin = margin(b = 7, t = 7), hjust = 0.5, size = 7),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.25, 'cm'),  # Adjust vertical spacing between legend items
        legend.spacing.x = unit(0.50, 'cm'),  # Adjust horizontal spacing between legend items
        legend.key.height = unit(0.25, 'cm'),  # Adjust the height of legend keys
        legend.key.width = unit(0.25, 'cm'),    # Adjust the width of legend keys
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.3)
    )

color_scheme <- palette.colors(palette = "Okabe-Ito")[c(
        "reddishpurple", # 1
        "vermillion", # 2
        "orange", # 3
        "yellow", # 4
        "bluishgreen", # 5
        "skyblue", # 6
        "blue", # 7
        "black", # 8
        "gray" # 9
    )
]
colors <- as.vector(color_scheme)

deltaPlot <- function(
    data, 
    vars, 
    dropouts = data %>% filter(dropout | excluded) %>% pull(subject_id),
    file = NULL,
    titles = c(),
    ytitles = c(),
    xtitles = c(),
    addpoints = FALSE,
    removeyaxisticks = FALSE,
    collectguides = TRUE,
    collectaxis = FALSE,
    legend = FALSE,
    plotsegments = c("all", "sign", "no"),
    display = c("none", "pval", "pstar", "effect", "mean"),
    ncol = 1,
    nrow = ceiling(length(vars) / ncol),
    dim = 2000,
    width = 1,
    colors = as.vector(palette.colors(palette = "Okabe-Ito")[c("reddishpurple", "vermillion", "orange", "yellow", "bluishgreen", "skyblue", "blue", "black", "gray")])
) {
    #' deltaPlot: Create Violin and Boxplots for Paired Data with Optional Significance Segments
    #'
    #' @description
    #' The `deltaPlot` function generates paired violin and boxplots for longitudinal or paired data, 
    #' highlighting differences between two timepoints (e.g., Pre vs Post). It includes options to 
    #' add significance segments with customizable annotations, and can combine multiple plots into 
    #' a grid layout.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param vars A character vector of variable names to plot. Each variable is plotted separately.
    #' @param dropouts A character vector of `subject_id` values to exclude from the analysis. 
    #'   Defaults to subjects who are marked as `dropout` or `excluded`.
    #' @param file A string specifying the output file path for saving the combined plot. Defaults to `NULL` (no saving).
    #' @param titles A character vector of titles for the plots. Defaults to an empty vector.
    #' @param ytitles A character vector of y-axis titles for the plots. Defaults to an empty vector.
    #' @param xtitles A character vector of x-axis titles for the plots. Defaults to an empty vector.
    #' @param addpoints Logical. If `TRUE`, adds individual data points to the plots. Defaults to `FALSE`.
    #' @param removeyaxisticks Logical. If `TRUE`, removes y-axis tick marks from the plots. Defaults to `FALSE`.
    #' @param collectguides Logical. If `TRUE`, collects the guides (e.g., legends) for all subplots. Defaults to `TRUE`.
    #' @param collectaxis Logical. If `TRUE`, aligns the axes across subplots. Defaults to `FALSE`.
    #' @param legend Logical. If `TRUE`, includes the legend in each subplot. Defaults to `FALSE`.
    #' @param plotsegments Character. Determines whether to plot significance segments. Options are:
    #'   - `"all"`: Always plot segments.
    #'   - `"sign"`: Only plot segments if the p-value is significant.
    #'   - `"no"`: Do not plot any segments.
    #'   Defaults to `"all"`.
    #' @param display Character vector. Specifies what information to display above the segments. Options are:
    #'   - `"none"`: No text above the segments.
    #'   - `"pval"`: Display p-values.
    #'   - `"pstar"`: Display significance stars.
    #'   - `"effect"`: Display effect sizes (e.g., mean differences).
    #'   - `"mean"`: Display mean values of each timepoint.
    #'   Defaults to `c("none", "pval", "pstar", "effect", "mean")`.
    #' @param ncol Integer. Number of columns in the combined plot layout. Defaults to `1`.
    #' @param nrow Integer. Number of rows in the combined plot layout. Defaults to `ceiling(length(vars) / ncol)`.
    #' @param dim Numeric. Dimensions (in pixels) for each plot when saving. Defaults to `2000`.
    #' @param width Numeric. Scaling factor for plot widths. Defaults to `1`.
    #' @param colors A character vector of colors for the timepoints. Defaults to the "Okabe-Ito" colorblind-friendly palette.
    #'
    #' @return A combined plot object (class `patchwork`) of all specified variables.
    #'
    #' @details
    #' - The function fits a linear mixed-effects model (`lmer`) for each variable with 
    #'   `exercise_timepoint` as a fixed effect and `subject_id` as a random effect.
    #' - Significance segments are drawn between timepoints based on the model's results.
    #' - Annotation text above segments can include p-values, significance stars, or effect sizes, 
    #'   depending on the `display` parameter.
    #' - If `file` is specified, the combined plot is saved as an image.
    #'
    #' @examples
    #' # Example usage:
    #' data <- data.frame(
    #'   subject_id = rep(1:10, each = 2),
    #'   exercise_timepoint = rep(c("pre", "post"), times = 10),
    #'   value = rnorm(20)
    #' )
    #' deltaPlot(
    #'   data = data,
    #'   vars = c("value"),
    #'   titles = c("Pre vs Post"),
    #'   ytitles = c("Value"),
    #'   xtitles = c("Timepoint"),
    #'   addpoints = TRUE,
    #'   plotsegments = "sign"
    #' )
    #'
    #' @import ggplot2
    #' @import lme4
    #' @import dplyr
    #' @import patchwork
    #' @export
    plots_list <- list()
    for (va in seq_along(vars)) {
        var <- vars[va]

        # Subset and prepare data
        data_sub <- data %>%
            ungroup() %>%
            dplyr::select(subject_id, exercise_timepoint, !!sym(var)) %>%
            mutate(exercise_timepoint = factor(dplyr::recode(exercise_timepoint, pre = "Pre", post = "Post"), levels = c("Pre", "Post"))) %>%
            filter(!(subject_id %in% dropouts))

        # Calculate max and min scores for y-axis limits
        max_score <- max(data_sub[[var]], na.rm = TRUE)
        min_score <- min(data_sub[[var]], na.rm = TRUE)
        delta <- max_score - min_score
        max_score <- max_score + delta * 0.1
        min_score <- min_score - delta * 0.1

        # Fit mixed model
        model <- lmer(pull(data_sub, var) ~ (1 | subject_id) + exercise_timepoint, data = data_sub)
        summary_model <- summary(model)
        p_value <- summary_model$coefficients[2, "Pr(>|t|)"]
        effect_size <- summary_model$coefficients[2, "Estimate"]

        # Generate annotation text
        annotation_text <- ""
        if ("effect" %in% display) {
            annotation_text <- paste0(annotation_text, "\u0394 = ", format(round(effect_size, 2), nsmall = 2))
        }
        if ("pval" %in% display) {
            annotation_text <- paste0(annotation_text, ", P = ", format(round(p_value, 3), nsmall = 3))
        }
        if ("pstar" %in% display) {
            annotation_text <- paste0(annotation_text, " ", pvalue_to_stars(p_value))
        }

        # Decide whether to plot segments based on significance
        segment_data <- NULL
        if (plotsegments == "all" || (plotsegments == "sign" && p_value < 0.05)) {
            segment_data <- data.frame(
                x = 1,
                xend = 2,
                y = max_score + delta * 0.05,
                yend = max_score + delta * 0.05  # Match height for horizontal segment
            )
        }

        # Generate plot
        plot <- ggplot(data_sub, aes(x = exercise_timepoint, y = !!sym(var), color = exercise_timepoint, fill = exercise_timepoint)) +
            mytheme +
            geom_violin(trim = TRUE, linewidth = 0.35 / width, alpha = 0.1, width = 1 / width) +
            geom_boxplot(width = 0.3 / width, linewidth = 0.35 / width, alpha = 0.1, outlier.size = 1) +
            geom_line(aes(group = subject_id), color = "grey", linewidth = 0.15) +
            geom_point(aes(shape = "Mean"), stat = "summary", fun = mean, size = 1, color = "black") +
            geom_line(aes(group = 1), stat = "summary", fun = mean, color = "black", linewidth = 0.25, show.legend = FALSE) + # Mean connecting line
            labs(
                title = if (length(titles) == 0) NULL else bquote(.(titles[[va]])), 
                y = if (length(ytitles) == 0) NULL else bquote(.(ytitles[[va]])),
                x = if (length(xtitles) == 0) NULL else bquote(.(xtitles[[va]])),
                color = "Timepoint",
                fill = "Timepoint"
            ) +
            scale_color_manual(values = colors[c(2, 7)]) +
            scale_fill_manual(values = colors[c(2, 7)])

        if (addpoints) {
            plot <- plot + geom_point(position = position_jitter(width = 0.2, seed = 24145))
        }

        if (!legend) {
            plot <- plot + theme(legend.position = "none")
        }

        # Plot the mean if required
        if ("mean" %in% display) {
            means <- data_sub %>% 
                group_by(exercise_timepoint) %>% 
                summarize(
                    mean_value = mean(!!sym(var), na.rm = TRUE),
                    y = max(!!sym(var), na.rm = TRUE)
                ) %>% 
                arrange(exercise_timepoint) %>% 
                mutate(
                    mean_value = round(mean_value, 2),
                    x = as.numeric(exercise_timepoint)  # Ensure x coordinates align with time points
                )
            
            plot <- plot +
                geom_text(data = means, aes(x = x, y = y, label = mean_value), vjust = -1, hjust = 0.5, size = 1.5, show.legend = FALSE)
            
            # set segment higher
            segment_data$y <- segment_data$y + delta*0.05

            
        }

        # Plot significance segment if required
        if (!is.null(segment_data)) {
            segment_data$y_end <- segment_data$y - 0.03 * delta
            plot <- plot <- plot + 
                    annotate("text", x = (segment_data$x + segment_data$xend) / 2, y = segment_data$y, label = annotation_text, hjust = 0.5, vjust = -0.5, size = 1.5, parse = FALSE) +
                    geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = y), inherit.aes = FALSE, color = "black", linewidth = 0.25) +
                    geom_segment(data = segment_data, aes(x = x, , y = y, yend = y_end), inherit.aes = FALSE, color = "black", linewidth = 0.25) +
                    geom_segment(data = segment_data, aes(x = xend, y = y, yend = y_end), inherit.aes = FALSE, color = "black", linewidth = 0.25)
        }

        # Adjust y-axis limits
        plot <- plot + ylim(min_score, max_score + 0.15*delta)

        if (removeyaxisticks) {
            plot <- plot + theme(axis.text.y = element_blank())
        }

        plots_list[[va]] <- plot
    }

    combined_plots <- wrap_plots(plots_list, ncol = ncol)

    if (collectguides) {
        combined_plots <- combined_plots + plot_layout(guides = "collect")
    }

    if (collectaxis) {
        combined_plots <- combined_plots + plot_layout(axis = "collect")
    }

    if (!is.null(file)) {
        ggsave(file, combined_plots, width = ncol * dim, height = nrow * dim, units = "px", bg = "white", limitsize = FALSE)
    }

    return(combined_plots)
}

correlationPlot <- function(
    data, 
    yvars, 
    xvars, 
    file = NULL, 
    sortplotsbyX = TRUE, 
    removexaxisticks = FALSE, 
    removeyaxisticks = FALSE, 
    ncol = 3,
    dim = 2000,
    titles = NULL, 
    xtitles = NULL, 
    ytitles = NULL, 
    labels = NULL, 
    labelxposition = "max", 
    labelyposition = "max", 
    color = as.vector(palette.colors(palette = "Okabe-Ito")[c("reddishpurple", "vermillion", "orange", "yellow", "bluishgreen", "skyblue", "blue", "black", "gray")])
) {
    #' correlationPlot: Create Correlation Plots with Automatic R² and p-value Calculation
    #'
    #' @description
    #' The `correlationPlot` function generates scatterplots for specified variable pairs, overlays a 
    #' linear regression line, and annotates the plots with R² and p-values. The function automatically 
    #' calculates R² and p-values for all combinations of specified variables.
    #'
    #' @param data A data frame containing the variables to plot.
    #' @param yvars A character vector of dependent variable names (y-axis).
    #' @param xvars A character vector of independent variable names (x-axis).
    #' @param file A string specifying the output file path for saving the combined plot. Defaults to `NULL` (no saving).
    #' @param sortplotsbyX Logical. If `TRUE`, the plots are sorted by `xvars`; otherwise, by `yvars`. Defaults to `TRUE`.
    #' @param removexaxisticks Logical. If `TRUE`, removes x-axis tick marks from the plots. Defaults to `FALSE`.
    #' @param removeyaxisticks Logical. If `TRUE`, removes y-axis tick marks from the plots. Defaults to `FALSE`.
    #' @param ncol Integer. Number of columns in the combined plot layout. Defaults to `3`.
    #' @param titles A character vector of plot titles. If not specified, titles are left blank.
    #' @param xtitles A character vector of x-axis titles. If only one is provided, it is replicated across plots.
    #' @param ytitles A character vector of y-axis titles. If only one is provided, it is replicated across plots.
    #' @param labels A character vector of custom labels to display on the plots. Defaults to `NULL`, in which case 
    #'   labels are automatically generated from the R² and p-values.
    #' @param labelxposition Character. Specifies the x-coordinate for labels. Options are `"max"` (default) 
    #'   or `"min"`. Determines whether the label is positioned at the maximum or minimum x-value.
    #' @param labelyposition Character. Specifies the y-coordinate for labels. Options are `"max"` (default) 
    #'   or `"min"`. Determines whether the label is positioned at the maximum or minimum y-value.
    #' @param color A character vector of colors for points in the plots. Defaults to the "Okabe-Ito" colorblind-friendly palette.
    #'
    #' @return A combined plot object (class `patchwork`) of all specified variable pairs.
    #'
    #' @details
    #' - The function calculates R² and p-values using linear regression for each variable pair and generates scatterplots.
    #' - Annotations display R² values and p-values, formatted automatically based on their magnitude.
    #' - Plots are arranged in a grid layout, with customizable titles, labels, and axis tick marks.
    #'
    #' @examples
    #' # Example usage:
    #' data <- data.frame(
    #'   var1 = rnorm(100),
    #'   var2 = rnorm(100),
    #'   var3 = rnorm(100)
    #' )
    #' correlationPlot(
    #'   data = data,
    #'   yvars = c("var1", "var2"),
    #'   xvars = c("var2", "var3"),
    #'   titles = c("Var1 vs Var2", "Var2 vs Var3"),
    #'   xtitles = c("Independent Variable"),
    #'   ytitles = c("Dependent Variable"),
    #'   ncol = 2
    #' )
    #'
    #' @import ggplot2
    #' @import patchwork
    #' @import dplyr
    #' @export
    # Internal function to perform linear regression and extract summary
    lmTest <- function(data, x, y) {
        formula <- as.formula(paste(y, "~", x))
        res <- summary(lm(formula, data = data))
        return(res)
    }

    # Ensure labelxposition and labelyposition are character vectors
    if (length(labelxposition) == 1) {
        labelxposition <- rep(as.character(labelxposition), length(xvars) * length(yvars))
    }
    if (length(labelyposition) == 1) {
        labelyposition <- rep(as.character(labelyposition), length(yvars) * length(xvars))
    }

    # Calculate R² and p-values for all combinations of yvars and xvars
    R2_list <- list()
    for (yvar in yvars) {
        R2_list[[yvar]] <- list()
        for (xvar in xvars) {
            R2_list[[yvar]][[xvar]] <- lmTest(data, xvar, yvar)
        }
    }

    # Create matrices for R² and p-values
    R2_matrix <- matrix(
        unlist(lapply(R2_list, function(x) sapply(x, function(y) y$r.squared))),
        ncol = length(yvars),
        dimnames = list(xvars, yvars)
    )
    p_values_matrix <- matrix(
        unlist(lapply(R2_list, function(x) sapply(x, function(y) y$coefficients[2, 4]))),
        ncol = length(yvars),
        dimnames = list(xvars, yvars)
    )

    data_R2 <- list(R2 = R2_matrix, p = p_values_matrix)

    # Generate plots
    plots_list <- list()
    if (length(labelxposition) == 1) {
        labelxposition <- rep(labelxposition, length(xvars) * length(yvars))
    }
    if (length(labelyposition) == 1) {
        labelyposition <- rep(labelyposition, length(yvars) * length(xvars))
    }
    if (!is.null(xtitles) && length(xtitles) == 1) {
        xtitles <- rep(xtitles, length(xvars) * length(yvars))
    }
    if (!is.null(ytitles) && length(ytitles) == 1) {
        ytitles <- rep(ytitles, length(xvars) * length(yvars))
    }
    if (sortplotsbyX) {
        combinations <- expand.grid(yvar = yvars, xvar = xvars)
    } else {
        combinations <- expand.grid(yvar = yvars, xvar = xvars) %>% arrange(yvar, xvar)
    }

    for (i in 1:nrow(combinations)) {
        yvar <- as.character(combinations$yvar[i])
        xvar <- as.character(combinations$xvar[i])
        data_sub <- data[complete.cases(data[c(xvar, yvar)]), c(xvar, yvar)]
        
        label_x <- as.numeric(if (labelxposition[i] == "max") max(data_sub[[xvar]], na.rm = TRUE) else min(data_sub[[xvar]], na.rm = TRUE))
        label_y <- as.numeric(if (labelyposition[i] == "max") max(data_sub[[yvar]], na.rm = TRUE) else min(data_sub[[yvar]], na.rm = TRUE))

        plot <- ggplot(data_sub, aes(x = !!sym(xvar), y = !!sym(yvar))) +
            mytheme +
            labs(
                title = bquote(.(titles[[i]])),
                x = bquote(.(xtitles[[i]])),
                y = bquote(.(ytitles[[i]]))
            ) +
            geom_smooth(method = "lm", se = FALSE, color = "darkgray", linetype = "dashed", linewidth = 0.35) +
            geom_point(size = 0.6, color = color[i])

        if (!is.null(labels)) {
            plot <- plot +
                annotate(
                    "text",
                    x = label_x,
                    y = label_y,
                    label = bquote(.(labels[i])),
                    hjust = ifelse(labelxposition[i] == "max", 1, 0),
                    vjust = 0.5,
                    color = "black",
                    size = 1.5
                )
        } else {
            plot <- plot +
                annotate(
                    "text",
                    x = label_x,
                    y = label_y,
                    label = {
                        p_value <- data_R2$p[xvar, yvar]
                        if (p_value < 0.001) {
                            p_formatted <- formatC(p_value, format = "e", digits = 7)
                            parts <- strsplit(p_formatted, "e")[[1]]
                            p1 <- as.numeric(parts[1])
                            p2 <- as.integer(parts[2])
                            bquote(R^2 == .(format(round(data_R2$R2[xvar, yvar], 2), nsmall = 2)) * "," ~ P == .(format(round(p1, 3), nsmall = 3)) %*% 10^.(p2))
                        } else {
                            bquote(R^2 == .(format(round(data_R2$R2[xvar, yvar], 2), nsmall = 2)) * "," ~ P == .(format(round(p_value, 3), nsmall = 3)))
                        }
                    },
                    hjust = ifelse(labelxposition[i] == "max", 1, 0),
                    vjust = ifelse(labelyposition[i] == "max", 1, 0),
                    color = "black",
                    size = 1.5
                )
        }

        if (removexaxisticks) {
            plot <- plot + theme(
                axis.text.x = element_blank()
            )
        }
        if (removeyaxisticks) {
            plot <- plot + theme(
                axis.text.y = element_blank()
            )
        }

        plots_list[[i]] <- plot
    }
    combined_plots <- wrap_plots(plots_list, ncol = ncol)

    if (!is.null(file)) {
        ggsave(file, combined_plots, width = ncol * dim, height = nrow * dim, units = "px", bg = "white", limitsize = FALSE)
    }

    return(combined_plots)
}


##############################
########    Figures    #######
##############################

######## EA ifo CA

pheno_plot <- pheno %>%
    mutate(exercise_timepoint = recode(exercise_timepoint, "pre" = "Pre", "post" = "Post")) %>%
    mutate(exercise_timepoint = factor(exercise_timepoint, levels = c("Pre", "Post")))
pos <- min(pheno[, "age"])

## GrimAge

clock <- "epigenetic_age_GrimAge"

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
    geom_smooth(method = "lm", se = FALSE, color = "darkgray", linetype = "dashed", linewidth = 0.35) +
    geom_abline(slope = 1, intercept = 0, color = "darkgray", linetype = "dotted", linewidth = 0.5) +
    geom_path(aes(group = subject_id), linewidth = 0.2, color = "black") +
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
    geom_point(aes(color = exercise_timepoint), size = 0.6) +
    scale_color_manual(
        name = "Timepoint",
        values = c("Pre" = colors[2], "Post" = colors[7])
    ) +
    mytheme +
    annotate( # trend line equation
        "text",
        x = max(pheno[, "age"]), y = min(pheno[, "age"]),
        label = paste("EA == ", format(round(coef(lm(pheno[, clock] ~ pheno[, "age"]))[2], 2), nsmall = 2), "%*% CA +", format(round(coef(lm(pheno[, clock] ~ pheno[, "age"]))[1], 2), nsmall = 2)),
        color = "black", hjust = 1, size = 1.5, parse = TRUE
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
        label = "R^2 == 0.84 * ',' ~ P < 2.092 %*% 10^{-20}", # mixed model anova & partial R2
        hjust = 0, vjust = 1, color = "black", size = 1.5, parse = TRUE
    ) +
    labs(
        x = bquote(.(expression(CA ~ (Years)))),
        y = bquote(.(expression(GrimAge ~ EA ~ (Years))))
    ) +
    geom_point(data = data.frame(), aes(x = pos, y = pos), color = "darkgray", size = 0.6) +
    geom_point(data = data.frame(), aes(x = pos + 0.5, y = pos + 0.5), color = "darkgray", size = 0.6) +
    geom_line(data = data.frame(), aes(x = c(pos, pos + 0.5), y = c(pos, pos + 0.5)), color = "darkgray", linewidth = 0.3) +
    theme(legend.position = "none")

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
    geom_smooth(method = "lm", se = FALSE, color = "darkgray", linetype = "dashed", linewidth = 0.35) +
    geom_abline(slope = 1, intercept = 0, color = "darkgray", linetype = "dotted", linewidth = 0.5) +
    geom_path(aes(group = subject_id), linewidth = 0.2, color = "black") +
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
    geom_point(aes(color = exercise_timepoint), size = 0.6) +
    scale_color_manual(
        name = "Timepoint",
        values = c("Pre" = colors[2], "Post" = colors[7])
    ) +
    mytheme +
    annotate( # trend line equation
        "text",
        x = max(pheno[, "age"]), y = min(pheno[, "age"]),
        label = paste("EA == ", format(round(coef(lm(pheno[, clock] ~ pheno[, "age"]))[2], 2), nsmall = 2), "%*% CA +", format(round(coef(lm(pheno[, clock] ~ pheno[, "age"]))[1], 2), nsmall = 2)),
        color = "black", hjust = 1, size = 1.5, parse = TRUE
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
        label = "R^2 == 0.71 * ',' ~ P == 1.333 %*% 10^{-13}",
        hjust = 0, vjust = 1, color = "black", size = 1.5, parse = TRUE
    ) +
    labs(
        x = bquote(.(expression(CA ~ (Years)))),
        y = bquote(.(expression(Horvath ~ EA ~ (Years))))
    ) +
    geom_point(data = data.frame(), aes(x = pos, y = pos), color = "darkgray", size = 0.6) +
    geom_point(data = data.frame(), aes(x = pos + 0.5, y = pos + 0.5), color = "darkgray", size = 0.6) +
    geom_line(data = data.frame(), aes(x = c(pos, pos + 0.5), y = c(pos, pos + 0.5)), color = "darkgray", linewidth = 0.3) +
    theme(legend.position = "none")

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
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 9))
ggsave("Images/Paper/Horvath_EA.png", combined_plot_horvath, dpi = 1200, width = 16, height = 6, units = "cm")

combined_plot_grimage <- (plot_grimage | dplot_grimage_acc) +
    plot_layout(widths = c(2, 1), guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 9))
ggsave("Images/Paper/GrimAge_EA.png", combined_plot_grimage, dpi = 1200, width = 16, height = 6, units = "cm")


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
plotdata <- d_bio
plotdata[,c("estimated_neutrophils")] <- d_bio[,c("estimated_neutrophils")]*100

Neu_plot <- correlationPlot(
    data = plotdata, 
    yvars = "epigenetic_age_acceleration_GrimAge", 
    xvars = "estimated_neutrophils", 
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
ggsave("Images/Paper/GrimAge_Correlations.png", combined_plot_corr, dpi = 1200, width = 16, height = 6, units = "cm")


######### (D) EAA correlations

plot_data <- bio %>%
    mutate(exercise_timepoint = factor(recode(exercise_timepoint, "pre" = "Pre", "post" = "Post"), levels = c("Pre", "Post"))) # set order right

EAA_plot <- ggplot(plot_data, aes(x = epigenetic_age_acceleration_Horvath, y = epigenetic_age_acceleration_GrimAge, color = exercise_timepoint)) +
    mytheme +
    labs(
        x = expression(Horvath ~ EAA ~ (Years)),
        y = expression(GrimAge ~ EAA ~ (Years))
    ) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgray", linetype = "dashed", linewidth = 0.35) +
    geom_path(aes(group = subject_id), color = "black", linewidth = 0.25) +
    geom_point(size = 0.6) +
    scale_color_manual(values = c("Pre" = colors[2], "Post" = colors[7])) +
    theme(legend.position = "none")

D_EAA_plot <- correlationPlot(
    data = d_bio, 
    yvars = "epigenetic_age_acceleration_GrimAge", 
    xvars = "epigenetic_age_acceleration_Horvath", 
    ncol = 1, 
    xtitles = expression(Horvath ~ Delta * EAA ~ (Years)),
    ytitles = expression(GrimAge ~ Delta * EAA ~ (Years)),
    labelxposition = "min",
    labelyposition = "max",
    color = colors[3]
)

combined_plot <- (EAA_plot | D_EAA_plot) +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 9))
ggsave("Images/Paper/Clocks_Association.png", combined_plot, dpi = 1200, width = 16, height = 6, units = "cm")

########## Delta Epigenetic age acceleration & bloodcells correlation

bloodcells <- c(
    "estimated_neutrophils",
    "estimated_CD4T_cells",
    "estimated_CD8T_cells",
    "estimated_natural_killer_cells",
    "estimated_monocytes",
    "estimated_B_cells"
)

bloodcell_names <- c(
    expression(Neutrophils),
    expression(CD4^{"+"} ~ T ~ Lymphocytes),
    expression(CD8^{"+"} ~ T ~ Lymphocytes),
    expression(Natural ~ Killer ~ Cells),
    expression(Monocytes),
    expression(B ~ Lymphocytes)
)

xtitles <- as.vector(lapply(bloodcell_names, function(name) {
    bquote(Delta * .(name) ~ ("%"))
}))

plotdata <- d_bio
plotdata[,bloodcells] <- d_bio[,bloodcells]

BC_grimage_plot <- correlationPlot(
    data = plotdata,
    yvars = "epigenetic_age_acceleration_GrimAge",
    xvars = bloodcells, 
    ncol = 3, 
    xtitles = xtitles,
    ytitles = rep(expression(GrimAge ~ Delta * EAA ~ (Years)), length(bloodcell_names)),
    color = colors[c(3, 4, 5, 6, 1, 2)],
    sortplotsbyX = FALSE
)

BC_horvath_plot <- correlationPlot(
    data = plotdata,
    yvars = "epigenetic_age_acceleration_Horvath",
    xvars = bloodcells, 
    ncol = 3, 
    xtitles = xtitles,
    ytitles = rep(expression(Horvath ~ Delta * EAA ~ (Years)), length(bloodcell_names)),
    color = colors[c(3, 4, 5, 6, 1, 2)],
    sortplotsbyX = FALSE
)

combined_plot <- BC_grimage_plot / BC_horvath_plot +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 9)
)
ggsave("Images/Paper/Clocks_Bloodcells_Associations.png", combined_plot, dpi = 1200, width = 16, height = 18, units = "cm")

######### Mean leukocyt composition pie chart

leukocytes <- bio %>%
    filter(!dropout, !excluded) %>%
    ungroup() %>%
    dplyr::select(exercise_timepoint, contains("estimated_")) %>%
    rename(
        "CD4+ T Lymphocytes" = estimated_CD4T_cells,
        "CD8+ T Lymphocytes" = estimated_CD8T_cells,
        "Natural Killer Cells" = estimated_natural_killer_cells,
        "B Lymphocytes" = estimated_B_cells,
        "Monocytes" = estimated_monocytes,
        "Neutrophils" = estimated_neutrophils
    ) %>%
    filter(exercise_timepoint %in% c("pre", "post")) %>%
    mutate(exercise_timepoint = factor(recode(exercise_timepoint, "pre" = "Pre", "post" = "Post"), levels = c("Pre", "Post"))) %>%
    pivot_longer(cols = -exercise_timepoint, names_to = "Leukocyt", values_to = "Concentration") %>%
    group_by(exercise_timepoint, Leukocyt) %>%
    summarize(across(everything(), mean, na.rm = TRUE)) %>%
    mutate(Concentration = (Concentration / sum(Concentration)) * 100)

leukocytes_pre <- leukocytes %>%
    filter(exercise_timepoint == "Pre") %>%
    ungroup() %>%
    dplyr::select(-exercise_timepoint) %>%
    dplyr::arrange(-Concentration) %>%
    mutate(Leukocyt = factor(Leukocyt, levels = Leukocyt))

leukocytes_post <- leukocytes %>%
    filter(exercise_timepoint == "Post") %>%
    ungroup() %>%
    dplyr::select(-exercise_timepoint) %>%
    dplyr::arrange(-Concentration) %>%
    mutate(Leukocyt = factor(Leukocyt, levels = Leukocyt))

# Pie charts

shared_fill_scale <- scale_fill_manual(values = colors[c(3, 4, 5, 6, 1, 2)])

piechart_pre <- ggplot(leukocytes_pre, aes(x = "", y = Concentration, fill = Leukocyt)) +
    geom_bar(stat = "identity", width = 1, alpha = 0.8) +
    coord_polar("y", start = 0) +
    shared_fill_scale +  
    theme_void() +
    theme(
        legend.position = "none",
        plot.title = element_text(margin = margin(b = 8, t = 8), hjust = 0.5, size = 8)
    ) +  
    geom_text(
        aes(label = paste0(round(Concentration, 2), "%")), 
        position = position_stack(vjust = 0.5),  # Center the text in each slice
        size = 2
    ) + 
    labs(title = "Pre-training")

piechart_post <- ggplot(leukocytes_post, aes(x = "", y = Concentration, fill = Leukocyt)) +
    geom_bar(stat = "identity", width = 1, alpha = 0.8) +
    coord_polar("y", start = 0) +
    shared_fill_scale +  
    theme_void() +
    geom_text(
        aes(label = paste0(round(Concentration, 2), "%")), 
        position = position_stack(vjust = 0.5),  # Center the text in each slice
        size = 2
    ) +
    theme(
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.spacing.x = unit(0.50, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'),
        plot.title = element_text(margin = margin(b = 8, t = 8), hjust = 0.5, size = 8)
    ) + 
    labs(title = "Post-training")

# Deltaplots

bloodcell_names <- c(
    expression(Neutrophils),
    expression(CD4^{"+"} ~ T ~ Lymphocytes),
    expression(CD8^{"+"} ~ T ~ Lymphocytes),
    expression(Natural ~ Killer ~ Cells),
    expression(Monocytes),
    expression(B ~ Lymphocytes)
)

plotdata <- bio
plotdata[,c("estimated_B_cells", "estimated_CD4T_cells", "estimated_CD8T_cells", "estimated_monocytes", "estimated_neutrophils", "estimated_natural_killer_cells")] <- bio[,c("estimated_B_cells", "estimated_CD4T_cells", "estimated_CD8T_cells", "estimated_monocytes", "estimated_neutrophils", "estimated_natural_killer_cells")]

ytitles <- sapply(bloodcell_names, function(name) {
    bquote(Delta * .(name) ~ ("%"))
})

BCs <- deltaPlot(
    data = plotdata, 
    vars = c("estimated_neutrophils", "estimated_CD4T_cells", "estimated_CD8T_cells", "estimated_natural_killer_cells", "estimated_monocytes", "estimated_B_cells"), 
    ncol = 3, 
    dropouts = bio %>% filter(dropout | excluded) %>% pull(subject_id),
    plotsegments = "all",
    display =  c("mean", "effect", "pval"),
    ytitles = ytitles,
    legend = FALSE,
    collectaxis = FALSE,
    width = 1.1
)

# add spacers to center top and bottom figure
combined_piecharts <- (plot_spacer() | piechart_pre | piechart_post | plot_spacer()) +
    plot_layout(widths = c(0.3, 1, 1, 0.3)) +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 9))

combined_leukocytes <- combined_piecharts / BCs + 
    plot_layout(heights = c(1, 2)) +
    plot_annotation(tag_levels = "A")  & 
    theme(plot.tag = element_text(size = 9))
ggsave("Images/Paper/Leukocyte_Composition.png", combined_leukocytes, dpi = 1200, width = 16, height = 16, units = "cm", bg = "white")


######### Training time

trainingtime <- read_xlsx("Data/Pheno/Training_Time.xlsx")

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
        linewidth = 0.25
    ) +
    # adherence is scaled by 12 to match the hours scale
    geom_errorbar( 
        data = filter(trainingtime_errors, variable == "adherence"),
        aes(x = week, ymin = (mean_value - error)/12, ymax = (mean_value + error)/12, color = variable),
        width = 0.2,
        linewidth = 0.25
    ) +
    geom_point(
        data = filter(trainingtime_means, variable %in% c("executed", "planned")),
        aes(x = week, y = mean_value, color = variable), 
        size = 0.5
    ) +
    geom_line(
        data = filter(trainingtime_means, variable %in% c("executed", "planned")),
        aes(x = week, y = mean_value, color = variable),
        linewidth = 0.25
    ) +
    geom_point(
        data = filter(trainingtime_means, variable == "adherence"),
        aes(x = week, y = mean_value / 12, color = variable), 
        size = 0.25
    ) +
    geom_line(
        data = filter(trainingtime_means, variable == "adherence"),
        aes(x = week, y = mean_value / 12, color = variable),
        linewidth = 0.25
    ) +
    geom_hline(yintercept = 100 / 12, linetype = "dotted", color = colors[7], linewidth = 0.25) +
    scale_y_continuous(sec.axis = sec_axis(~ . * 12, name = bquote(.(expression(Adherence ~ ("%")))), breaks = c(0, 50, 100, 150))) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
    guides(color = guide_legend(title = NULL)) +
    scale_color_manual(
        values = c("executed" = colors[5], "planned" = colors[3], "adherence" = colors[7]),
        breaks = c("Planned", "Executed", "Adherence")
    ) +
    labs(
        y = bquote(.(expression(Planned ~ "&" ~ Executed ~ Training ~ (Hours)))),
        x = bquote(.(expression(Week)))
    ) +
    mytheme
ggsave("Images/Paper/Training_Time.png", plot, dpi = 1200, width = 16, height = 6, units = "cm", bg = "white")
