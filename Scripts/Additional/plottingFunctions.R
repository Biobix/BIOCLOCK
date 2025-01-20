
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
        axis.line = element_line(colour = "black", size = 0.3)
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

# deltaPlot: Create Violin and Boxplots for Paired Data with Optional Significance Segments
#
# @description
# The `deltaPlot` function generates paired violin and boxplots for longitudinal or paired data, 
# highlighting differences between two timepoints (e.g., Pre vs Post). It includes options to 
# add significance segments with customizable annotations, and can combine multiple plots into 
# a grid layout.
#
# @param data A data frame containing the data to be plotted.
# @param vars A character vector of variable names to plot. Each variable is plotted separately.
# @param dropouts A character vector of `subject_id` values to exclude from the analysis. 
#   Defaults to subjects who are marked as `dropout` or `excluded`.
# @param file A string specifying the output file path for saving the combined plot. Defaults to `NULL` (no saving).
# @param titles A character vector of titles for the plots. Defaults to an empty vector.
# @param ytitles A character vector of y-axis titles for the plots. Defaults to an empty vector.
# @param xtitles A character vector of x-axis titles for the plots. Defaults to an empty vector.
# @param addpoints Logical. If `TRUE`, adds individual data points to the plots. Defaults to `FALSE`.
# @param removeyaxisticks Logical. If `TRUE`, removes y-axis tick marks from the plots. Defaults to `FALSE`.
# @param collectguides Logical. If `TRUE`, collects the guides (e.g., legends) for all subplots. Defaults to `TRUE`.
# @param collectaxis Logical. If `TRUE`, aligns the axes across subplots. Defaults to `FALSE`.
# @param legend Logical. If `TRUE`, includes the legend in each subplot. Defaults to `FALSE`.
# @param plotsegments Character. Determines whether to plot significance segments. Options are:
#   - `"all"`: Always plot segments.
#   - `"sign"`: Only plot segments if the p-value is significant.
#   - `"no"`: Do not plot any segments.
#   Defaults to `"all"`.
# @param display Character vector. Specifies what information to display above the segments. Options are:
#   - `"none"`: No text above the segments.
#   - `"pval"`: Display p-values.
#   - `"pstar"`: Display significance stars.
#   - `"effect"`: Display effect sizes (e.g., mean differences).
#   - `"mean"`: Display mean values of each timepoint.
#   Defaults to `c("none", "pval", "pstar", "effect", "mean")`.
# @param ncol Integer. Number of columns in the combined plot layout. Defaults to `1`.
# @param nrow Integer. Number of rows in the combined plot layout. Defaults to `ceiling(length(vars) / ncol)`.
# @param dim Numeric. Dimensions (in pixels) for each plot when saving. Defaults to `2000`.
# @param width Numeric. Scaling factor for plot widths. Defaults to `1`.
# @param colors A character vector of colors for the timepoints. Defaults to the "Okabe-Ito" colorblind-friendly palette.
#
# @return A combined plot object (class `patchwork`) of all specified variables.
#
# @details
# - The function fits a linear mixed-effects model (`lmer`) for each variable with 
#   `exercise_timepoint` as a fixed effect and `subject_id` as a random effect.
# - Significance segments are drawn between timepoints based on the model's results.
# - Annotation text above segments can include p-values, significance stars, or effect sizes, 
#   depending on the `display` parameter.
# - If `file` is specified, the combined plot is saved as an image.
#
# @examples
# # Example usage:
# data <- data.frame(
#   subject_id = rep(1:10, each = 2),
#   exercise_timepoint = rep(c("pre", "post"), times = 10),
#   value = rnorm(20)
# )
# deltaPlot(
#   data = data,
#   vars = c("value"),
#   titles = c("Pre vs Post"),
#   ytitles = c("Value"),
#   xtitles = c("Timepoint"),
#   addpoints = TRUE,
#   plotsegments = "sign"
# )
#

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
            geom_violin(trim = TRUE, size = 0.35 / width, alpha = 0.1, width = 1 / width) +
            geom_boxplot(width = 0.3 / width, size = 0.35 / width, alpha = 0.1, outlier.size = 1) +
            geom_line(aes(group = subject_id), color = "grey", size = 0.15) +
            geom_point(aes(shape = "Mean"), stat = "summary", fun = mean, size = 1, color = "black") +
            geom_line(aes(group = 1), stat = "summary", fun = mean, color = "black", size = 0.25) + # Mean connecting line
            labs(
                title = if (length(titles) == 0) NULL else bquote(.(titles[[va]])), 
                y = if (length(ytitles) == 0) NULL else bquote(.(ytitles[[va]])),
                x = if (length(xtitles) == 0) NULL else bquote(.(xtitles[[va]])),
                color = "Timepoint",
                shape = ""
            ) +
            guides(
                fill = "none"
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
                    geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = y), inherit.aes = FALSE, color = "black", size = 0.25) +
                    geom_segment(data = segment_data, aes(x = x, , y = y, xend = x, yend = y_end), inherit.aes = FALSE, color = "black", size = 0.25) +
                    geom_segment(data = segment_data, aes(x = xend, y = y, xend = xend, yend = y_end), inherit.aes = FALSE, color = "black", size = 0.25)
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

# correlationPlot: Create Correlation Plots with Automatic R² and p-value Calculation
#
# @description
# The `correlationPlot` function generates scatterplots for specified variable pairs, overlays a 
# linear regression line, and annotates the plots with R² and p-values. The function automatically 
# calculates R² and p-values for all combinations of specified variables.
#
# @param data A data frame containing the variables to plot.
# @param yvars A character vector of dependent variable names (y-axis).
# @param xvars A character vector of independent variable names (x-axis).
# @param file A string specifying the output file path for saving the combined plot. Defaults to `NULL` (no saving).
# @param sortplotsbyX Logical. If `TRUE`, the plots are sorted by `xvars`; otherwise, by `yvars`. Defaults to `TRUE`.
# @param removexaxisticks Logical. If `TRUE`, removes x-axis tick marks from the plots. Defaults to `FALSE`.
# @param removeyaxisticks Logical. If `TRUE`, removes y-axis tick marks from the plots. Defaults to `FALSE`.
# @param ncol Integer. Number of columns in the combined plot layout. Defaults to `3`.
# @param titles A character vector of plot titles. If not specified, titles are left blank.
# @param xtitles A character vector of x-axis titles. If only one is provided, it is replicated across plots.
# @param ytitles A character vector of y-axis titles. If only one is provided, it is replicated across plots.
# @param labels A character vector of custom labels to display on the plots. Defaults to `NULL`, in which case 
#   labels are automatically generated from the R² and p-values.
# @param labelxposition Character. Specifies the x-coordinate for labels. Options are `"max"` (default) 
#   or `"min"`. Determines whether the label is positioned at the maximum or minimum x-value.
# @param labelyposition Character. Specifies the y-coordinate for labels. Options are `"max"` (default) 
#   or `"min"`. Determines whether the label is positioned at the maximum or minimum y-value.
# @param color A character vector of colors for points in the plots. Defaults to the "Okabe-Ito" colorblind-friendly palette.
#
# @return A combined plot object (class `patchwork`) of all specified variable pairs.
#
# @details
# - The function calculates R² and p-values using linear regression for each variable pair and generates scatterplots.
# - Annotations display R² values and p-values, formatted automatically based on their magnitude.
# - Plots are arranged in a grid layout, with customizable titles, labels, and axis tick marks.
#
# @examples
# # Example usage:
# data <- data.frame(
#   var1 = rnorm(100),
#   var2 = rnorm(100),
#   var3 = rnorm(100)
# )
# correlationPlot(
#   data = data,
#   yvars = c("var1", "var2"),
#   xvars = c("var2", "var3"),
#   titles = c("Var1 vs Var2", "Var2 vs Var3"),
#   xtitles = c("Independent Variable"),
#   ytitles = c("Dependent Variable"),
#   ncol = 2
# )
#
# @import ggplot2
# @import patchwork
# @import dplyr
# @export

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
            geom_smooth(method = "lm", se = FALSE, color = "darkgray", linetype = "dashed", size = 0.35) +
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
