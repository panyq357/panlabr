library(ggplot2)

#' Plot bar, errorbar and jitter
#'
#' Plot bar, errorbar and jitter dots in one function call with proper default
#' settings.
#'
#' @param bar_args Arguments for barplot.
#' @param errorbar_args Arguments for errorbar.
#' @param jitter_args Arguments for jitter.
#'
#' @export
#'
geom_triad <- function(bar_args = NULL,
                       errorbar_args = NULL,
                       jitter_args = NULL) {
  args <- list(
    bar = list(
      geom = "col",
      fun = "mean",
      color = "black",
      na.rm = TRUE,
      width = 0.6
    ),
    errorbar = list(
      geom = "errorbar",
      fun.data = "mean_sdl",
      fun.args = list("mult" = 1),
      na.rm = TRUE,
      width = 0.4
    ),
    jitter = list(
      height = 0,
      width = 0.2,
      size = 1,
      alpha = 0.5,
      stroke = 0
    )
  )

  if (length(bar_args) > 0) {
    for (name in names(bar_args)) {
      args$bar[[name]] <- bar_args[[name]]
    }
  }
  if (length(errorbar_args) > 0) {
    for (name in names(errorbar_args)) {
      args$errorbar[[name]] <- errorbar_args[[name]]
    }
  }
  if (length(jitter_args) > 0) {
    for (name in names(jitter_args)) {
      args$jitter[[name]] <- jitter_args[[name]]
    }
  }
  return(list(
    do.call(stat_summary, args$bar),
    do.call(stat_summary, args$errorbar),
    do.call(geom_jitter, args$jitter)
  ))
}

#' @export
theme_simple <- function(y_top_mult = 0.2, y_bottom_mult = 0) {
  return(list(
    scale_y_continuous(expand = expansion(mult = c(y_bottom_mult, y_top_mult))),
    theme_classic(),
    theme(legend.position = "none")
  ))
}

#' @export
theme_x_label_45 <- function() {
  return(theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}

GeomTestWithRef <- ggproto("GeomTestWithRef", GeomText,
  required_aes = c("x", "y"),
  extra_params = c("na.rm", "ref_group", "label_y_fun", "test"),
  setup_data = function(data, params) {
    ref_y <- data$y[data$x == params$ref]

    new_label <- split(data$y, data$group) |> sapply(function(grp_y) {
      p <- get(params$test)(ref_y, grp_y)$p.value
      if (p < 0.001) {
        return("***")
      } else if (p < 0.01) {
        return("**")
      } else if (p < 0.05) {
        return("*")
      } else {
        return("")
      }
    })

    new_y <- params$label_y_fun(data)

    data$label <- new_label[match(data$group, names(new_label))]
    data$y <- new_y[match(data$group, names(new_y))]

    data <- dplyr::distinct(data)

    return(data)
  }
)

#' @export
label_y_fun_max_0.1 <- function(data) {
  split(data$y, data$group) |>
    sapply(function(x) max(x, na.rm = T) + max(data$y, na.rm = T) * 0.1)
}

#' @export
label_y_fun_mean_0.1 <- function(data) {
  split(data$y, data$group) |>
    sapply(function(x) mean(x, na.rm = T) + max(data$y, na.rm = T) * 0.1)
}


#' Add signif stars (test with ref)
#'
#' Test all groups with ref group and plot signif stars.
#'
#' @param label_y_fun A function for signif star y pos.
#' @param ref_group Index of reference group level.
#' @param test Name of function for run test.
#' @export
#'
geom_test_with_ref <- function(mapping = NULL, data = NULL,
                               stat = "identity", position = "identity",
                               na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
                               label_y_fun = label_y_fun_max_0.1,
                               ref_group = 1, test = "t.test", ...) {
  layer(
    stat = stat,
    data = data,
    mapping = mapping,
    geom = GeomTestWithRef,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      ref_group = ref_group,
      label_y_fun = label_y_fun,
      test = test,
      na.rm = na.rm,
      ...
    )
  )
}

GeomANOVALSDGroups <- ggproto("GeomANOVALSDGroups", GeomText,
  required_aes = c("x", "y"),
  extra_params = c("na.rm", "label_y_fun"),
  setup_data = function(data, params) {
    model <- aov(y ~ group, data = data)
    out <- agricolae::LSD.test(model, trt = "group")
    new_label <- setNames(out$groups$groups, row.names(out$groups))

    new_y <- params$label_y_fun(data)

    data$label <- new_label[match(data$group, names(new_label))]
    data$y <- new_y[match(data$group, names(new_y))]

    data <- dplyr::distinct(data)

    return(data)
  }
)

#' Add signif stars (ANOVA LSD)
#'
#' Perform ANOVA and LSD tests, then plot LSD letters.
#'
#' @param label_y_fun A function for signif star y pos.
#' @export
#'
geom_anova_lsd_groups <- function(mapping = NULL, data = NULL,
                                  stat = "identity", position = "identity",
                                  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
                                  label_y_fun = label_y_fun_max_0.1, ...) {
  layer(
    stat = stat,
    data = data,
    mapping = mapping,
    geom = GeomANOVALSDGroups,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      label_y_fun = label_y_fun,
      na.rm = na.rm,
      ...
    )
  )
}
