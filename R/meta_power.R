#' Power surface for meta-analysis designs
#'
#' Computes the frequentist power of fixed-effects and random-effects
#' meta-analysis z-tests over a grid of true mean effects (`mu`) and
#' between-study heterogeneity (`tau`), based on study-level standard errors.
#'
#' You can supply:
#' - a numeric vector of standard errors,
#' - a data frame with a `se` column, or
#' - a fitted [baggr] object (from which `se` values are extracted).
#'
#' If a [baggr] object is supplied and hyper-parameter summaries are available,
#' the plot includes a highlighted point at (`hypermean`, `hypersd`) to show
#' where the fitted model sits on the power surface.
#'
#' @param x Input containing study-level standard errors. One of: numeric vector,
#'   data frame with `se` column, or a [baggr] object.
#' @param max_mu Maximum value of `mu` in the grid. If `NULL`, defaults to
#'   `2 * max(se)` (or `2 * hypermean(x)` when available for baggr objects).
#' @param max_tau Maximum value of `tau` in the grid. If `NULL`, defaults to
#'   `2 * max(se)` (or `2 * hypersd(x)` when available for baggr objects).
#' @param n Number of grid points for each axis (`mu` and `tau`).
#' @param alpha Significance level for the z-test.
#' @param sided Number of test sides: `1` (one-sided) or `2` (two-sided).
#' @param contour_breaks Power levels used to draw contour lines.
#' @param add_contours Logical; if `TRUE`, overlays contour lines.
#' @param print_plot Logical; if `TRUE`, prints the plot.
#'
#' @return An invisible list with:
#' - `plot`: a `ggplot2` object,
#' - `values$grid_wide`: one row per `(mu, tau)` with `power_FE` and `power_RE`,
#' - `values$grid_long`: long-format version used for plotting.
#'
#' @examples
#' # 8 schools example
#' data(schools)
#' out <- meta_power(schools, n = 25, contour_breaks = c(0.5, 0.8))
#' head(out$values$grid_wide)
#'
#' # You can also pass a baggr object directly.
#' \dontrun{
#' bg <- baggr(schools, iter = 1000)
#' meta_power(bg, n = 25)
#' }
#'
#' @export
meta_power <- function(x,
                       max_mu = NULL,
                       max_tau = NULL,
                       n = 10,
                       alpha = 0.05,
                       sided = 2,
                       contour_breaks = c(0.5, 0.8),
                       add_contours = TRUE,
                       print_plot = TRUE) {

  stopifnot(n >= 2, sided %in% c(1, 2), alpha > 0, alpha < 1)

  extract_se <- function(x) {
    if (is.numeric(x)) return(as.numeric(x))

    if (is.data.frame(x)) {
      if (!("se" %in% names(x))) {
        stop("If x is a data frame it must have a column named `se`.")
      }
      return(as.numeric(x$se))
    }

    if (inherits(x, "baggr")) {
      if (!is.null(x$data) && !is.null(x$data$se)) {
        return(as.numeric(x$data$se))
      }

      if (!is.null(x$summary_data)) {
        sd <- x$summary_data
        if (is.data.frame(sd) && ("se" %in% names(sd))) return(as.numeric(sd$se))
        if (is.list(sd) && !is.null(sd$se)) return(as.numeric(sd$se))
      }

      stop("Could not find SEs in baggr object (expected bg$data$se or bg$summary_data$se).")
    }

    stop("x must be a numeric vector of SEs, a data frame with column `se`, or a baggr object.")
  }

  as_scalar <- function(z) {
    if (is.null(z) || inherits(z, "try-error")) return(NA_real_)
    zv <- as.numeric(z)
    if (length(zv) < 1) return(NA_real_)
    zv[1]
  }

  se <- extract_se(x)
  se <- se[is.finite(se)]
  if (length(se) < 2 || any(se <= 0)) {
    stop("SEs must be finite, positive, and length >= 2.")
  }

  K <- length(se)

  baggr_point <- NULL
  if (inherits(x, "baggr")) {
    hm <- try(baggr::hypermean(x)[["mean"]], silent = TRUE)
    hs <- try(baggr::hypersd(x)[["mean"]], silent = TRUE)
    hm <- as_scalar(hm)
    hs <- as_scalar(hs)
    if (is.finite(hm) && is.finite(hs) && hs >= 0) {
      baggr_point <- c(mu = hm, tau = hs)
    }
  }

  if (is.null(max_mu) || is.null(max_tau)) {
    if (!is.null(baggr_point)) {
      if (is.null(max_mu)) max_mu <- 2 * baggr_point[["mu"]]
      if (is.null(max_tau)) max_tau <- 2 * baggr_point[["tau"]]
    }

    if (is.null(max_mu)) max_mu <- 2 * max(se)
    if (is.null(max_tau)) max_tau <- 2 * max(se)
  }

  if (!is.finite(max_mu) || max_mu < 0) stop("max_mu must be finite and >= 0.")
  if (!is.finite(max_tau) || max_tau < 0) stop("max_tau must be finite and >= 0.")

  mu_grid <- seq(0, max_mu, length.out = n)
  tau_grid <- seq(0, max_tau, length.out = n)

  zcrit <- stats::qnorm(1 - alpha / sided)

  z_power <- function(ncp) {
    if (sided == 2) {
      stats::pnorm(-zcrit - ncp) + (1 - stats::pnorm(zcrit - ncp))
    } else {
      1 - stats::pnorm(zcrit - ncp)
    }
  }

  power_fe <- function(mu, tau) {
    w <- 1 / se^2
    denom <- sum(w)
    var_hat <- sum((w^2) * (se^2 + tau^2)) / (denom^2)
    z_power(mu / sqrt(var_hat))
  }

  power_re_known_tau <- function(mu, tau) {
    w <- 1 / (se^2 + tau^2)
    var_hat <- 1 / sum(w)
    z_power(mu / sqrt(var_hat))
  }

  grid <- expand.grid(mu = mu_grid, tau = tau_grid)
  grid$power_FE <- mapply(power_fe, grid$mu, grid$tau)
  grid$power_RE <- mapply(power_re_known_tau, grid$mu, grid$tau)

  grid_fe <- grid[, c("mu", "tau", "power_FE")]
  names(grid_fe)[3] <- "power"
  grid_fe$model <- "Fixed-effects z-test"

  grid_re <- grid[, c("mu", "tau", "power_RE")]
  names(grid_re)[3] <- "power"
  grid_re$model <- "Random-effects z-test (known tau)"

  grid_long <- rbind(grid_fe, grid_re)

  subtitle_txt <- paste0(
    "K = ", K,
    ", alpha = ", alpha,
    if (sided == 2) ", two-sided" else ", one-sided",
    "; SE range = [", signif(min(se), 3), ", ", signif(max(se), 3), "]"
  )

  p <- ggplot2::ggplot(grid_long, ggplot2::aes(x = mu, y = tau, fill = power)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::facet_wrap(~ model, nrow = 1) +
    ggplot2::scale_fill_gradient(limits = c(0, 1)) +
    ggplot2::labs(
      x = expression(mu ~ "(true mean effect)"),
      y = expression(tau ~ "(heterogeneity SD)"),
      fill = "Power",
      title = "Meta-analysis power over (mu, tau)",
      subtitle = subtitle_txt
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )

  if (add_contours && length(contour_breaks) > 0) {
    p <- p +
      ggplot2::geom_contour(
        ggplot2::aes(z = power),
        breaks = contour_breaks,
        color = "red",
        linewidth = 1.0
      )
  }

  if (!is.null(baggr_point)) {
    point_df <- data.frame(
      mu = baggr_point[["mu"]],
      tau = baggr_point[["tau"]],
      model = c("Fixed-effects z-test", "Random-effects z-test (known tau)"),
      stringsAsFactors = FALSE
    )

    p <- p +
      ggplot2::geom_point(
        data = point_df,
        mapping = ggplot2::aes(x = mu, y = tau),
        shape = 21,
        size = 5,
        stroke = 1.4,
        fill = "yellow",
        color = "black",
        inherit.aes = FALSE
      )
  }

  if (isTRUE(print_plot)) print(p)

  invisible(list(
    plot = p,
    values = list(
      grid_wide = grid,
      grid_long = grid_long
    )
  ))
}
