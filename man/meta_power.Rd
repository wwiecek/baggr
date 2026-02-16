% Generated manually
\name{meta_power}
\alias{meta_power}
\title{Power surface for random-effects meta-analysis designs}
\usage{
meta_power(
  x,
  max_mu = NULL,
  max_tau = NULL,
  n = 10,
  alpha = 0.05,
  sided = 2,
  contours = c(0.8),
  print_plot = TRUE
)
}
\arguments{
\item{x}{Input containing study-level standard errors. One of: numeric vector,
data frame with \code{se} column, or a \code{baggr} object.}

\item{max_mu}{Maximum value of \code{mu} in the grid. If \code{NULL}, defaults to
\code{2 * max(se)} (or \code{2 * hypermean(x)} when available for baggr objects).}

\item{max_tau}{Maximum value of \code{tau} in the grid. If \code{NULL}, defaults to
\code{2 * max(se)} (or \code{2 * hypersd(x)} when available for baggr objects).}

\item{n}{Number of grid points for each axis (\code{mu} and \code{tau}).}

\item{alpha}{Significance level for the z-test.}

\item{sided}{Number of test sides: \code{1} (one-sided) or \code{2} (two-sided).}

\item{contours}{Numeric vector of power levels used to draw contour lines.
Use \code{c()} to suppress contours.}

\item{print_plot}{Logical; if \code{TRUE}, prints the plot.}
}
\value{
An invisible list with:
\itemize{
\item \code{plot}: a \code{ggplot2} object,
\item \code{values$grid_wide}: one row per \code{(mu, tau)} with \code{power_RE},
\item \code{values$grid_long}: long-format version used for plotting.
}
}
\description{
Computes frequentist power of a random-effects meta-analysis z-test,
over a grid of true mean effects (\code{mu}) and between-study heterogeneity
(\code{tau}), based on study-level standard errors.
}
\details{
For each grid point \code{(mu, tau)}, power is computed for a random-effects
z-test using inverse-variance weighting with study variance terms
\code{se_k^2 + tau^2} (i.e. treating heterogeneity as known).

This assumes known heterogeneity (\code{tau}), which is never literally true in
practice. Therefore this assessment should be treated as a supporting,
approximate diagnostic and used alongside simulation-based power analysis.

You can supply a numeric vector of standard errors, a data frame with a
\code{se} column, or a fitted \code{baggr} object (from which \code{se} values are extracted).

If a \code{baggr} object is supplied and hyper-parameter summaries are available,
the plot includes a highlighted point at \code{(hypermean, hypersd)} to show
where the fitted model sits on the power surface.
}
\examples{
# 8 schools example
data(schools)
out <- meta_power(schools, n = 25, contours = c(0.5, 0.8))
head(out$values$grid_wide)

# You can also pass a baggr object directly.
\dontrun{
bg <- baggr(schools, iter = 1000)
meta_power(bg, n = 25)
}

# No contours:
meta_power(schools, n = 20, contours = c())
}
