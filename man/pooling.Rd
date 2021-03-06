% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pooling_metrics.R
\name{pooling}
\alias{pooling}
\alias{heterogeneity}
\title{Pooling metrics for baggr}
\usage{
pooling(bg, type = c("groups", "total"), summary = TRUE)

heterogeneity(bg, summary = TRUE)
}
\arguments{
\item{bg}{output of a baggr() function}

\item{type}{In \code{pooling} calculation is done for each of the \code{"groups"}
(default) or for \code{"total"} hypereffect(s).
See Details section for how calculation is done.}

\item{summary}{logical; if \code{FALSE} a whole vector of pooling values is returned,
otherwise only the means and intervals}
}
\value{
Matrix with mean and intervals for chosen pooling metric,
each row corresponding to one meta-analysis group.
}
\description{
Compute statistics relating to \code{heterogeneity} (whole model) and
\code{pooling} (for each group) given a \link{baggr} meta-analysis model.
The statistics are the pooling metric by Gelman & Pardoe (2006) or its
complement, the \emph{I-squared} statistic.
}
\details{
Pooling statistic describes the extent to which group-level estimates of treatment
effect are "pooled" (or pulled!) toward average treatment effect in the meta-analysis model.
If \code{pooling = "none"} or "full" in \link{baggr}, then the returned values are always 0 or 1, respectively.
If \code{pooling = "partial"}, the value is somewhere between 0 and 1.

\strong{Formulae for the calculations below are provided in main package vignette.} See \verb{vignette("baggr").}

#' \strong{Estimate of pooling in a group}: this is the calculation done by \code{pooling()}
if \code{type = "groups"} (default).

In a partial pooling model (see \link{baggr}), group \emph{k} (e.g. study) has a treatment effect
estimate, with some SE around the real treatment effect (TE).
Each TE itself is distributed with mean and variance.

The quantity of interest is ratio of variability in \eqn{\tau} to total variability.
By convention, we subtract it from 1, to obtain a \emph{pooling metric} \emph{p}.

\deqn{p = 1 - (\sigma(\tau)^2 / (\sigma_(\tau)^2 + se_k^2))}
\itemize{
\item If \eqn{p < 0.5}, that means the variation across studies is higher than variation within studies.
\item Values close to 1 indicate nearly full pooling. Variation across studies dominates.
\item Values close to 0 -- no pooling. Variation within studies dominates.
}

Note that, since \eqn{\sigma_{\tau}^2} is a Bayesian parameter (rather than a single fixed value)
\emph{p} is also a parameter. It is typical for \emph{p} to have very high dispersion, as in many cases we
cannot precisely estimate \eqn{\sigma_{\tau}}. To obtain the whole distribution of_p_
(rather than summarised values), set \code{summary=FALSE}.

\strong{Overall pooling (in the model)}

Typically it is a single measure of heterogeneity that is of interest to researchers.
This is calculated by setting \code{type = "total"} or simply writing \code{heterogeneity(mymodel)}

In many contexts, i.e. medical statistics, it is typical to report \emph{1-P}, called \eqn{I^2}
(see Higgins \emph{et al}, 2003). Higher values of \emph{I-squared} indicate higher heterogeneity.
Von Hippel (2015) provides useful details for \emph{I-squared} calculations.

Same as for group-specific estimates, \emph{P} is a Bayesian parameter and its dispersion can be high.

\strong{Relationship to R-squared statistic}

See Gelman & Pardoe (2006) Section 1.1 for a short explanation of how \eqn{R^2}
statistic relates to the pooling metric.
}
\references{
Gelman, Andrew, and Iain Pardoe.
"Bayesian Measures of Explained Variance and Pooling in Multilevel (Hierarchical) Models."
\emph{Technometrics 48, no. 2 (May 2006): 241-51}. \url{https://doi.org/10.1198/004017005000000517}.

Higgins, Julian P T, Simon G Thompson, Jonathan J Deeks, and Douglas G Altman.
"Measuring Inconsistency in Meta-Analyses."
\emph{British Medical Journal 327, no. 7414 (September 6, 2003): 557-60.}

Hippel, Paul T von. "The Heterogeneity Statistic I2 Can Be Biased in Small Meta-Analyses."
\emph{BMC Medical Research Methodology 15 (April 14, 2015).} \url{https://doi.org/10.1186/s12874-015-0024-z}.
}
