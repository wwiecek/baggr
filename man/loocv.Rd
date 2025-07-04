% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loocv.R
\name{loocv}
\alias{loocv}
\title{Leave one group out cross-validation for \code{baggr} models}
\usage{
loocv(data, return_models = FALSE, ...)
}
\arguments{
\item{data}{Input data frame - same as for \link{baggr} function.}

\item{return_models}{logical; if FALSE, summary statistics will be returned and the
models discarded;
if TRUE, a list of models will be returned alongside summaries}

\item{...}{Additional arguments passed to \link{baggr}.}
}
\value{
log predictive density value, an object of class \code{baggr_cv};
full model, prior values and \emph{lpd} of each model are also returned.
These can be examined by using \code{attributes()} function.
}
\description{
Performs exact leave-one-group-out cross-validation on a baggr model.
}
\details{
The values returned by \code{loocv()} can be used to understand how excluding
any one group affects the overall result, as well as how well the model
predicts the omitted group. LOO-CV approaches are a good general practice
for comparing Bayesian models, not only in meta-analysis.

To learn about cross-validation see Gelman et al 2014.

This function automatically runs \emph{K} baggr models, where \emph{K} is number of groups (e.g. studies),
leaving out one group at a time. For each run, it calculates
\emph{expected log predictive density} (ELPD) for that group (see Gelman et al 2013).
(In the logistic model, where the proportion in control group is unknown, each of
the groups is divided into data for controls, which is kept for estimation, and data for
treated units, which is not used for estimation but only for calculating predictive density.
This is akin to fixing the baseline risk and only trying to infer the odds ratio.)

The main output is the cross-validation
information criterion, or -2 times the ELPD summed over \emph{K} models.
(We sum the terms as we are working with logarithms.)
This is related to, and often approximated by, the Watanabe-Akaike
Information Criterion. When comparing models, smaller values mean
a better fit.

For running more computation-intensive models, consider setting the
\code{mc.cores} option before running loocv, e.g. \code{options(mc.cores = 4)}
(by default baggr runs 4 MCMC chains in parallel).
As a default, rstan runs "silently" (\code{refresh=0}).
To see sampling progress, please set e.g. \code{loocv(data, refresh = 500)}.
}
\examples{
\dontrun{
# even simple examples may take a while
cv <- loocv(schools, pooling = "partial")
print(cv)      # returns the lpd value
attributes(cv) # more information is included in the object
}

}
\references{
Gelman, Andrew, Jessica Hwang, and Aki Vehtari.
'Understanding Predictive Information Criteria for Bayesian Models.'
Statistics and Computing 24, no. 6 (November 2014): 997–1016.
}
\seealso{
\link{loo_compare} for comparison of many LOO CV results; you can print and plot
output via \link{plot.baggr_cv} and \link{print.baggr_cv}
}
\author{
Witold Wiecek
}
