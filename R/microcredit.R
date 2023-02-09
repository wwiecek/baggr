#' 7 studies on effect of microcredit supply
#'
#' This dataframe contains the data used in Meager (2019) to estimate hierarchical
#' models on the data from 7 randomized controlled trials of expanding access to microcredit.
#'
#' The columns include the group indicator which gives the name of the lead author
#' on each of the respective studies, the value of the 6 outcome variables of most
#' interest (consumer durables spending, business expenditures, business profit,
#' business revenues, temptation goods spending and consumption spending) all of
#' which are standardised to USD PPP in 2009 dollars per two weeks (these are flow variables),
#' and finally a treatment assignment status indicator.
#'
#' The dataset has not otherwise been cleaned and therefore includes NAs and other
#' issues common to real-world datasets.
#'
#' For more information on how and why these variables were chosen and standardised,
#' see Meager (2019) or consult the associated code repository which includes the
#' standardisation scripts:
#' [link](https://bitbucket.org/rmeager/aggregate-average-impacts-microcredit/src/master/)
#'
#' @format A data frame with 40267 rows, 7 study identifiers and 7 outcomes
#' @references Meager, Rachael (2019) Understanding the average impact of microcredit expansions:
#' A Bayesian hierarchical analysis of seven randomized experiments.
#' American Economic Journal: Applied Economics, 11(1), 57-91.
"microcredit"



#' Simplified version of the microcredit dataset.
#'
#' This dataframe contains the data used in Meager (2019) to estimate hierarchical
#' models on the data from 7 randomized controlled trials of expanding access to microcredit.
#'
#' The columns include the group indicator which gives the name of the lead author on
#' each of the respective studies, the value of the household consumption
#' spending standardised to USD PPP in 2009 dollars per two weeks (these are flow variables),
#' and finally a treatment assignment status indicator.
#'
#' The dataset has not otherwise been cleaned and therefore includes NAs and other
#' issues common to real data.
#'
#' For more information on how and why these variables were chosen and standardised,
#' see Meager (2019) or consult the associated code repository:
#' [link](https://bitbucket.org/rmeager/aggregate-average-impacts-microcredit/src/master/)
#'
#' This dataset includes only complete cases and only the consumption outcome variable.
#'
#' @format A data frame with 14224 rows, 7 study identifiers and 1 outcome
#' @references Meager, Rachael (2019) Understanding the average impact of microcredit expansions:
#' A Bayesian hierarchical analysis of seven randomized experiments. American Economic Journal:
#' Applied Economics, 11(1), 57-91.
"microcredit_simplified"


#' 8 schools example
#'
#' A classic example of aggregate level continuous data in Bayesian hierarchical modelling.
#' This dataframe contains a column of estimated treatment effects of an SAT prep
#' program implemented in 8 different schools in the US, and a column of estimated standard errors.
#'
#' See Gelman et al (1995), Chapter 5, for context and applied example.
#'
#' @references Gelman, Andrew, John B. Carlin, Hal S. Stern, and Donald B. Rubin.
#' Bayesian Data Analysis. Taylor & Francis, 1995.
"schools"


#' Chickens: impact of electromagnetic field on calcium ion efflux in chicken brains
#'
#' An experiment conducted by Blackman et al. (1988) and documented in the following
#' [GitHub repository](https://github.com/VMatthijs/Slamming-the-sham) by Vakar
#' and Gelman. The dataset consists of a large number of experiments (`tau`, `se.tau`)
#' repeated at varying wave frequencies. Sham experiments (`mu`, `se.mu`) are also
#' included, allowing us to compare performance of models with and without
#' control measurements.
#'
#' @references Blackman, C. F., S. G. Benane, D. J. Elliott, D. E. House, and
#' M. M. Pollock.
#' “Influence of Electromagnetic Fields on the Efflux of Calcium Ions from Brain
#' Tissue in Vitro: A Three-Model Analysis Consistent with the Frequency
#' Response up to 510 Hz.” Bioelectromagnetics 9, no. 3 (1988): 215–27.
"chicks"

#' Yusuf et al: beta-blockers and heart attacks
#'
#' This replicates Table 6 from the famous Yusuf et al. (1985), removing one trial (Snow)
#' that had NA observations only. The paper is notable for application of rare-event
#' corrections, which we discuss more in package vignette `baggr_binary`.
#'
#' @references Yusuf, S., Peto, R., Lewis, J., Collins, R., & Sleight, P. (1985).
#' Beta blockade during and after myocardial infarction:
#' An overview of the randomized trials.
#' Progress in Cardiovascular Disease, 27(5), 335–371.
"yusuf"

#' Spike & slab example dataset
#'
"data_spike"
