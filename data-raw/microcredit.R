load("data-raw/microcredit_project_data.RData")

# Rachael Meager's code from the original project
# for consumer durables outcome


# This preps data so that the model can be written efficiently -
# it's not necessary but it's the best way I know to code it in STAN
site <- c( attanasio_indicator, augsberg_indicator,
           banerjee_indicator, crepon_indicator)
site <- site - 1
consumerdurables <- c( attanasio_consumerdurables,
                       augsberg_consumerdurables,banerjee_consumerdurables,
                       crepon_consumerdurables)
treatment <- c( attanasio_treatment, augsberg_treatment,
                banerjee_treatment, crepon_treatment)

# Now we have to standardise any variables which are in
# local currency units to USD PPP per fortnight in 2009 dollars
expanded_standardiser_USD_PPP_per_fortnight <- c(
  rep(the_consumerdurables_standardiser_USD_PPP_per_fortnight[1],
      length(attanasio_indicator)),
  rep(the_consumerdurables_standardiser_USD_PPP_per_fortnight[2],
      length(augsberg_indicator)),
  rep(the_consumerdurables_standardiser_USD_PPP_per_fortnight[3],
      length(banerjee_indicator)),
  rep(the_consumerdurables_standardiser_USD_PPP_per_fortnight[4],
      length(crepon_indicator)))

consumerdurables <- consumerdurables*expanded_standardiser_USD_PPP_per_fortnight

# bind everything into a data frame
data <- data.frame(site, consumerdurables, treatment)

# We gotta remove the NA values for analyses
data <- data[complete.cases(data),]

# WW: now export for use in the package
microcredit_simplified <- data
# names(microcredit) <- c("site", "outcome", "treatment")
devtools::use_data(microcredit_simplified)

