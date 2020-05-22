load("data-raw/microcredit_project_data.RData")

# Rachael Meager's code from the original project
# for consumer durables outcome


# This preps data so that the model can be written efficiently -
# it's not necessary but it's the best way I know to code it in STAN
group <- c( angelucci_indicator, attanasio_indicator, augsberg_indicator,
           banerjee_indicator, crepon_indicator)

consumption <- c(angelucci_consumption, attanasio_consumption,
                       augsberg_consumption,banerjee_consumption,
                       crepon_consumption)
treatment <- c( angelucci_treatment, attanasio_treatment, augsberg_treatment,
                banerjee_treatment, crepon_treatment)

# Now we have to standardise any variables which are in
# local currency units to USD PPP per fortnight in 2009 dollars
expanded_standardiser_USD_PPP_per_fortnight <- c(
  rep(the_consumption_standardiser_USD_PPP_per_fortnight[1],
      length(angelucci_indicator)),
  rep(the_consumption_standardiser_USD_PPP_per_fortnight[2],
      length(attanasio_indicator)),
  rep(the_consumption_standardiser_USD_PPP_per_fortnight[3],
      length(augsberg_indicator)),
  rep(the_consumption_standardiser_USD_PPP_per_fortnight[4],
      length(banerjee_indicator)),
  rep(the_consumption_standardiser_USD_PPP_per_fortnight[5],
      length(crepon_indicator)))

consumption <- consumption*expanded_standardiser_USD_PPP_per_fortnight

# bind everything into a data frame
data <- data.frame(group, consumption, treatment)

# We gotta remove the NA values for analyses
data <- data[complete.cases(data),]

# WW: now export for use in the package
microcredit_simplified <- data
# names(microcredit) <- c("group", "outcome", "treatment")
usethis::use_data(microcredit_simplified, overwrite = TRUE)

