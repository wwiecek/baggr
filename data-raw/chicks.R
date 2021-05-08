# data from https://github.com/VMatthijs/Slamming-the-sham 12 Mar 2021

# Reading in chicken data
chicks <- read.table("data-raw/chickens.dat", header = TRUE)
names(chicks) <- c("frequency", "mu_n", "mu", "se.mu", "tau_n", "tau", "se.tau")
usethis::use_data(chicks, overwrite=TRUE)
