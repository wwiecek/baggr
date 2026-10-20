# data from https://github.com/VMatthijs/Slamming-the-sham 12 Mar 2021

# Reading in chicken data
chicks <- read.table("data-raw/chickens.dat", header = TRUE)
names(chicks) <- c("frequency", "n0i", "y0i", "se0i", "ni", "yi", "sei")
usethis::use_data(chicks, overwrite=TRUE)
