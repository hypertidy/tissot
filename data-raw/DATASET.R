## code to prepare `DATASET` dataset goes here

w <- do.call(cbind, maps::map(plot = FALSE)[1:2])
w[!w[,1] > -180, 1] <- -179.9
w[!w[,1] < 180, 1] <- 179.9
world <- w
usethis::use_data(world, overwrite = FALSE, compress = "xz")
