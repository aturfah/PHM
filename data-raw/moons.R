## code to prepare `moons` dataset goes here
set.seed(2)
N <- 1000
R <- 4
theta <- runif(N) * 2 * pi
theta <- sort(theta)
r <- R + runif(N, max=1.5)

moons <- cbind(
  r * cos(theta), r * sin(theta)
)
moons[, 1] <- moons[, 1] + c(
  rep(R/3, N/2), rep(-R/3, N/2)
)

usethis::use_data(moons, overwrite = TRUE)
