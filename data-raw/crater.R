## code to prepare `crater` dataset goes here

set.seed(1)

core_radius <- 4
ring_radius <- 5
ring_width <- 1.5

## Inner core points
N <- 1000
theta <- runif(N) * 2 * pi
r <- core_radius * runif(N)
inner_core <- cbind(
  r * cos(theta),
  r * sin(theta)
)


## Outer ring points
set.seed(2)
N <- 750
theta <- runif(N) * 2 * pi
r <- ring_radius + runif(N) * ring_width
outer_ring <- cbind(
  r * cos(theta),
  r * sin(theta)
)

## Combine
crater <- rbind(
  inner_core,
  outer_ring
)

usethis::use_data(crater, overwrite = TRUE)
