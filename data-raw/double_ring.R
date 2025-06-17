## code to prepare `double_ring` dataset goes here
N <- 5000
NPTS <- ceiling(N/2)
R1 <- 3
R2 <- 7
THETA <- runif(NPTS) * 2 * pi
R_SCALE <- runif(NPTS, 0, 1)

circle_inner <- cbind((R1 - R_SCALE) * cos(THETA),
                    (R1 - R_SCALE) * sin(THETA))
circle_outer <- cbind((R2 - R_SCALE) * cos(THETA),
                    (R2 - R_SCALE) * sin(THETA))

double_ring <- rbind(circle_inner, circle_outer)

usethis::use_data(double_ring, overwrite = TRUE)
