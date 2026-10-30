generateRing <- function(N, R_inner, W, center=c(0, 0)) {
  R_outer <- R_inner + W
  theta <- runif(N, min=0, max=2*pi)
  
  u <- runif(N)
  r <- sqrt(u * (R_outer^2 - R_inner^2) + R_inner^2)
  
  x <- r * cos(theta) + center[1]
  y <- r * sin(theta) + center[2]
  
  cbind(x, y)
}

set.seed(20260930)
CENTER <- c(10, 0)
N <- 500
density_example <- rbind(
  ## Cross at 0, 0
  rmvnorm(N, sigma=matrix(c(1, -0.85, -0.85, 1), byrow=T, nrow=2)),
  rmvnorm(N, sigma=matrix(c(1, 0.85, 0.85, 1), byrow=T, nrow=2)),
  ## Ring w/ Cross in middle
  generateRing(2*N, 4, 0.75, center=CENTER),
  # rmvnorm(n=N, mean=CENTER, sigma=diag(c(0.5, 0.05))),
  # rmvnorm(n=N, mean=CENTER, sigma=diag(c(0.05, 0.5)))
  rmvnorm(n=2*N, mean=CENTER, sigma=diag(c(0.5, 0.5)))
)

usethis::use_data(density_example, overwrite = TRUE)