## code to prepare `nested_diamonds` dataset goes here
set.seed(2)
Nk <- 100
TOTAL <- 64
nested_diamonds <- Reduce(
  rbind,
  lapply((1:TOTAL)-1, function(idx) {
    OFFSET1 <- 8
    OFFSET2 <- 4
    OFFSET3 <- 1.5

    ## Overall Position
    flipX <- 1
    flipY <- 1
    if (idx < TOTAL * 1/4) {
      flipX <- -1
    } else if (idx < TOTAL * 1/2) {
      flipY <- -1
    } else if (idx < TOTAL * 3/4) {
      flipX <- -1
      flipY <- -1
    }
    
    ## Position within cluster
    scale <- 1
    if (idx %% 4 == 0) {
      offset <- c(0, OFFSET2)
    } else if (idx %% 4 == 1) {
      offset <- c(0, -OFFSET2)
    } else if (idx %% 4 == 2) {
      offset <- c(OFFSET2, 0)
    } else {
      offset <- c(-OFFSET2, 0)
    }
    offsetX1 <- offset[1]
    offsetY1 <- offset[2]
    
    ## Position even more nested
    if (floor(idx/4) %% 4 == 0) {
      offsetX2 <- 0
      offsetY2 <- OFFSET3
    } else if (floor(idx/4) %% 4 == 1) {
      offsetX2 <- OFFSET3
      offsetY2 <- 0
    } else if (floor(idx/4) %% 4 == 2) {
      offsetX2 <- 0
      offsetY2 <- -OFFSET3
    } else {
      offsetX2 <- -OFFSET3
      offsetY2 <- 0
    } 

    cbind(
      rnorm(scale * Nk, OFFSET1*flipX + offsetX1 + offsetX2, 0.5),
      rnorm(scale * Nk, OFFSET1*flipY + offsetY1 + offsetY2, 0.5)
    )
  })
)

usethis::use_data(nested_diamonds, overwrite = TRUE)
