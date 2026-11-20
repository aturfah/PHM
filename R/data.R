#' Simulated two-dimensional datasets
#'
#' @description
#' Small simulated datasets used in the PHM examples and user guide. Each is a
#' plain numeric matrix with one row per point and two columns (x, y), without
#' cluster labels. The scripts that generated them are in `data-raw/` of the
#' source repository.
#'
#' - `crater`: 1750 points. 1000 points in a disk of radius 4 (the radius of
#'   each point is uniform on 0 to 4, so the center is denser), surrounded by
#'   750 points in a ring with radius between 5 and 6.5.
#' - `density_example`: 3000 points. Two crossing elongated Gaussians at the
#'   origin (500 points each, correlations -0.85 and 0.85), and, centered at
#'   (10, 0), a ring with inner radius 4 and width 0.75 (1000 points spread
#'   uniformly over its area) around a round Gaussian (1000 points, variance
#'   0.5 in each direction).
#' - `double_ring`: 5000 points on two concentric rings centered at the
#'   origin, with radius between 2 and 3 (rows 1 to 2500) and between 6 and 7
#'   (rows 2501 to 5000). The two rings were generated from the same random
#'   angles and radial offsets.
#' - `moons`: 1000 points on a ring with radius between 4 and 5.5, cut into
#'   two halves by angle, with the halves shifted 4/3 to the right (rows 1 to
#'   500) and to the left (rows 501 to 1000).
#' - `nested_diamonds`: 6400 points in 64 round Gaussian clusters of 100
#'   points each (standard deviation 0.5), arranged in three nested levels:
#'   four groups at (±8, ±8), each made of four diamonds offset by 4, each
#'   made of four clusters offset by 1.5.
#'
#' @format Numeric matrices with two columns: `crater` 1750 x 2,
#'   `density_example` 3000 x 2, `double_ring` 5000 x 2, `moons` 1000 x 2,
#'   `nested_diamonds` 6400 x 2.
#'
#' @examples
#' dim(double_ring)
#' plot(double_ring, pch = ".", asp = 1)
#'
#' @name phm_datasets
#' @aliases crater density_example double_ring moons nested_diamonds
#' @docType data
#' @keywords datasets
NULL
