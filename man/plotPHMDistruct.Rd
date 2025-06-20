% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis.R
\name{plotPHMDistruct}
\alias{plotPHMDistruct}
\title{Generate the distruct plot from the posterior matrix}
\usage{
plotPHMDistruct(
  phm,
  K = length(phm),
  colors = NULL,
  labels = NULL,
  axisTextSize = 6,
  partition = F
)
}
\arguments{
\item{phm}{Output from the \code{\link[=PHM]{PHM()}} function}

\item{K}{Number of clusters for which to generate the distruct plot}

\item{colors}{Optinal vector with colors for the mixture components}

\item{labels}{Ground truth class labels for the observations (ordered factor vector)}

\item{axisTextSize}{Size for axis labels}

\item{partition}{Whether to visualize from the posterior matrix or partition labels}
}
\value{
A ggplot object
}
\description{
Visualize a distruct plot based on either a partition or the posterior cluster probabilities.
}
\details{
In the case of visualizing for a partition, the posterior probabilities are set to 1 if it is the cluster the obesrvation is assigned to and 0 otherwise.
}
