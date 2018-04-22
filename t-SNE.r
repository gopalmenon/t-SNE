library("MASS")

INITIAL_SOLUTION_VARIANCE = 0.0001

## Get normally distributed initial solution
## number_of_points: size of data
## dimensionality: dimensionality of each data point
## variance: variqnce of initial solution
##
## Will return the normally distributed initial solution
get_initial_solution <- function(number_of_points, dimensionality, variance=INITIAL_SOLUTION_VARIANCE) {
  
  # Create covariance matrix
  covariance_matrix <- variance * diag(dimensionality)
  
  return(mvrnorm(n=number_of_points, mu=rep(0, dimensionality), Sigma=variance * diag(dimensionality)))

}
  