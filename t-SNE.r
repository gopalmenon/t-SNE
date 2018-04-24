library("MASS")

INITIAL_SOLUTION_VARIANCE = 0.0001
INITIAL_DENSITY = 1
DENSITY_PRECISION = 0.1

## Get normally distributed initial solution
## number_of_points: size of data
## dimensionality: dimensionality of each data point
## variance: variance of initial solution used to compute covariance as variance x I
##
## Will return the normally distributed initial solution
get_initial_solution <- function(number_of_points, dimensionality, variance=INITIAL_SOLUTION_VARIANCE) {
  
  return(mvrnorm(n=number_of_points, mu=rep(0, dimensionality), Sigma=variance * diag(dimensionality)))

}

## Get similarity of point j to point i
## point_j_index: index of point j in the data array
## point_i_index: index of point i in the data array
## high_dimensional_data: data to be visualized
## density_at_i: density of point at i
##
## Will return the un-normalized similarity of point j to point i
get_un_normalized_similarity_j_to_i <- function(point_j_index, point_i_index, high_dimensional_data, density_at_i) {
  
  return(exp(-1 * sum((high_dimensional_data[point_i_index, ] - high_dimensional_data[point_j_index, ])^2)/(2 * density_at_i ^ 2)))
  
}

## Get entropy at point i for given density
## point_i_index: index of point i in the data array
## density: density at point i 
## high_dimensional_data: data to be visualized
##
## Will return the entropy at point i for given density
get_entropy_at_i <- function(point_i_index, density, high_dimensional_data) {
  
  un_normalized_similarity_at_i <- lapply(seq(1: length(high_dimensional_data[, 1]))[-point_i_index], 
                                          get_un_normalized_similarity_j_to_i, 
                                          point_i_index=point_i_index, 
                                          high_dimensional_data=high_dimensional_data, 
                                          density_at_i=density)
  
  normalized_similarity_at_i <- un_normalized_similarity_at_i / sum(exp(-1 * colSums((high_dimensional_data[point_i_index, ] - high_dimensional_data[-point_i_index, ])^2)/(2 * density ^ 2)))

  entropy_at_i <- sum(normalized_similarity_at_i * log2(normalized_similarity_at_i))

}

## Get density range corresponding to perplexity at point i 
## perplexity: required perplexity 
## point_i_index: index of point i in the data array
## high_dimensional_data: data to be visualized
##
## Will return the density range corresponding to perplexity
get_density_range <- function(perplexity, point_i_index, high_dimensional_data) {
  
  # Find the search direction
  log_perplexity <- log2(perplexity)
  entropy_at_i <- get_entropy_at_i(point_i_index, density=INITIAL_DENSITY, high_dimensional_data)
  if (entropy_at_i > log_perplexity) {
    search_ascending <- FALSE
    density_range_end <- INITIAL_DENSITY
  } else {
    search_ascending <- TRUE
    density_range_start <- INITIAL_DENSITY
  }
  
  # Search for the other end of the density range
  repeat {
    
    # Search for range end
    if (isTRUE(search_ascending)) {
      density_range_end <- density_range_start * 2
      entropy_at_i <- get_entropy_at_i(point_i_index, density=density_range_end, high_dimensional_data)
      if (entropy_at_i > log_perplexity) {
        break
      } else {
        density_range_start <- density_range_end
      }
    } else {
      density_range_start <- density_range_end / 2
      entropy_at_i <- get_entropy_at_i(point_i_index, density=density_range_start, high_dimensional_data)
      if (entropy_at_i < log_perplexity) {
        break
      } else {
        density_range_end <- density_range_start
      }
    }
    
  }
  
  density_range <- c(density_range_start, density_range_end)
  names(density_range) <- c("Density Range Start", "Density Range End")
  return(density_range)
  
}


## Get density corresponding to perplexity at point i 
## perplexity: required perplexity 
## point_i_index: index of point i in the data array
## high_dimensional_data: data to be visualized
##
## Will return the density corresponding to perplexity
get_density_at_i <- function(perplexity, point_i_index, high_dimensional_data) {
  
  # Get the density range
  density_range <- get_density_range(perplexity, point_i_index, high_dimensional_data)
  density_range_low <- density_range[1]
  density_range_high <- density_range[2]
  next_test_density <- 0.0
  
  # Do binary search to get density corresponding to perplexity
  repeat {
    
    next_test_density <- (density_range_low + density_range_high) / 2
    
    # Get entropy at test density
    entropy_at_i <- get_entropy_at_i(point_i_index, density=next_test_density, high_dimensional_data)
    if (abs(entropy_at_i - log_perplexity) > DENSITY_PRECISION) {
      if (entropy_at_i > log_perplexity) {
        density_range_high <- next_test_density
      } else {
        density_range_low <- next_test_density
      }
    } else {
      break
    }
    
  }
  
  return(next_test_density)
  
}