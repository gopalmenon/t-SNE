library("MASS")

INITIAL_SOLUTION_VARIANCE = 0.0001
INITIAL_DENSITY = 3162
DENSITY_PRECISION = 0.1
DEFAULT_PERPLEXITY = 5
TARGET_LOW_DIMENSIONALITY = 3
NUMBER_OF_ITERATIONS = 1000
INITIAL_LEARNING_RATE = 100
LEARNING_RATE_NORMALIZER = 200
DECAY_LEARNING_RATE = TRUE
INITIAL_MOMENTUM = 0.5
FINAL_MOMENTUM = 0.8

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
  
  if (point_j_index == point_i_index) {
    return(0.0)
  } else {
    return(exp(-1 * sum((high_dimensional_data[point_i_index, ] - high_dimensional_data[point_j_index, ])^2)/(2 * density_at_i ^ 2)))
  }
  
}

## Get similarity of point j to point i in low dimensions
## point_j_index: index of point j in the data array
## point_i_index: index of point i in the data array
## low_dimensional_data: data that will be visualized
##
## Will return the un-normalized similarity of point j to point i
get_un_normalized_similarity_j_to_i_low_dimension <- function(point_j_index, point_i_index, low_dimensional_data, do_computation) {
  
  if (point_j_index == point_i_index || !isTRUE(do_computation)) {
    return(0.0)
  } else {
    return(1 / (1 + sum((low_dimensional_data[point_i_index, ] - low_dimensional_data[point_j_index, ])^2)))
  }
  
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
  
  normalized_similarity_at_i <- unlist(un_normalized_similarity_at_i) / sum(unlist(un_normalized_similarity_at_i))

  entropy_at_i <- sum(-normalized_similarity_at_i * ifelse(normalized_similarity_at_i != 0, log2(normalized_similarity_at_i), 0))

  return(entropy_at_i)
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
  log_perplexity <- log2(perplexity)
  
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

## Get density corresponding to perplexity at each point
## perplexity: required perplexity 
## high_dimensional_data: data to be visualized
##
## Will return a density vector corresponding to perplexity
get_density_at_each_point <- function(perplexity, high_dimensional_data) {
  
  return(unlist(lapply(seq(1: length(high_dimensional_data[, 1])), 
                       get_density_at_i, 
                       perplexity=perplexity,
                       high_dimensional_data=high_dimensional_data)))
  
}

## Get high dimensional pairwise affinities with perplexity
## perplexity: required perplexity 
## high_dimensional_data: data to be visualized
##
## Will return high dimensional pairwise affinities
get_high_dimensional_pairwise_affinities <- function(perplexity, high_dimensional_data) {
  
  point_i_index <- matrix(data=rep(seq(1:length(high_dimensional_data[,1])),length(high_dimensional_data[,1])),
                          nrow=length(high_dimensional_data[,1]),
                          ncol=length(high_dimensional_data[,1]),
                          byrow=FALSE)
  
  point_j_index <- matrix(data=rep(seq(1:length(high_dimensional_data[,1])),length(high_dimensional_data[,1])),
                          nrow=length(high_dimensional_data[,1]),
                          ncol=length(high_dimensional_data[,1]),
                          byrow=TRUE)
  
  density_at_i <- get_density_at_each_point(perplexity, high_dimensional_data)
  density_at_i_matrix <- matrix(rep(density_at_i, length(high_dimensional_data[,1])), 
                                nrow=length(high_dimensional_data[,1]), 
                                ncol=length(high_dimensional_data[,1]),
                                byrow=FALSE)
  
  un_normalized_pairwise_affinities <- matrix(mapply(get_un_normalized_similarity_j_to_i, 
                                                     point_j_index=c(t(point_j_index)), 
                                                     point_i_index=c(t(point_i_index)),
                                                     high_dimensional_data=list(high_dimensional_data),
                                                     c(t(density_at_i_matrix))),
                                              nrow=length(high_dimensional_data[,1]),
                                              ncol=length(high_dimensional_data[,1]))
  
  normalized_pairwise_affinities <- un_normalized_pairwise_affinities/rowSums(un_normalized_pairwise_affinities)
  return(get_symmetrized_high_dimensional_pairwise_affinities(high_dimensional_pairwise_affinities=normalized_pairwise_affinities,
                                                              number_of_points=length(high_dimensional_data[,1])))
  
}

## Get symmetrized affinity for a pair of points
## affinity_i_j: high dimensional affinity between point i and j
## affinity_j_i: high dimensional affinity between point j and i
##
## Will return symmetrized symmetrized affinity for a pair of points defined as Pij = (Pj|i + Pi|j)/2n
get_symmetrized_affinity <- function(affinity_i_j, affinity_j_i, number_of_points) {
  
  if (is.na(affinity_i_j) || is.na(affinity_j_i)) {
    return(0.0)
  } else {
    return((affinity_i_j + affinity_j_i)/(2 * number_of_points))
  }
  
}

## Get symmetrized high dimensional pairwise affinities
## high_dimensional_pairwise_affinities: high dimensional pairwise affinities
## number_of_points: number of data points
##
## Will return symmetrized high dimensional pairwise affinities defined as Pij = (Pj|i + Pi|j)/2n
get_symmetrized_high_dimensional_pairwise_affinities <- function(high_dimensional_pairwise_affinities, number_of_points) {
  
  # Get transpose of pairwise affinities and make it upper triangular so that symmetrized computation
  # can be done for corresponding elements.
  corr_high_dimensional_pairwise_affinities = t(high_dimensional_pairwise_affinities)
  
  
  lower_triangular <- matrix(rep(1, (number_of_points * number_of_points)), 
                             nrow=number_of_points, 
                             ncol=number_of_points)
  
  lower_triangular <- lower.tri(lower_triangular, diag = FALSE)
  
  high_dimensional_pairwise_affinities[lower_triangular] <- NA
  corr_high_dimensional_pairwise_affinities[lower_triangular] <- NA
  
  return(matrix(mapply(get_symmetrized_affinity, 
                       affinity_i_j=c(t(high_dimensional_pairwise_affinities)), 
                       affinity_j_i=c(t(corr_high_dimensional_pairwise_affinities)), 
                       number_of_points=number_of_points),
                nrow=number_of_points,
                ncol=number_of_points,
                byrow=TRUE))
  
}

## Get low dimensional pairwise affinities
## low_dimensional_data: data that will be visualized
##
## Will return low dimensional pairwise affinities
get_low_dimensional_pairwise_affinities <- function(low_dimensional_data) {
  
  point_i_index <- matrix(data=rep(seq(1:length(low_dimensional_data[,1])),length(low_dimensional_data[,1])),
                          nrow=length(low_dimensional_data[,1]),
                          ncol=length(low_dimensional_data[,1]),
                          byrow=FALSE)
  
  point_j_index <- matrix(data=rep(seq(1:length(low_dimensional_data[,1])),length(low_dimensional_data[,1])),
                          nrow=length(low_dimensional_data[,1]),
                          ncol=length(low_dimensional_data[,1]),
                          byrow=TRUE)
  
  upper_triangular <- matrix(rep(1, (length(low_dimensional_data[,1]) * length(low_dimensional_data[,1]))), 
                             nrow=length(low_dimensional_data[,1]), 
                             ncol=length(low_dimensional_data[,1]))
  
  upper_triangular <- upper.tri(upper_triangular, diag = FALSE)
    
  un_normalized_pairwise_affinities <- matrix(mapply(get_un_normalized_similarity_j_to_i_low_dimension, 
                                                     point_j_index=c(t(point_j_index)), 
                                                     point_i_index=c(t(point_i_index)),
                                                     low_dimensional_data=list(low_dimensional_data),
                                                     do_computation=c(t(upper_triangular))),
                                              nrow=length(low_dimensional_data[,1]),
                                              ncol=length(low_dimensional_data[,1]),
                                              byrow=TRUE)
  
  return(un_normalized_pairwise_affinities/sum(un_normalized_pairwise_affinities))
  
}

## Get gradient of KL Divergence between high and low dimensional pairwise affinities at a particular point
## point_i_index: point at which gradient is to be computed
## high_dimensional_pairwise_affinities: upper triangular matrix of high dimensional symmetrized pairwise affinities
## low_dimensional_pairwise_affinities: upper triangular matrix of low dimensional symmetrized pairwise affinities
## low_dimensional_points: current estimate of low dimensional points
##
## Will return the gradient of KL Divergence between high and low dimensional pairwise affinities at a particular point
get_gradient_at_point_i <- function(point_i_index, high_dimensional_pairwise_affinities, low_dimensional_pairwise_affinities, low_dimensional_points) {
  
  number_of_points <- length(low_dimensional_pairwise_affinities[1, ])
  
  pairwise_affinity_points <- matrix(rep(NA, (number_of_points * number_of_points)), 
                             nrow=number_of_points, 
                             ncol=number_of_points)
  
  # Pairwise affinity point to be used for gradient at point i
  if (point_i_index < number_of_points) {
    pairwise_affinity_points[point_i_index, seq(point_i_index + 1, number_of_points)] <- 1
    affinities_used_2 <- seq(point_i_index + 1, number_of_points)
  } else {
    affinities_used_2 <- numeric()
  }
  
  if (point_i_index > 1) {
    pairwise_affinity_points[seq(1, point_i_index - 1), point_i_index] <- 1
    affinities_used_1 <- seq(1, point_i_index - 1)
   } else {
    affinities_used_1 <- numeric()
   }
  
  affinities_used <- c(affinities_used_1, affinities_used_2)
  
  # Compute gradient
  Pij_Qij <- high_dimensional_pairwise_affinities[!is.na(pairwise_affinity_points)] - 
    low_dimensional_pairwise_affinities[!is.na(pairwise_affinity_points)]
  
  Yi_Yj <- low_dimensional_points[rep(point_i_index, number_of_points - 1),] - 
    low_dimensional_points[affinities_used, ]
  
  gradient_at_point_i <- Pij_Qij * Yi_Yj / (1 + rowSums(Yi_Yj^2))
  
  return(4 * colSums(gradient_at_point_i))
  
}

## Get the learning rate for the current iteration
## initial_learning_rate: learning rate for first iteration
## iteration_number: iteration number
## learning_rate_normalizer: how fast to decay the learning rate
## decay_learning_rate: should the frate decay or stay constant
##
## Will return the learning rate for the current iteration
get_learning_rate <- function(initial_learning_rate, iteration_number, learning_rate_normalizer, decay_learning_rate) {
  
  if (isTRUE(decay_learning_rate)) {
    return(initial_learning_rate/(1 + initial_learning_rate * (iteration_number - 1) / learning_rate_normalizer))
  } else {
    return(initial_learning_rate)
  }
  
}

## Get the next estimate for the low dimensional point by doing a gradient descent
## previous_low_dimensional_point_estimate: previous estimate for the low dimensional point
## current_learning_rate: learning rate at the current step
## kl_divergence_gradient_at_point: gradient at the point
## momentum_term: momentum for convergence
## before_previous_low_dimensional_point_estimate: estimate for the low dimensional point two steps before
##
## Will return the next estimate for the low dimensional point
get_next_low_dimensional_point_estimate <- function(previous_low_dimensional_point_estimates, previous_low_dimensional_point_estimate_index, current_learning_rate, kl_divergence_gradient_at_point, momentum_term, before_previous_low_dimensional_point_estimates) {
  
  return(previous_low_dimensional_point_estimates[previous_low_dimensional_point_estimate_index, ] +
           current_learning_rate * kl_divergence_gradient_at_point[previous_low_dimensional_point_estimate_index, ] +
           momentum_term * (ifelse(all(before_previous_low_dimensional_point_estimates[previous_low_dimensional_point_estimate_index, ] == 0),
                                   rep(0, length(previous_low_dimensional_point_estimates[previous_low_dimensional_point_estimate_index, ])),
                                   (previous_low_dimensional_point_estimates[previous_low_dimensional_point_estimate_index, ] - 
                                      before_previous_low_dimensional_point_estimates[previous_low_dimensional_point_estimate_index, ]))))
  
}

run_t_sne <- function(high_dimensional_data, low_dimensionality=TARGET_LOW_DIMENSIONALITY, perplexity=DEFAULT_PERPLEXITY, number_of_iterations=NUMBER_OF_ITERATIONS, initial_learning_rate=INITIAL_LEARNING_RATE, learning_rate_normalizer=LEARNING_RATE_NORMALIZER, initial_momentum=INITIAL_MOMENTUM, final_momentum=FINAL_MOMENTUM) {
  
  # Compute symmetrized pairwise affinities with perplexity
  high_dimensional_pairwise_affinities <- get_high_dimensional_pairwise_affinities(perplexity=perplexity, high_dimensional_data=high_dimensional_data)
  
  # Sample initial solution from a normal distribution
  previous_low_dimensional_point_estimate <- get_initial_solution(number_of_points=length(high_dimensional_data[, 1]), dimensionality=low_dimensionality)
  before_previous_low_dimensional_point_estimate <- matrix(rep(0, length(high_dimensional_data[, 1]) * low_dimensionality),
                                                           nrow=length(high_dimensional_data[, 1]),
                                                           ncol=low_dimensionality,
                                                           byrow=TRUE)
  
  # Find low dimensional representation using gradient descent
  for (iteration_counter in 1:number_of_iterations) {
  
    # Compute low dimensional affinities
    low_dimensional_pairwise_affinities <- get_low_dimensional_pairwise_affinities(low_dimensional_data=previous_low_dimensional_point_estimate)
  
    # Compute gradient at each point in low dimensions
    gradient_at_each_low_dimensional_point <- matrix(mapply(get_gradient_at_point_i,
                                                            point_i_index=rep(1:length(high_dimensional_data[, 1])), 
                                                            high_dimensional_pairwise_affinities=list(high_dimensional_pairwise_affinities),
                                                            low_dimensional_pairwise_affinities=list(low_dimensional_pairwise_affinities),
                                                            low_dimensional_points=list(previous_low_dimensional_point_estimate)),
                                                     nrow=length(high_dimensional_data[, 1]),
                                                     ncol=low_dimensionality,
                                                     byrow=TRUE)
    

    # Take a step by doing gradient descent
    learning_rate <- get_learning_rate(initial_learning_rate=initial_learning_rate, iteration_number=iteration_counter, learning_rate_normalizer=learning_rate_normalizer, decay_learning_rate=DECAY_LEARNING_RATE)
    next_low_dimensional_point_estimate <- matrix(mapply(get_next_low_dimensional_point_estimate,
                                                         previous_low_dimensional_point_estimates=list(previous_low_dimensional_point_estimate),
                                                         previous_low_dimensional_point_estimate_index=seq(1:length(high_dimensional_data[, 1])),
                                                         current_learning_rate=learning_rate,
                                                         kl_divergence_gradient_at_point=list(gradient_at_each_low_dimensional_point),
                                                         momentum_term=ifelse(iteration_counter < 250, initial_momentum, final_momentum),
                                                         before_previous_low_dimensional_point_estimates=list(before_previous_low_dimensional_point_estimate)),
                                                  nrow=length(high_dimensional_data[, 1]),
                                                  ncol=low_dimensionality,
                                                  byrow=TRUE)

    # Save previous values
    before_previous_low_dimensional_point_estimate <- previous_low_dimensional_point_estimate
    previous_low_dimensional_point_estimate <- next_low_dimensional_point_estimate

  }
  
  return(next_low_dimensional_point_estimate)
  
}
