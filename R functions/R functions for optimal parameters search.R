############################################################
# Helper function 1:
# Stage-I power under the least favorable configuration
############################################################

stage1_power <- function(n1,
                         gamma,
                         delta1,
                         delta2,
                         k,
                         sigma) {
  
  lower <- -sqrt(n1) * (delta2 - gamma) / sigma
  
  integrand <- function(z) {
    pnorm(
      z + sqrt(n1) * (delta2 - delta1) / sigma
    )^(k - 1) * dnorm(z)
  }
  
  integrate(
    integrand,
    lower = lower,
    upper = Inf
  )$value
}


############################################################
# Helper function 2:
# Stage-II sample size
############################################################

stage2_sample_size <- function(power2,
                               alpha,
                               delta2,
                               sigma) {
  
  ceiling(
    sigma^2 *
      (qnorm(1 - alpha) + qnorm(power2))^2 /
      delta2^2
  )
}


############################################################
# Helper function 3:
# Expected sample sizes, cost, and maximum sample size
############################################################

sampling_quantities <- function(n1,
                                n2,
                                gamma,
                                delta1,
                                delta2,
                                k,
                                sigma) {
  
  # Probability of terminating after Stage I
  # under the global null
  p_terminate_null <-
    pnorm(
      gamma / (sigma / sqrt(n1))
    )^k
  
  
  # Expected sample size under the global null
  EN_null <-
    k * n1 +
    n2 * (1 - p_terminate_null)
  
  
  # Probability of terminating after Stage I
  # under the least favorable configuration
  p_terminate_LFC <-
    pnorm(
      (gamma - delta1) / (sigma / sqrt(n1))
    )^(k - 1) *
    pnorm(
      (gamma - delta2) / (sigma / sqrt(n1))
    )
  
  
  # Expected sample size under the LFC
  EN_LFC <-
    k * n1 +
    n2 * (1 - p_terminate_LFC)
  
  
  # Sampling-cost criterion
  cost <-
    0.5 * EN_null +
    0.5 * EN_LFC
  
  
  # Maximum possible sample size
  Nmax <- k * n1 + n2
  
  
  c(
    p_terminate_null = p_terminate_null,
    EN_null = EN_null,
    EN_LFC = EN_LFC,
    cost = cost,
    Nmax = Nmax
  )
}


############################################################
# Main function 1:
# Only overall power is specified
############################################################

optimize_design1 <- function(n1_range,
                             gamma_range,
                             overall_power,
                             alpha,
                             delta1,
                             delta2,
                             k,
                             sigma) {
  
  results <- list()
  count <- 0
  
  
  for (n1 in n1_range) {
    
    for (gamma in gamma_range) {
      
      # Compute Stage-I power
      power1 <- stage1_power(
        n1 = n1,
        gamma = gamma,
        delta1 = delta1,
        delta2 = delta2,
        k = k,
        sigma = sigma
      )
      
      
      # Stage-I power must be strictly larger
      # than the target overall power
      if (power1 > overall_power) {
        
        # Required Stage-II power
        power2 <- overall_power / power1
        
        
        # Compute Stage-II sample size
        n2 <- stage2_sample_size(
          power2 = power2,
          alpha = alpha,
          delta2 = delta2,
          sigma = sigma
        )
        
        
        # Compute sampling quantities
        q <- sampling_quantities(
          n1 = n1,
          n2 = n2,
          gamma = gamma,
          delta1 = delta1,
          delta2 = delta2,
          k = k,
          sigma = sigma
        )
        
        
        count <- count + 1
        
        results[[count]] <- c(
          n1 = n1,
          n2 = n2,
          gamma = gamma,
          p_terminate_null = q[["p_terminate_null"]],
          EN_null = q[["EN_null"]],
          EN_LFC = q[["EN_LFC"]],
          cost = q[["cost"]],
          Nmax = q[["Nmax"]],
          power1 = power1,
          power2 = power2
        )
      }
    }
  }
  
  
  if (length(results) == 0) {
    stop("No feasible design was found.")
  }
  
  
  results <- do.call(rbind, results)
  
  
  # Select the combination minimizing the sampling cost
  best <- which.min(results[, "cost"])
  
  result <- results[best, ]
  
  
  names(result) <- c(
    "n1",
    "n2",
    "gamma",
    "P_terminate_null",
    "EN_null",
    "EN_LFC",
    "cost",
    "Nmax",
    "power_stage1",
    "power_stage2"
  )
  
  
  return(result)
}


############################################################
# Main function 2:
# Overall power and Stage-II power are specified
############################################################

optimize_design2 <- function(n1_range,
                             gamma_range,
                             overall_power,
                             stage2_power,
                             alpha,
                             delta1,
                             delta2,
                             k,
                             sigma) {
  
  # Stage-II power must be strictly larger
  # than overall power
  if (stage2_power <= overall_power) {
    stop(
      "Stage-II power must be strictly larger than overall power."
    )
  }
  
  
  # Compute Stage-II sample size first
  n2 <- stage2_sample_size(
    power2 = stage2_power,
    alpha = alpha,
    delta2 = delta2,
    sigma = sigma
  )
  
  
  # Required minimum Stage-I power
  required_power1 <-
    overall_power / stage2_power
  
  
  results <- list()
  count <- 0
  
  
  for (n1 in n1_range) {
    
    for (gamma in gamma_range) {
      
      # Compute Stage-I power
      power1 <- stage1_power(
        n1 = n1,
        gamma = gamma,
        delta1 = delta1,
        delta2 = delta2,
        k = k,
        sigma = sigma
      )
      
      
      # Retain feasible designs
      if (power1 >= required_power1) {
        
        # Compute sampling quantities
        q <- sampling_quantities(
          n1 = n1,
          n2 = n2,
          gamma = gamma,
          delta1 = delta1,
          delta2 = delta2,
          k = k,
          sigma = sigma
        )
        
        
        count <- count + 1
        
        results[[count]] <- c(
          n1 = n1,
          n2 = n2,
          gamma = gamma,
          p_terminate_null = q[["p_terminate_null"]],
          EN_null = q[["EN_null"]],
          EN_LFC = q[["EN_LFC"]],
          cost = q[["cost"]],
          Nmax = q[["Nmax"]],
          power1 = power1,
          power2 = stage2_power
        )
      }
    }
  }
  
  
  if (length(results) == 0) {
    stop("No feasible design was found.")
  }
  
  
  results <- do.call(rbind, results)
  
  
  # Select the combination minimizing the sampling cost
  best <- which.min(results[, "cost"])
  
  result <- results[best, ]
  
  
  names(result) <- c(
    "n1",
    "n2",
    "gamma",
    "P_terminate_null",
    "EN_null",
    "EN_LFC",
    "cost",
    "Nmax",
    "power_stage1",
    "power_stage2"
  )
  
  
  return(result)
}



 