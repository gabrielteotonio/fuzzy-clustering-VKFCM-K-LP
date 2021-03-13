# Title: Variable-wise kernel fuzzy c-means clustering algorithms with kernalization of the metric
# Author: Gabriel Teotonio
# Date: 2021/01/29


# Aux. functions -----
kernel_gaussian <- function(x, v) {
  
  dif <- (x - v)
  percentiles <- quantile(abs(dif)^2, names = FALSE, probs = c(.1, .9))
  two_sigma <- sum(percentiles)/2
  
  res <- c(length(x))
  for (k in 1:length(x)) {
    res[k] <- exp(-(dif[k]^2)/two_sigma)
  }
  
  return(res)
  
}

metric_membership <- function(u, m, x, v) {
  
  kernel_value <- kernel_gaussian(x, v)
  #print(kernel_value)
  res <- sum(u^m * (2 * (1 - kernel_value)))
  
  return(res)
  
}

metric_weight <- function(l, x, v) {
  
  kernel_value <- kernel_gaussian(x, v)
  res <- sum(l * (2 * (1 - kernel_value)))
  
  return(res)
  
}

quality_measures <- function(membership, label_prior) {
  
  label_posterior <- 
    t(membership) %>% 
    as.data.frame() %>% 
    rename("class_1" = V1, "class_2" = V2) %>% 
    transmute("class" = ifelse(class_1 > class_2, 0, 1)) %>% 
    mutate(class = factor(class, levels = c(0, 1)))
  
  # F-measure -----
  conf_mat <- confusionMatrix(label_posterior$class, label_prior$class)
  precision <- conf_mat$table[1,1] / (conf_mat$table[1,1] + conf_mat$table[1,2])
  recall <- conf_mat$table[1,1] / (conf_mat$table[1,1] + conf_mat$table[2,1])
  Fmeasure <- 2 * precision * recall / (precision + recall)
  
  # Adjusted Rand index -----
  Rand_index <- adjustedRandIndex(label_posterior$class, label_prior$class)
  
  # Modified partition coefficient -----
  Vpc <- (1/length(label_prior$class)) * sum(diag(membership %*% t(membership)))
  Mod_Vpc <- 1 - (2 * (1 - Vpc))
  
  # Partition entropy -----
  Vpe <- -(1/length(label_prior$class)) * sum(diag(membership %*% t(log(membership))))
  
  return(list("Fmeasure"= Fmeasure,
              "Rand_index" = Rand_index,
              "Mod_Vpc" = Mod_Vpc,
              "Vpe" = Vpe
  )
  )
  
}

# Algorithm -----
vkfcm_k_lp <- function(x, c = 2, m, T_limit = 150, error_cond = 10^(-10)) {
  
  # Initiate objects -----
  V <- matrix(0L, nrow = c, ncol = ncol(x))  
  U <- replicate(nrow(x), sample(c(0.5, 0.5)))
  L <- matrix(1, nrow = c, ncol = ncol(x))
  J <- 10^(10)
  J_aux <- 10^(10)
  t <- 1L
  aux_prod <- c(1, 1)
  aux_sum <- 0
  error <- 10^(10)
  
  initials_kick <- sample(1:nrow(x), c, replace = FALSE)
  V <- x[initials_kick, ]
  #print(V)
  while(error > error_cond && t <= T_limit) {
    
    # Update prototype -----
    for (i in 1:c) {
      for (j in 1:ncol(x)) {
        
        V[i, j] <- sum(U[i, ]^m * kernel_gaussian(x[, j], V[i, j]) * x[, j])/sum(U[i, ]^m * kernel_gaussian(x[, j], V[i, j]))
        
        aux_prod[i] <- aux_prod[i] * metric_membership(U[i, ], m, x[, j], V[i, j])
        
      }
    }
    
    # Update  weight ----
    for (i in 1:c) {
      for (j in 1:ncol(x)) {
        #print(metric_membership(U[i, ], m, x[, j], V[i, j]))
        L[i, j] <- ((aux_prod[i])^(1/ncol(x)))/metric_membership(U[i, ], m, x[, j], V[i, j])
        
      }
    }
    
    # Update membership ----
    for (i in 1:c) {
      for (k in 1:nrow(x)) {
        for (h in 1:c) {
          
          aux_sum <- aux_sum + (metric_weight(L[i, ], x[k, ], V[i, ])/metric_weight(L[i, ], x[k, ], V[h, ]))^(1/(m-1)) 
        } 
        
        U[i, k] <- aux_sum^(-1) 
        aux_sum <- 0
        
      }
    }
    
    J_aux <- J
    J <- 0
    
    # Update objective ----
    for (i in 1:c) {
      for (j in 1:ncol(x)) {
        
        J <- J + L[i, j] * metric_membership(U[i, ], m, x[, j], V[i, j])
        
      }
    }
    
    error <- abs(J_aux - J)
    t <- t + 1
    aux_prod <- c(1, 1)
    
  }
  
  return(list("initials"= initials_kick,
              "fuzziness_param" = m,
              "J_objective" = J,
              "Prototype" = V,
              "Weight" = L,
              "Membership" = U
  )
  )
}