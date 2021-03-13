# Title: Simulation for Variable-wise kernel fuzzy c-means clustering algorithms with kernalization of the metric
# Author: Gabriel Teotonio
# Date: 2021/01/29

# Packages -----
library(tidyverse)
library(mclust)

# Data loading -----
data <- read_csv("data/data_banknote_authentication.txt", col_names = FALSE) %>% 
  rename("variance" = X1, "skewness" = X2, "curtosis" = X3, "entropy" = X4, "class" = X5) %>% 
  mutate(class = factor(class, levels = c(0, 1)))

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

data <- data %>% 
  mutate_at(c("variance", "skewness", "curtosis", "entropy"), ~normalize(.))

label_prior <- data %>% select("class") %>% as.vector()
x <- unname(as.matrix(data[,-5]))

# Monte Carlo -----
nrep <- 100
params_fuzzy <- c(1.01, 1.1, 1.6, 2)
result <- list()
best_J_objective <- 10^(10)
best <- list()
best_register <- list()
i <- 1

for (params in params_fuzzy) {
  for (n in 1:nrep) {
    
    result <- vkfcm_k_lp(x, m = params)
    
    measures <- quality_measures(result$Membership, label_prior)
    measures$int <- n
    measures$param <- params
    write.table(measures, "data/results/measures.csv", sep = ",",
                col.names = !file.exists("data/results/measures.csv"),
                append = TRUE, row.names = FALSE)
    
    if (result$J_objective < best_J_objective) {
      best_J_objective <- result$J_objective
      best <- result
    }
  }
  
  best_register[[i]] <- best
  i <- i + 1
  best_J_objective <- 10^(10)
  best <- list()
  
}

quality_measures(best_register[[2]]$Membership, label_prior)