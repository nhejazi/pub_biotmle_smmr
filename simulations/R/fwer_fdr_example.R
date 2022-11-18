# using numerous hypothesis tests
many_tests <- function(n, k, mu, alpha = 0.05, correct_mc = FALSE) {
  test_stats <- rnorm(n, mean = 0)
  test_stats[seq_len(k)] <- test_stats[seq_len(k)] + mu

  if (correct_mc) {
    rejections <- sum(test_stats > qnorm(1 - alpha / n))
  } else {
    rejections <- sum(test_stats > qnorm(1 - alpha))
  }

  # Family-wise Error (FWE)
  fwe <- as.integer(any(rejections > k))

  # False Discovery Proportion (FDP)
  fdp <- sum(rejections > k) / max(1, length(rejections))

  # output
  return(c(fwe, fdp))
}

# test example
set.seed(74591)
B <- 10000
n <- 150
k <- 0   # global null
mu <- 1
alpha <- 0.05

# use uncorrected tests
results_raw <- t(replicate(B, many_tests(n = n, k = k, mu = mu, alpha = alpha,
                                         correct_mc = FALSE)))
colnames(results_raw) <- c("FWER", "FDR")
colMeans(results_raw)

# use corrected tests
results_mc <- t(replicate(B, many_tests(n = n, k = k, mu = mu, alpha = alpha,
                                        correct_mc = TRUE)))
colnames(results_mc) <- c("FWER", "FDR")
colMeans(results_mc)

