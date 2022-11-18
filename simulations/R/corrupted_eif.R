# function from Alan to corrupt estimates of the efficient influence function
# in a manner from which variance moderation is able to recover
corrupt_eif_sim <- function(eif, delta = 0.2) {
  n <- length(eif)
  rr <- rbinom(1, 1, 0.5)
  uu <- delta * runif(n, 0, abs(eif))
  eif_out <- rr * ((eif < 0) * (eif + uu) + (eif >= 0) * (eif - uu)) +
      (1 - rr) * ((eif < 0) * (eif - uu) + (eif > 0) * (eif + uu))
  eif_out <- eif_out - mean(eif_out)
  return(eif_out)
}
