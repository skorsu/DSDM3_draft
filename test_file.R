zi_test <- (1:5) * 10
g <- rbinom(5, 1, 0.5)
bk <- (1:5)/5

log_marginal(zi = zi_test, gamma_ik = g, beta_k = bk)

log(factorial(sum(zi_test[1:4]))) - sum(log(factorial(zi_test[1:4])))
log(factorial(zi_test))
