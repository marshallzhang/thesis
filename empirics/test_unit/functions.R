q.unit = function(theta, r.diffusion, theta.tilde, r.nu, M, steps, exact = F) {
  set.seed(1)
  theta1 = theta[1]
  theta2 = theta[2]
  expectation = 0
  bridges = if (exact) {
    exact.bridges(1, r.nu, M, r.diffusion, theta.tilde, steps)
  } else {
    gen.nu.bridge(r.nu, r.diffusion, theta.tilde, M, steps)
  }
  us = sample(1:steps, M, replace = T)
  expectation = sum(rowSums(rowDiffs(bridges) * (theta1 - theta2 * bridges[,1:(steps-1)]))) - sum((1/2) * (theta1 - theta2 * bridges[cbind(1:M, us)])^2)
  -expectation / M
}
