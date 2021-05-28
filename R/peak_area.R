peak_area <- function(x, y) {
  PA <- sapply(1:(length(x)-1), function(i) {
    (x[i + 1] - x[i])*(y[i + 1] + y[i])
  })
  sum_PA <- 0.5*sum(PA)
  return(sum_PA)
}
