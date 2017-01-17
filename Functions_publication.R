DeBois <- function(Ht, Wt){
  0.007184 * Ht^0.725 * Wt^0.425
}



rmse <- function(x, y)
  sqrt(mean( (x - y)^2, na.rm = TRUE))