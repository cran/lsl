.vec <- function(x) {
  x <- as.matrix(as.vector(x))
  return(x)
}

.vech <- function(x) {
  x <- as.matrix(x[lower.tri(x, diag = T)])
  return(x)
}

.ltri <- function(x) {
  x[upper.tri(x == 1)] = 0
  return(x)
}

.is_est <- function(x){
  x <- is.na(x) | x == 1
  return(x)
}

.is_one <- function(x){
  x <- !is.na(x) & x == 1
  return(x)
}
