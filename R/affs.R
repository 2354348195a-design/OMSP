#' Construct Affinity Matrix from Input Data
#'
#' This function computes an affinity (similarity) matrix based on pairwise
#' distances between samples. The affinity is calculated using a Gaussian kernel
#' with locally adaptive scaling for each sample.
#'
#' @param x A numeric matrix or data frame where each row represents a sample.
#' @param na.action A function defining how to handle missing values.
#'   Default is \code{na.omit}, which removes rows containing NA.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Removes missing values using \code{na.action}.
#'   \item Computes the pairwise squared Euclidean distance matrix using
#'         fast matrix operations.
#'   \item For each sample, estimates a local scale parameter as the median
#'         of the 5 smallest non-zero distances.
#'   \item Constructs an affinity matrix using a Gaussian kernel:
#'         \deqn{K_{ij} = \exp \left( - \frac{d_{ij}}{s_i s_j} \right)}
#'         where \eqn{d_{ij}} is the squared Euclidean distance and
#'         \eqn{s_i} is the local scaling factor for sample \eqn{i}.
#' }
#'
#' This approach is similar to the adaptive Gaussian kernel used in
#' self-tuning spectral clustering.
#'
#' @return An m × m affinity matrix, where m is the number of samples.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' x <- matrix(rnorm(50), nrow = 10)
#' K <- affs(x)
#' }
#'
#' @export
affs <- function(x, na.action = na.omit) {



  x <- na.action(x)
  rown <- rownames(x)
  x <- as.matrix(x)
  m <- nrow(x)



  s <- rep(0, m)
  dota <- rowSums(x * x) / 2
  dis <- crossprod(t(x))  # dis 是一个 m x m 的矩阵





  for (i in 1:m) {
    dis[i, ] <- 2 * (-dis[i, ] + dota + rep(dota[i], m))
  }



  dis[dis < 0] <- 0



  for (i in 1:m) {
    s[i] <- median(sort(sqrt(dis[i, ]))[1:5])
  }





  km <- exp(-dis / (s %*% t(s)))  # 确保维度正确，避免广播错误

  return(km)
}
