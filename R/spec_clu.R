#' Spectral Clustering
#'
#' This function performs spectral clustering on a given affinity matrix using
#' one of several normalized Laplacian formulations. The algorithm extracts the
#' first \code{K} eigenvectors of the normalized Laplacian and applies
#' discretization to obtain cluster assignments.
#'
#' @param affinity A square affinity matrix.
#' @param K Integer. The number of clusters.
#' @param type Integer (1, 2, or 3). Specifies the normalization strategy for the
#' Laplacian matrix:
#' \itemize{
#'   \item \code{1}: Unnormalized Laplacian \eqn{L = D - A}
#'   \item \code{2}: Random-walk normalized Laplacian \eqn{D^{-1} L}
#'   \item \code{3}: Symmetric normalized Laplacian \eqn{D^{-1/2} L D^{-1/2}}
#' }
#'
#' @return A vector of cluster labels for each sample.
#' @keywords spectral-clustering Laplacian eigenvectors affinity
#' @export
#'
#' @examples
#' # Suppose A is an affinity matrix
#' # labels <- spec.clu(A, K = 4, type = 3)


#' Discretization of Eigenvectors for Spectral Clustering
#'
#' Internal function used to convert continuous spectral embedding
#' (eigenvectors) into discrete cluster indicator matrix. This follows the
#' algorithm described in spectral clustering literature where a rotation
#' matrix is iteratively refined to minimize the Ncut objective.
#'
#' @param eigenVectors A matrix of continuous eigenvectors (embedding space).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{discrete}: Discrete cluster indicator matrix.
#'   \item \code{continuous}: The normalized continuous embedding.
#' }
#'
#' @keywords internal spectral-clustering discretization



spec.clu <- function(affinity, K, type=3) {
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps
  D = diag(d)
  L = D - affinity
  if (type == 1) {
    NL = Di %*% L %*% Di
  }
  eig = eigen(NL)
  res = sort(abs(eig$values),index.return = TRUE)
  U = eig$vectors[,res$ix[1:K]]
  normalize <- function(x) x / sqrt(sum(x^2))
  if (type == 3) {
    U = t(apply(U,1,normalize))
  }
  eigDiscrete = .discretisation(U)
  eigDiscrete = eigDiscrete$discrete
  labels = apply(eigDiscrete,1,which.max)
  return(labels)
}

.discretisation <- function(eigenVectors) {
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))

  n = nrow(eigenVectors)
  k = ncol(eigenVectors)

  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])

  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }

  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }

  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)

    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]

    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps)
      break

    lastObjectiveValue = NcutValue
    R = V %*% t(U)
  }

  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}
