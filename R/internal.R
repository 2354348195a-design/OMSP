#' Label Propagation for Semi-supervised Subtype Prediction
#'
#' This function predicts labels (subtypes) for new samples using graph-based
#' semi-supervised learning. It supports two commonly used methods:
#' local and global consistency (method = 0), and classical label propagation
#' (method > 0).
#'
#' @param W A similarity or affinity matrix.
#' @param Y0 A label indicator matrix. Rows correspond to samples;
#'   columns correspond to classes. Labeled samples should have one-hot
#'   rows; unlabeled samples contain all zeros.
#' @param method An integer specifying the algorithm to use:
#'   \itemize{
#'     \item 0: Local and Global Consistency (Zhou et al.)
#'     \item >0: Classical Label Propagation
#'   }
#'
#' @details
#' For method = 0, the closed-form solution is:
#' \deqn{Y = (1 - \alpha) (I - \alpha P)^{-1} Y_0}
#' where \eqn{P} is the row-normalized affinity matrix.
#'
#' For method > 0, the method iteratively updates:
#' \deqn{Y^{(t+1)} = P Y^{(t)}}
#' with the labels of known samples clamped at each iteration.
#'
#' @return A matrix of predicted soft labels.
#'
#' @keywords internal
.csPrediction <- function(W,Y0,method){
  alpha=0.9;
  P= W/rowSums(W)
  if(method==0){
    Y= (1-alpha)* solve( diag(dim(P)[1])- alpha*P)%*%Y0;
  } else {
    NLabel=which(rowSums(Y0)==0)[1]-1;
    Y=Y0;
    for (i in 1:1000){
      Y=P%*%Y;
      Y[1:NLabel,]=Y0[1:NLabel,];
    }
  }
  return(Y);
}

#' Discretize Continuous Eigenvectors in Spectral Clustering
#'
#' Internal function for converting continuous eigenvectors to discrete
#' cluster indicators using Yu & Shi (2003) algorithm.
#'
#' @param eigenVectors A matrix of continuous eigenvectors.
#'
#' @return A list containing:
#' \item{discrete}{Discrete cluster indicator matrix}
#' \item{continuous}{Row-normalized eigenvectors}
#'
#' @keywords internal
.discretisation <- function(eigenVectors){
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  mini <- function(x) { i = which(x == min(x)); return(i[1]) }
  c = matrix(0,n,1)
  for (j in 2:k){
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  lastObjectiveValue = 0
  for (i in 1:20){
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) break
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
  }
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}

#' Convert Eigenvectors to One-hot Cluster Assignment
#'
#' Internal helper function for discretization.
#'
#' @param eigenVector A numeric matrix of eigenvectors.
#'
#' @return A binary matrix representing cluster assignments.
#'
#' @keywords internal
.discretisationEigenVectorData <- function(eigenVector){
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) { i = which(x == max(x)); return(i[1]) }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1
  return(Y)
}

#' Construct Dominant Neighborhood Graph
#'
#' Internal function to retain only the top KK largest elements from each row
#' of a similarity matrix and normalize rows.
#'
#' @param xx A similarity matrix.
#' @param KK Number of dominant neighbors to keep per row.
#'
#' @return A row-normalized matrix with only top-KK neighbors preserved.
#'
#' @keywords internal
.dominateset <- function(xx,KK=20){
  zero <- function(x) {
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  normalize <- function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){ A[i,] = zero(xx[i,]); }
  return(normalize(A))
}

#' Compute Mutual Information Between Two Label Vectors
#'
#' Internal function for calculating mutual information.
#'
#' @param x A vector of discrete labels.
#' @param y A second vector of discrete labels.
#'
#' @return Numeric value of mutual information (bits).
#'
#' @keywords internal
.mutualInformation <- function(x, y){
  classx <- unique(x)
  classy <- unique(y)
  nx <- length(x)
  ncx <- length(classx)
  ncy <- length(classy)
  probxy <- matrix(NA, ncx, ncy)
  for (i in 1:ncx){ for (j in 1:ncy){ probxy[i, j] <- sum((x == classx[i]) & (y == classy[j])) / nx } }
  probx <- matrix(rowSums(probxy), ncx, ncy)
  proby <- matrix(colSums(probxy), ncx, ncy, byrow=TRUE)
  result <- sum(probxy * log(probxy / (probx * proby), 2), na.rm=TRUE)
  return(result)
}

#' Compute Shannon Entropy of a Label Vector
#'
#' @param x A vector of categorical labels.
#'
#' @return Shannon entropy (base-2).
#'
#' @export
entropy <- function(x){
  class <- unique(x)
  nx <- length(x)
  nc <- length(class)
  prob <- rep.int(NA, nc)
  for (i in 1:nc){ prob[i] <- sum(x == class[i])/nx }
  result <- -sum(prob * log(prob, 2))
  return(result)
}

#' Replicate Matrix (MATLAB \code{repmat} Equivalent)
#'
#' Internal helper function for matrix replication.
#'
#' @param X A scalar, vector, or matrix to be replicated.
#' @param m Number of row replications.
#' @param n Number of column replications.
#'
#' @return Replicated matrix.
#'
#' @keywords internal
.repmat = function(X,m,n){
  if (is.null(dim(X))) { mx = length(X); nx = 1 } else { mx = dim(X)[1]; nx = dim(X)[2] }
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

#' Convert Values to Proportions
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector where each element is divided by sum(x).
#'
#' @export
perc<-function(x){
  sapply(1:length(x), function(i){ x[i]/sum(x) })
}

