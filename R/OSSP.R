#' Multi-omics Cancer Subtyping Pipeline
#'
#' This function performs a complete multi-omics subtyping workflow:
#' data standardization, similarity matrix construction, self-diffusion denoising,
#' feature extraction via Laplacian eigenvectors, clustering, and survival analysis.
#'
#' @param mydatGE Numeric matrix of gene expression data (samples in rows, genes in columns).
#' @param mydatME Numeric matrix of methylation data.
#' @param mydatMI Numeric matrix of microRNA data.
#' @param survival Data frame with at least two columns: 'Survival' (time) and 'Death' (event).
#' @param K Optional integer. Number of clusters for k-means. Default is 7.
#' @param lap_type Integer. Type of Laplacian: 1=unnormalized, 2=random walk, 3=symmetric normalized. Default is 3.
#' @param diffusion_iter Integer. Number of iterations for self-diffusion denoising. Default is 4.
#' @param seed Integer. Random seed for reproducibility. Default is 11111.
#'
#' @return A list containing:
#' \item{lab}{Cluster assignments for each sample.}
#' \item{rd}{Combined feature matrix used for clustering.}
#' \item{sil}{Silhouette object for clustering quality visualization.}
#'
#' @examples
#' \dontrun{
#' result <- multiomics_pipeline(mydatGE, mydatME, mydatMI, survival)
#' }
#'
#' @export
OSSP <- function(mydatGE, mydatME, mydatMI, survival,
                                K = 7, lap_type = 3, diffusion_iter = 4,
                                seed = 11111) {

  # ------------------------------
  # 1. 数据标准化和 NA 处理
  # ------------------------------
  mydatGE2 <- apply(mydatGE, 2, function(x) (x - mean(x)) / sd(x))
  mydatME2 <- apply(mydatME, 2, function(x) (x - mean(x)) / sd(x))
  mydatMI2 <- apply(mydatMI, 2, function(x) (x - mean(x)) / sd(x))

  mydatGE2[is.na(mydatGE2)] <- 0
  mydatMI2[is.na(mydatMI2)] <- 0

  d <- list(mydatGE2, mydatME2, mydatMI2)

  # ------------------------------
  # 2. 主程序：构建特征矩阵
  # ------------------------------
  ptm <- proc.time()

  ul <- lapply(1:3, function(i) {
    x <- affs(d[[i]])
    affinity <- self.diffusion(x, diffusion_iter)

    # 构建拉普拉斯矩阵
    d_vec <- rowSums(affinity)
    d_vec[d_vec == 0] <- .Machine$double.eps
    D <- diag(d_vec)
    L <- D - affinity
    if (lap_type == 1) {
      NL <- L
    } else if (lap_type == 2) {
      Di <- diag(1 / d_vec)
      NL <- Di %*% L
    } else if (lap_type == 3) {
      Di <- diag(1 / sqrt(d_vec))
      NL <- Di %*% L %*% Di
    }

    # 特征分解
    eig <- eigen(NL)
    res <- sort(abs(eig$values), index.return = TRUE)
    e <- res$x[1:15]
    sa <- sapply(1:14, function(j) abs(e[j+1] - e[j]))
    w <- which(sa == max(sa))
    K_feat <- w + 3  # 特征向量个数

    U <- eig$vectors[, res$ix[1:K_feat]]

    # 对特征向量归一化
    normalize <- function(x) x / sqrt(sum(x^2))
    if (lap_type == 3) {
      U <- t(apply(U, 1, normalize))
    }

    return(U)
  })

  rd <- do.call(cbind, ul)
  cat("Feature extraction done. Time elapsed:", proc.time() - ptm, "\n")

  # ------------------------------
  # 3. 聚类分析
  # ------------------------------
  set.seed(seed)
  labx <- kmeans(rd, K)
  lab <- labx$cluster

  # ------------------------------
  # 4. 生存分析
  # ------------------------------
  library(survcomp)
  clins <- Surv(as.numeric(survival$Survival), as.numeric(survival$Death))

  plot.KM(
    clins,
    as.integer(lab),
    palette = "lanonc",
    xlab = "Follow up(Months)",
    ylab = "OS(pro.)",
    risk.table = TRUE
  )

  # ------------------------------
  # 5. 聚类质量评估
  # ------------------------------
  x <- affs(rd)
  sil <- silhouette_SimilarityMatrix(lab, x)
  plot(sil)

  # 返回结果
  return(list(
    lab = lab,
    rd = rd,
    sil = sil
  ))
}
