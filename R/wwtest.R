#' Wilcoxon--Wigner homogeneity Test
#'
#' This function tests whether the independent entries in a symmetric matrix are identically distributed.
#' It first transforms the entry values to the corresponding normalized ranks,
#' then computes the first eigenvalue or eigenvector of the transformed matrix. The test statistic is constructed based on
#' the eigenvalue or eigenvector and compared to the standard normal distribution. The function then returns test result
#' including the test statistic, alternative hypothesis and p-value.
#'
#' @param data_matrix A numeric symmetric matrix with \code{nrow > 2}.
#' @param method A character string, either \code{"eigenvalue"} or \code{"eigenvector"}, indicating which test to use.
#' @return An object of class \code{"htest"} containing the test statistic, leading eigenvalue of the transformed matrix, matrix size, and p-value.
#' \item{data.name}{The name of the input data matrix.}
#' \item{statistic}{The test statistic.}
#' \item{parameter}{The size of the matrix (number of rows/columns).}
#' \item{p.value}{The p-value associated with the test statistic.}
#' \item{method}{The name of the test. Specify which alternative (eigenvalue test or eigenvector test) is used.}
#' \item{alternative}{The alternative hypothesis: "The entries in the matrix do not come from the same distribution".}
#' @importFrom stats pnorm
#' @export
#' @examples
#' mat <- matrix(c(1, 2, 3, 2, 1, 1, 3, 1, 1), 3, 3)
#' wwtest(mat)
wwtest <- function(data_matrix, method = 'eigenvalue') {
  # 1. Check if the input is a square matrix and contains no NA values
  if (!is.matrix(data_matrix) || nrow(data_matrix) != ncol(data_matrix)) {
    stop("The input must be a square matrix.")
  }
  if (any(is.na(data_matrix))) {
    stop("The input matrix contains NA values.")
  }

  # 2. Check if the matrix is symmetric
  if (!all(data_matrix == t(data_matrix))) {
    stop("The input must be a symmetric matrix.")
  }

  n = nrow(data_matrix)
  N = n*(n-1)/2

  # 2. Transform the upper triangular entries to ranks
  upper_tri_values <- data_matrix[upper.tri(data_matrix, diag = FALSE)]
  ranked_values <- rank(upper_tri_values)
  transformed_matrix <- matrix(0, nrow = n, ncol = n)
  transformed_matrix[upper.tri(data_matrix, diag = FALSE)] <- ranked_values/(N+1)
  transformed_matrix = transformed_matrix + t(transformed_matrix)

  # 3. Compute the first eigenvalue
  if (method == 'eigenvalue'){
    first_eigenvalue <- eigen(transformed_matrix, only.values = TRUE)$values[1]
  } else {
    first_eigenvector <- eigen(transformed_matrix)$vectors[,1]
  }

  # 4. Construct the test statistic
  sigma2 = 1/12 - 1/(6*(N+1))
  mu1 = 1/2*(n-1) + 2*sigma2
  tilde_sigma2 = 8*sigma2^2/n

  if (method == 'eigenvalue'){
    test_statistic <- (first_eigenvalue - mu1)/(sqrt(tilde_sigma2))
  } else {
    test_statistic <- n*(sum(first_eigenvector)/sqrt(n) - 1 + 1/(6*n))/(sqrt(tilde_sigma2))
  }

  if (method == 'eigenvalue'){
    met_name <- "Wilcoxon--Wigner homogeneity test: eigenvalue test"
  } else {
    met_name <- "Wilcoxon--Wigner homogeneity test: eigenvector test"
  }

  # 5. Calculate the p-value
  p_value <- 2 * pnorm(abs(test_statistic), mean = 0, sd = 1, lower.tail = FALSE)


  # 6. Return the test result in a format similar to other R tests
  result <- list(
    statistic = c(TestStatistic = test_statistic),
    parameter = c(MatrixSize = nrow(data_matrix)), # matrix size (nrow or ncol)
    p.value = p_value,
    method = met_name,
    alternative = "The entries in the matrix are not identically distributed",
    data.name = deparse(substitute(data_matrix))
  )
  class(result) <- "htest"

  return(result)
}
