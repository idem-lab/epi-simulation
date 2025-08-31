#' check_contact
#' @title Validate a contact/mixing matrix
#' @description Ensures `C` is square, finite, and non-negative.
#' @param C A numeric matrix.
#' @return TRUE (invisibly) or an error if invalid.
#' @examples
#' # check_contact(diag(3))
check_contact <- function(C) {
  if (!is.matrix(C)) stop("C must be a matrix.")
  if (nrow(C) != ncol(C)) stop("C must be square (same rows and cols).")
  if (any(!is.finite(C))) stop("C contains non-finite values.")
  if (any(C < 0)) stop("C contains negative entries.")
  invisible(TRUE)
}

#' normalize_contact
#' @title Normalize a contact matrix by rows or columns
#' @description Divides each row (or column) by its sum so rows (or cols) sum to 1.
#' @param C A numeric matrix (will be validated).
#' @param method One of "row" or "col".
#' @return A normalized matrix of the same size as `C`.
#' @examples
#' # Cn <- normalize_contact(matrix(1,3,3), "row"); rowSums(Cn)
normalize_contact <- function(C, method = c("row", "col")) {
  method <- match.arg(method)
  check_contact(C)
  if (method == "row") {
    rs <- rowSums(C); if (any(rs == 0)) stop("Row with zero sum cannot be normalized.")
    return(C / rs)
  } else {
    cs <- colSums(C); if (any(cs == 0)) stop("Column with zero sum cannot be normalized.")
    return(t(t(C) / cs))
  }
}
