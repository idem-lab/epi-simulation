#' @title Build a within/between-group contact (mixing) matrix
#' 
#' @description
#' Constructs a \eqn{P \times P} contact matrix where each row gives the
#' contact *distribution* (probabilities) for one group. A fraction `within`
#' of contacts stay inside the same group (diagonal), and the remaining
#' \eqn{1 - \text{within}} is shared equally across all other groups
#' (off-diagonals). Each row sums to 1.
#'
#' @param P Integer \eqn{\ge 1}. Number of groups (populations).
#' @param within Numeric in \eqn{[0,1]}. Fraction of contacts that stay within
#'   the same group (placed on the diagonal).
#'
#' @return A \eqn{P \times P} numeric matrix. Rows represent probability
#'   distributions over contact targets and therefore sum to 1.
#'
#' @examples
#' # Three groups with strong within-group contact
#' C1 <- make_contact(P = 3, within = 0.9)
#' C1
#' rowSums(C1)   # all 1
#'
#' # More mixing across groups
#' C2 <- make_contact(P = 4, within = 0.6)
#' C2
#'
#' # Edge cases
#' make_contact(P = 1)        # -> 1x1 matrix [1]
#' make_contact(P = 3, within = 1)  # identity matrix
#' make_contact(P = 3, within = 0)  # all off-diagonal and equal
#'
#' @export
make_contact <- function(P, within = 0.9) {
  # Create a P x P contact/mixing matrix for P groups.
  # - "within" = fraction of contacts that stay within the same group.
  # - The remainder (1 - within) is spread equally across other groups.
  # - Rows are normalized to sum to 1 (each group’s contact pattern is a probability distribution).
  
  stopifnot(P >= 1, within >= 0, within <= 1)
  
  if (P == 1) {
    # Special case: only one group, so all contact is within itself.
    return(matrix(1, 1, 1))
  }
  
  # Start with equal share of "between-group" contact in each off-diagonal entry
  C <- matrix((1 - within) / (P - 1), nrow = P, ncol = P)
  
  # Set diagonals = within-group contact fraction
  diag(C) <- within
  
  # Row-normalize (each row sums to 1)
  C <- sweep(C, 1, rowSums(C), "/")
  
  return(C)
}


