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

