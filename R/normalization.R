#' Convert counts to counts per million
#'
#' Normalize a CRISPR count matrix so each sample column sums to one million by
#' default. Set `margin = 1` to normalize rows instead.
#'
#' @param counts Numeric matrix or data frame with guides/features in rows and
#'   samples in columns.
#' @param pseudocount Non-negative numeric scalar added to every entry before
#'   normalization.
#' @param margin Integer margin to normalize. Use `2` for columns and `1` for
#'   rows.
#'
#' @return A numeric matrix scaled so each selected margin sums to one million.
#' @examples
#' counts <- matrix(
#'   c(100, 300, 50, 50),
#'   nrow = 2,
#'   dimnames = list(c("guide_1", "guide_2"), c("sample_a", "sample_b"))
#' )
#'
#' counts_per_million(counts)
#' @export
counts_per_million <- function(counts, pseudocount = 0, margin = 2) {
  if (
    !is.numeric(pseudocount) ||
      length(pseudocount) != 1 ||
      is.na(pseudocount) ||
      pseudocount < 0
  ) {
    stop("pseudocount must be a single non-negative number.", call. = FALSE)
  }
  if (!margin %in% c(1, 2)) {
    stop("margin must be 1 for rows or 2 for columns.", call. = FALSE)
  }

  values <- as.matrix(counts)
  if (!is.numeric(values)) {
    stop("counts must be numeric.", call. = FALSE)
  }

  values <- values + pseudocount
  totals <- if (margin == 2) colSums(values) else rowSums(values)
  if (any(totals == 0)) {
    stop("Cannot normalize a margin with total count zero.", call. = FALSE)
  }

  if (margin == 2) {
    sweep(values, 2, totals, "/") * 1e6
  } else {
    sweep(values, 1, totals, "/") * 1e6
  }
}
