# Validate and align the recorded four-age, five-replicate design.
read_age_samples <- function(matrix_file, manifest_file, tissue) {
  x <- read.delim(matrix_file, row.names = 1, check.names = FALSE)
  m <- read.delim(manifest_file, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("sample", "tissue", "age")
  if (!all(required %in% names(m))) stop("Sample manifest requires sample, tissue and age")
  m <- m[m$tissue == tissue, , drop = FALSE]
  if (nrow(m) != 20L || anyNA(m[, required]) || anyDuplicated(m$sample)) stop("Require 20 unique samples for the selected tissue")
  if (!setequal(m$age, c(20, 50, 80, 100)) || any(table(m$age) != 5L)) stop("Require five samples at each of 20, 50, 80 and 100 weeks")
  if (anyDuplicated(colnames(x)) || !all(m$sample %in% colnames(x))) stop("Expression columns do not match the sample manifest")
  m <- m[order(m$age, m$sample), , drop = FALSE]
  x <- as.matrix(x[, m$sample, drop = FALSE])
  if (!is.numeric(x)) stop("Expression matrix must be numeric")
  rownames(m) <- m$sample
  list(matrix = x, samples = m)
}
