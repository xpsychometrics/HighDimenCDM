#' @importFrom stats runif aggregate
#' @importFrom utils getFromNamespace
NULL

gdina_l2m <- function(...) {
	getFromNamespace("l2m", "GDINA")(...)
}

gdina_lik_nr <- function(...) {
	getFromNamespace("LikNR", "GDINA")(...)
}

gdina_match_matrix <- function(...) {
	getFromNamespace("matchMatrix", "GDINA")(...)
}
