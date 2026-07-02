#' Simulated CDM data with 10 attributes
#'
#' Simulated responses for 2000 examinees generated with the script in
#' `data-raw/simN2000K10.R` using `GDINA::simGDINA()`.
#'
#' @format A list with components:
#' \describe{
#' \item{\code{data}}{A \eqn{2000 \times 30} response matrix.}
#' \item{\code{Q}}{The \eqn{30 \times 10} Q-matrix used for simulation.}
#' \item{\code{true.params}}{The true item category response probabilities used to generate the data.}
#' }
#'
#'
#' @examples
#' data(simN2000K10)
#' str(simN2000K10)
"simN2000K10"

#' Simulated CDM data with 15 attributes
#'
#' Simulated responses for 2000 examinees generated with the script in
#' `data-raw/simN2000K15.R` using `GDINA::simGDINA()`.
#'
#' @format A list with components:
#' \describe{
#' \item{\code{data}}{A \eqn{2000 \times 45} response matrix.}
#' \item{\code{Q}}{The \eqn{45 \times 15} Q-matrix used for simulation.}
#' \item{\code{true.params}}{The true item category response probabilities used to generate the data.}
#' }
#'
#'
#' @examples
#' data(simN2000K15)
#' str(simN2000K15)
"simN2000K15"

#' Simulated CDM data with 20 attributes
#'
#' Simulated responses for 2000 examinees generated with the script in
#' `data-raw/simN2000K20.R` using `GDINA::simGDINA()`.
#'
#' @format A list with components:
#' \describe{
#' \item{\code{data}}{A \eqn{2000 \times 60} response matrix.}
#' \item{\code{Q}}{The \eqn{60 \times 20} Q-matrix used for simulation.}
#' \item{\code{true.params}}{The true item category response probabilities used to generate the data.}
#' }
#'
#'
#' @examples
#' data(simN2000K20)
#' str(simN2000K20)
"simN2000K20"