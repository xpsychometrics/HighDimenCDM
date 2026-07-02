#' Mental health CDM data (Tan2023)
#'
#' In this study, we analyze the responses of 719 college students
#' (34.6\% men, 83.9\% White, and 16.3\% first-year students) to 40 items
#' measuring mental health problems such as alcohol-related problems,
#' anxiety, hostility, and depression.
#'
#' The data were previously analyzed in Tan et al. (2022). The items
#' measure four attributes: alcohol-related problems, anxiety, hostility, and depression.
#'
#' @format A list of responses and Q-matrix with components:
#' \describe{
#' \item{\code{dat}}{Responses of 719 participants to 40 items.}
#' \item{\code{Q}}{The \eqn{40 \times 4} Q-matrix.}
#' }
#'
#' @author Yuxuan Mei, The University of Minnesota, \email{mei00060@umn.edu}
#'
#' @examples
#' \dontrun{
#' data(Tan2023)
#' str(Tan2023)
#' fit <- stEM(
#'   dat = Tan2023$dat,
#'   Q = Tan2023$Q
#'   )
#'}
#'
#' @references
#' Tan, Z., de la Torre, J., Ma, W., Huh, D., Larimer, M. E., & Mun, E.-Y. (2022). A tutorial  on cognitive diagnosis modeling for characterizing mental health symptom profiles using existing item responses. Prevention Science.
"Tan2023"
