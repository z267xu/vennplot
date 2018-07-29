#' Data on human encounters with great white sharks.
#'
#' @format A dataset with 65 rows and 11 columns.
#' \describe{
#'   \item{Year}{Years encounter sharks}
#'   \item{Sex}{Sex of victims}
#'   \item{Age}{Age of victims}
#'   \item{Time}{Encounter sharks in AM or PM}
#'   \item{Australia}{Encounter in Australia}
#'   \item{USA}{Encounter in the United States}
#'   \item{Surfing}{Surfing incident}
#'   \item{Scuba}{Scuba-diving incident}
#'   \item{Fatality}{Whether or not there was a fatality}
#'   \item{Injury}{Whether or not there was an injury}
#'   \item{Length}{The length of great white sharks}
#' }
#' @source \url{http://sharkattackinfo.com/shark_attack_news_sas.html}.  Data collected by Professor Pierre-Jerome Bergeron, University of Ottawa.
#' @examples
#' vennplot(disjoint.combinations = sharks, vars = c("Au","USA","Fa","Ti"))
"sharks"
