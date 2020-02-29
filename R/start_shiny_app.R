#' Interactive pathway curation
#'
#' Start the interactive shiny app fo pathway cluster curation.
#'
#' @examples
#' \dontrun{start_shiny_app()}
#'
#' @include zzz.R
#'
#' @export
start_shiny_app <- function() {
    shiny::runApp(
        glue("{find.package('pwcuratr')}/shiny"),
        launch.browser = TRUE
    )
}
