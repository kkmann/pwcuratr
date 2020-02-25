start_shiny_app <- function() {
    shiny::runApp(
        glue("{find.package('pwcuratr')}/shiny"),
        launch.browser = TRUE
    )
}
