# Synthetic presentation-only harness. Does not load data sources or call APIs.
source("src/interpretation_module.R")
source("src/interpretation_backend.R")
source("tests/fixtures/summary_fixture.R")
shiny::addResourcePath("summary-assets", normalizePath("www"))
ui <- shiny::fluidPage(shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "summary-assets/frontend.css")),
  shiny::h3("Variant interpretation · Evidence summary"),
  shiny::radioButtons("scenario", "Synthetic validation scenario", c("Steps not run" = "before", "All steps complete" = "after"), inline = TRUE),
  shiny::uiOutput("steps"), shiny::radioButtons("view", "Evidence view", c("Table", "Map"), inline = TRUE),
  shiny::div(id = "summary-test-container", style = "width:1260px;max-width:100%;margin-inline:auto", shiny::uiOutput("evidence")), shiny::tags$script(type = "application/json", id = "fixture-expected",
    shiny::HTML(jsonlite::toJSON(lapply(c("before", "after"), function(s) {
      states <- if (s == "before") interp_request_states() else interp_request_states(list(), list(), list())
      f <- summary_stress_fixture(states)
      list(state = s, text = unlist(lapply(f$items, `[[`, "text")), chars = sum(nchar(unlist(f$summary))),
        finding_lengths = nchar(f$summary$finding), evidence_lengths = unlist(lapply(f$items, function(x) nchar(x$text))))
    }), auto_unbox = TRUE))))
server <- function(input, output, session) {
  states <- shiny::reactive(if (input$scenario == "before") interp_request_states() else interp_request_states(list(), list(), list()))
  output$steps <- shiny::renderUI(shiny::div(class = "regshiny-request-states", interpretation_step_tags(states())))
  output$evidence <- shiny::renderUI({
    f <- summary_stress_fixture(states())
    if (input$view == "Map") interpretation_summary_map(f$summary, f$items, states(), f$identifiers) else
      interpretation_summary_cards(f$summary, items = f$items, states = states(), identifiers = f$identifiers,
        baseline = if (input$scenario == "after") summary_stress_fixture()$summary else NULL)
  })
}
shiny::runApp(shiny::shinyApp(ui, server), host = "127.0.0.1", port = as.integer(Sys.getenv("SUMMARY_TEST_PORT", "8765")), launch.browser = FALSE)
