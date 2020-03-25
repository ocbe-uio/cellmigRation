library(shiny)
library(DT)

# ==============================================================================
# User interface
# ==============================================================================
ui <- pageWithSidebar(
    headerPanel("CellMigRation"),

    # --------------------------------------------------------------------------
    # Sidebar panel for inputs
    # --------------------------------------------------------------------------
    sidebarPanel(
        fileInput("import", "Import trackedCells file")
    ),

    # --------------------------------------------------------------------------
    # Main panel for displaying outputs
    # --------------------------------------------------------------------------
    mainPanel(
        tabsetPanel(type="tabs",
            tabPanel(
                "Data preview",
                
                DT::dataTableOutput("data")
            )
        )
    )
)

# ==============================================================================
# Server logic
# ==============================================================================
server <- function(input, output) {
    # --------------------------------------------------------------------------
    # Show imported data
    # --------------------------------------------------------------------------
    output$data <- DT::renderDataTable({
        req(input$import)
        df <- get(load(input$import$datapath))
        if (!is(df, "trackedCells")) {
            stop(
                "Wrong data format.",
                " Please load a file containing an object",
                " of class 'trackedCells'.")
        }
        return(df@images$images[[1]])
    })
}

shinyApp(ui, server)