library(shiny)
library(tiff)
library(png)

# ==============================================================================
# User interface
# ==============================================================================
ui <- fluidPage(
    titlePanel("CellMigRation"),

    # --------------------------------------------------------------------------
    # Sidebar panel for inputs
    # --------------------------------------------------------------------------
    sidebarPanel(
        fileInput("imported_tiff", "Import TIFF file")
    ),

    # --------------------------------------------------------------------------
    # Main panel for displaying outputs
    # --------------------------------------------------------------------------
    mainPanel(
        h1("Cell tracking"),
        imageOutput("tiff")
    )
)

# ==============================================================================
# Server logic
# ==============================================================================
options(shiny.maxRequestSize=500*1024^2)  # file limit: 500 MB
server <- function(input, output) {
    # --------------------------------------------------------------------------
    # Load imported data
    # --------------------------------------------------------------------------
    output$tiff <- renderImage({
        req(input$imported_tiff)
        filename <- normalizePath(file.path(input$imported_tiff$datapath))
        list(
            src=filename,
            alt="there should be an image here",
            width=400,
            height=400
        )
    })
}

shinyApp(ui, server)