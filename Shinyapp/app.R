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
        filepath <- gsub(
            x=input$imported_tiff$datapath,
            pattern=".\\.tif$",
            replacement=""
        )
        split_tiff <- readTIFF(filename, all=TRUE)
        split_png <- list()
        for (i in seq_along(split_tiff)) {
            writePNG(split_tiff[[i]], paste0(filepath, i, '.png'))
        }
        file_list <- list.files(filepath, pattern="*.png")
        list(
            src=paste0(filepath, file_list[1]),
            alt="there should be an image here",
            width="50%"
        )
    })
}

shinyApp(ui, server)