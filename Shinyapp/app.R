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
        h1("Cell tracking"),
        sliderInput("frame", "Frame:", min=1, max=10, value=1), # TODO: custom max
        h1("Data summary"),
        "Number of images:", textOutput("tot_frames"),
        "Dimensions of images:", textOutput("dimdata"),
        "Data structure:", DT::dataTableOutput("data")
    )
)

# ==============================================================================
# Server logic
# ==============================================================================
server <- function(input, output) {
    # --------------------------------------------------------------------------
    # Load imported data
    # --------------------------------------------------------------------------
    data <- reactive({
        req(input$import)
        data <- get(load(input$import$datapath))
        if (!is(data, "trackedCells")) {
            stop(
                "Wrong data format.",
                " Please load a file containing an object",
                " of class 'trackedCells'.")
        }
        data <- data@images$images
        return(data)
    })
    tot_frames <- reactive(length(data()))
    output$tot_frames <- renderText(tot_frames())
    output$dimdata <- renderText(dim(data()[[1]]))
    output$data <- DT::renderDataTable(data.frame(data()[[1]]))
    output$images <- renderPlot({
        # if (!is(data, "trackedCells")) {
        #     browser()
        # } else {
            plot(rnorm(100), rnorm(100))
        # }
    })
}

shinyApp(ui, server)