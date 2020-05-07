# ==============================================================================
# Loading necessary packages
# ==============================================================================
library(shiny)
library(tiff)
library(png)
# ==============================================================================
# Defining the user interface
# ==============================================================================
ui <- fluidPage(
	titlePanel("CellMigRation"),

	# --------------------------------------------------------------------------
	# Sidebar panel for inputs
	# --------------------------------------------------------------------------
	sidebarPanel(
		fileInput("imported_tiff", "Import TIFF file"),
		uiOutput("slider"),
		fixedRow(
			column(width=2, uiOutput("prev")),
			column(width=4, uiOutput("play")),
			column(width=1, uiOutput("nxt"))
		)
	),

	# --------------------------------------------------------------------------
	# Main panel for displaying outputs
	# --------------------------------------------------------------------------
	mainPanel(
		h1("Cell tracking"),
		imageOutput("image_frame")
	)
)
# ==============================================================================
# Defining the server logic
# ==============================================================================
options(shiny.maxRequestSize=500*1024^2)  # file limit: 500 MB
server <- function(input, output) {
	frame <- reactiveValues(out=1, autoplay=FALSE)
	# --------------------------------------------------------------------------
	# Load imported data
	# --------------------------------------------------------------------------
	image <- reactive({
		req(input$imported_tiff)
		filename <- normalizePath(file.path(input$imported_tiff$datapath))
		filepath <- gsub(
			x=input$imported_tiff$datapath,
			pattern=".\\.tif$",
			replacement=""
		)
		split_tiff <- readTIFF(filename, all=TRUE, convert=TRUE)
		split_png <- list()
		for (i in seq_along(split_tiff)) {
			writePNG(
				image = split_tiff[[i]],
				target = paste0(
					filepath, formatC(i, flag="0", width=5), '.png'
				)
			)
		}
		file_list <- list.files(filepath, pattern="*.png")
		return(list(path = filepath, name = file_list))
	})
	# --------------------------------------------------------------------------
	# Create image controls
	# --------------------------------------------------------------------------
	tot_frames <- reactive(length(image()$name))
	output$tot_frames <- renderText(tot_frames())
	output$slider <- renderUI(
		sliderInput(
			"frameSelector", "Frame:", min=1, max=tot_frames(), value=frame$out,
			step=1
		)
	)
	output$prev <- renderUI({
		if (is.null(tot_frames())) return()
		actionButton("prev", "<")
	})
	output$nxt <- renderUI({
		if (is.null(tot_frames())) return()
		actionButton("nxt", ">")
	})
	output$play <- renderUI({
		if (is.null(tot_frames())) return()
		actionButton("play", "Autoplay")
	})
	# --------------------------------------------------------------------------
	# Determining slide to show
	# --------------------------------------------------------------------------
	observeEvent(input$prev, {
		frame$out <- max(input$frameSelector - 1, 1)
	})
	observeEvent(input$nxt, {
		frame$out <- min(input$frameSelector + 1, tot_frames())
	})
	observeEvent(input$play, {
		frame$autoplay <- TRUE
	})
	src_output <- reactive({
		filename <- image()$name[input$frameSelector]
		out <- paste(paste0(image()$path, filename))
		return(out)
	})
	# --------------------------------------------------------------------------
	# Render imported data
	# --------------------------------------------------------------------------
	output$image_frame <- renderImage({
		req(input$imported_tiff)
		list(
			src=src_output(),
			alt="image not found",
			width="60%"
		)
	}, deleteFile=FALSE)
}
# ==============================================================================
# Running the server
# ==============================================================================
shinyApp(ui, server)