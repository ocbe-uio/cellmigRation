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
			column(width=4, uiOutput("prev")),
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
options(shiny.maxRequestSize=1024*1024^2)  # file limit: 1 GB
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
				image=split_tiff[[i]],
				target=paste0(
					filepath, formatC(i, flag="0", width=5), '.png'
				)
			)
		}
		file_list <- list.files(filepath, pattern="*.png")
		return(list(path=filepath, name=file_list))
	})
	# --------------------------------------------------------------------------
	# Creating image controls
	# --------------------------------------------------------------------------
	tot_frames <- reactive(length(image()$name))
	output$tot_frames <- renderText(tot_frames())
	output$slider <- renderUI(
		sliderInput(
			inputId="frameSelector", label="Frame select:",
			min=1, max=tot_frames(), value=frame$out,
			step=1, animate=animationOptions(
				interval=200,
				playButton="Autoplay",
				pauseButton="Pause"
			),
		)
	)
	output$prev <- renderUI({
		if (is.null(tot_frames())) return()
		actionButton("prev", "Previous frame")
	})
	output$nxt <- renderUI({
		if (is.null(tot_frames())) return()
		actionButton("nxt", "Next frame")
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
	src_output <- reactive({
		filename <- image()$name[input$frameSelector]
		out <- paste0(image()$path, filename)
		return(out)
	})
	# --------------------------------------------------------------------------
	# Render imported data
	# --------------------------------------------------------------------------
	output$image_frame <- renderImage(
		expr={
			req(input$imported_tiff)
			list(
				src=src_output(),
				alt="image not found",
				width="60%"
			)
		},
		deleteFile=FALSE
	)
}
# ==============================================================================
# Running the server
# ==============================================================================
shinyApp(ui, server)
