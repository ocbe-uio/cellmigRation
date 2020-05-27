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
	titlePanel(
		title = "CellMigRation Shiny App",
		windowTitle = "CellMigRation Shiny App"
	),
	# --------------------------------------------------------------------------
	# Sidebar panel for inputs
	# --------------------------------------------------------------------------
	sidebarPanel(
		# ----------------------------------------------------------------------
		# Loading TIFF image
		# ----------------------------------------------------------------------
		h3("1. Data loading"),
		fileInput("imported_tiff", "Import TIFF file"),
		uiOutput("slider"),
		fluidRow(
			column(width = 6, uiOutput("prev")),
			column(width = 6, uiOutput("nxt")),
		),
		# ----------------------------------------------------------------------
		# Metadata
		# ----------------------------------------------------------------------
		conditionalPanel(
			condition = "output.slider",
			hr(),
			textInput(
				inputId = "project_name",
				label = "Project name",
				placeholder = "Unnamed project"
			),
			fluidRow(
				column(
					width = 6,
					textInput(
						inputId = "project_condition",
						label = "Condition",
						placeholder = "Treated, dose, etc."
					)
				),
				column(
					width = 6,
					numericInput(
						inputId = "replicate",
						label = "Replicate",
						value = 1,
						min = 1
					)
				)
			),
			fluidRow(
				column(
					width = 8,
					numericInput(
						inputId = "pixel_size",
						label = "Pixel size",
						value = 1,
						min = 0
					)
				),
				column(
					width = 4,
					selectInput(
						inputId = "pixel_unit",
						label = "unit",
						choices = c("nm", "µm", "mm", "cm")
					)
				)
			),
			fluidRow(
				column(
					width = 6,
					numericInput(
						inputId = "frame_duration",
						label = "Frame duration",
						value = 1,
						min = 0

					)
				),
				column(
					width = 4,
					selectInput(
						inputId = "frame_unit",
						label = "unit",
						choices = c("s", "min", "h")
					)
				)
			),
		),
		# ----------------------------------------------------------------------
		# Model parameters
		# ----------------------------------------------------------------------
		conditionalPanel(
			condition = "output.slider",
			hr(),
			h3("2. Model estimation"),
			radioButtons(
				inputId = "who_estimates_parms",
				label = "Model parameters",
				choices = list(
					"Estimate with CellMigRation" = "auto",
					"Use values below" = "user"
				)
			),
			numericInput("parm1", "Parm1", 0),
			numericInput("parm2", "Parm2", 0),
			actionButton("fit_model", "Submit")
		),
		# ----------------------------------------------------------------------
		# Cell tracking
		# ----------------------------------------------------------------------
		conditionalPanel(
			condition = "input.fit_model",
			hr(),
			h3("3. Cell tracking"),
			actionButton("track_cells", "Track cells")
		)
	),
	# --------------------------------------------------------------------------
	# Main panel for displaying outputs
	# --------------------------------------------------------------------------
	mainPanel(
		conditionalPanel(
			condition = "!output.slider",
			img(
				src = "https://raw.githubusercontent.com/ocbe-uio/CellMigRation/master/CellMigRationLogo.png",
				width = "30%"
			)
		),
		conditionalPanel(
			condition = "output.slider",
			imageOutput("image_frame")
		),
		br(),br(),br(),br(),br(),br(),br(),br(),
		conditionalPanel(
			condition = "input.track_cells",
			tabsetPanel(
				tabPanel(
					"Trajectories"
				),
				tabPanel(
					"Calculations"
				)
			)
		)
	)
)
# ==============================================================================
# Defining the server logic
# ==============================================================================
options(shiny.maxRequestSize = 1024*1024^2)  # file limit: 1 GB
server <- function(input, output) {
	# --------------------------------------------------------------------------
	# Reactive values
	# --------------------------------------------------------------------------
	frame <- reactiveValues(out = 1)
	# --------------------------------------------------------------------------
	# Load imported data
	# --------------------------------------------------------------------------
	image <- reactive({
		req(input$imported_tiff)
		filename <- normalizePath(file.path(input$imported_tiff$datapath))
		filepath <- gsub(
			x = input$imported_tiff$datapath,
			pattern = ".\\.tif$",
			replacement = ""
		)
		split_tiff <- readTIFF(filename, all = TRUE, convert = TRUE)
		split_png <- list()
		for (i in seq_along(split_tiff)) {
			writePNG(
				image = split_tiff[[i]],
				target = paste0(
					filepath, formatC(i, flag = "0", width = 5), '.png'
				)
			)
		}
		file_list <- list.files(filepath, pattern = "*.png")
		return(list(path = filepath, name = file_list))
	})
	# --------------------------------------------------------------------------
	# Creating image controls
	# --------------------------------------------------------------------------
	tot_frames <- reactive(length(image()$name))
	output$tot_frames <- renderText(tot_frames())
	output$slider <- renderUI(
		sliderInput(
			inputId = "frameSelector", label = "Frame select:",
			min = 1, max = tot_frames(), value = frame$out, step = 1,
			animate = animationOptions(
				interval = 200,
				playButton = "Autoplay",
				pauseButton = "Pause"
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
		expr = {
			req(input$imported_tiff)
			list(
				src = src_output(),
				alt = "image not found",
				height = 500
			)
		},
		deleteFile = FALSE
	)
	# --------------------------------------------------------------------------
	# Fitting model
	# --------------------------------------------------------------------------
	eventReactive(input$fit_model, {
		# TODO: fit model using CellMigRation functions
	})
	# --------------------------------------------------------------------------
	# Tracking cells
	# --------------------------------------------------------------------------
	eventReactive(input$track_cells, {
		# TODO: track cells using CellMigRation functions
	})
}
# ==============================================================================
# Running the server
# ==============================================================================
shinyApp(ui, server)
