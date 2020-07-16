# ==============================================================================
# Loading necessary packages
# ==============================================================================
library(shiny)
library(tiff)
library(png)
library(cellmigRation)
# ==============================================================================
# Defining the user interface
# ==============================================================================
ui <- fluidPage(
	titlePanel(
		title = div(
			img(
				src = "https://raw.githubusercontent.com/ocbe-uio/CellMigRation/master/cell_migration_logo.png",
				width = "5%"
			),
			"cellmigRation"
		),
		windowTitle = "cellmigRation Shiny"
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
			checkboxInput("invert_background", "Invert background"),
			hr(),
			h4("Metadata"),
			textInput(
				inputId = "project_name",
				label = "Project name",
				placeholder = "Unnamed project"
			),
			fluidRow(
				column(
					width = 7,
					textInput(
						inputId = "project_condition",
						label = "Condition",
						placeholder = "Treated, dose, etc."
					)
				),
				column(
					width = 5,
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
					width = 3,
					numericInput(
						inputId = "frame_duration",
						label = "Frame duration",
						value = 1,
						min = 0

					)
				),
				column(
					width = 3,
					selectInput(
						inputId = "frame_unit",
						label = "unit",
						choices = c("s", "min", "h")
					)
				),
				column(
					width = 3,
					numericInput(
						inputId = "pixel_size",
						label = "Pixel size",
						value = 1,
						min = 0
					)
				),
				column(
					width = 3,
					selectInput(
						inputId = "pixel_unit",
						label = "unit",
						choices = c("nm", "µm", "mm", "cm")
					)
				)
			)
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
					"Automated parameter estimation" = "auto",
					"Use values below" = "user"
				)
			),
			fluidRow(
				column(4, numericInput("lnoise", "L-noise", 0, step = .1)),
				column(4, numericInput("diamenter", "Diameter", 0, step = .1)),
				column(4, numericInput("threshold", "Threshold", 0, step = .1))
			),
			numericInput("num_threads", "Number of CPU threads to use", 2, 1),
			actionButton("fit_model", "Submit")
		),
		# ----------------------------------------------------------------------
		# Cell tracking
		# ----------------------------------------------------------------------
		conditionalPanel(
			condition = "input.fit_model",
			hr(),
			h3("3. Cell tracking"),
			numericInput(
				inputId = "max_disp",
				label = "Max displacement",
				value = 16
			),
			actionButton("track_cells", "Track cells")
		),
		# ----------------------------------------------------------------------
		# Output data
		# ----------------------------------------------------------------------
		conditionalPanel(
			condition = "input.track_cells",
			hr(),
			h3("4. Output data"),
			actionButton("extract_trajectories", "Extract Trajectories"),
			actionButton("extract_summary", "Extract Summary")
		)
	),
	# --------------------------------------------------------------------------
	# Main panel for displaying outputs
	# --------------------------------------------------------------------------
	mainPanel(
		conditionalPanel(
			condition = "!output.slider",
			img(
				src = "https://raw.githubusercontent.com/ocbe-uio/CellMigRation/master/cell_migration_logo.png",
				width = "50%"
			),
			h4("Welcome to the CellMigRation Shiny app!"), p(),
			"Please load a proper TIFF file using the 'Browse' button on the",
			"left. After the file is loaded, you will be presented with more",
			"options and a help page."
		),
		conditionalPanel(
			condition = "output.slider",
			tabsetPanel(
				tabPanel("Original image", imageOutput("image_frame")),
				tabPanel("Processed image", plotOutput("processed_image")),
				tabPanel("Model estimation",
					h1("Matrix image"), p(),
					plotOutput("VisualizeImg")
				),
				tabPanel("Help!",
					img(
						src = "https://raw.githubusercontent.com/ocbe-uio/CellMigRation/master/cell_migration_logo.png",
						width = "30%"
					),
					h4("Welcome to the CellMigRation Shiny app!"), p(),
					"- Frame selection on the slider and autoplay is disabled", "for the processed image for performance purposes.",
					"Use the 'Previous/Next frame' buttons instead.",
					"You can also select a frame on the slider and then press",
					"one of the aforementioned buttons to load a slide."
				),
				id = "post_load"
			)

		),
		br(),br(),br(),br(),br(),br(),br(),br(),
		conditionalPanel(
			condition = "input.track_cells",
			tabsetPanel(
				tabPanel(
					"Trajectories"
				),
				tabPanel(
					"Summary"
				)
			)
		)
	)
)
# ==============================================================================
# Defining the server logic
# ==============================================================================
options(shiny.maxRequestSize = 1024*1024^2)  # file limit: 1 GB
server <- function(input, output, session) {
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
	# Displaying data
	# --------------------------------------------------------------------------
	output$processed_image <- renderPlot({
		req(input$imported_tiff)
		x1 <- CellMigRation::LoadTiff(
			tiff_file  = input$imported_tiff$datapath,
			experiment = input$project_name,
			condition  = input$project_condition,
			replicate  = input$replicate
		)
		# X1 <- readRDS(...) # ASK: what is this about?
		# Store in variables for now
		time_var <- input$frame_duration
		res_var <- input$pixel_size
		invert_background <- input$invert_background
		VisualizeImg(x1@images$images[[frame$out]], las = 1, main = paste("Stack num.", frame$out))
	})
	# --------------------------------------------------------------------------
	# Fitting model
	# --------------------------------------------------------------------------
	observeEvent(input$fit_model, {
		updateTabsetPanel(
			session,
			inputId = "post_load",
			selected = "Model estimation"
		) #FIXME: not selecting anymore
		# FIXME: text below not printing to UI
		h1("Estimating parameters. This often takes some minutes. Please wait.")
		# Automated parameter optimization
		x1 <- LoadTiff(
			tiff_file  =  input$imported_tiff$datapath,
			experiment = input$project_name,
			condition  = input$project_condition,
			replicate  = input$replicate
		) # TODO: DRY: move this and L:291 to one reactive function
		# TODO: move x1, b, pk and cnt to renderText?
		x1 <- OptimizeParams(tc_obj = x1, threads = input$num_threads)
		# Retrieve optimized values
		lnoise    <- x1@optimized$auto_params$lnoise
		diameter  <- x1@optimized$auto_params$diameter
		threshold <- x1@optimized$auto_params$threshold
	# 	# Visualize Centroids
		# TODO: add UI feedback (R output or message)
		b <- CellMigRation:::bpass(
			image_array = x1@images$images[[frame$out]],
			lnoise = lnoise,
			lobject = diameter,
			threshold = threshold
		)
		pk <- cellmigRation:::pkfnd(
			im = b,
			th = threshold,
			sz = cellmigRation:::NextOdd(diameter)
		)
		cnt <- cellmigRation:::cntrd(
			im = b,
			mx = pk,
			sz = cellmigRation:::NextOdd(diameter)
		)
		# TODO: return the output of the following to the user
		output$VisualizeImg <- renderPlot({
			VisualizeImg(
				img_mtx = b, las = 1, main = paste0("Stack num. ", frame$out)
			)
			cellmigRation:::VisualizeCntr(
				centroids = cnt, width_px = ncol(b), height_px = nrow(b)
			)
		})
	})
	# --------------------------------------------------------------------------
	# Tracking cells
	# --------------------------------------------------------------------------
	eventReactive(input$track_cells, {
		# TODO: track cells using cellmigRation functions
	})
}
# ==============================================================================
# Running the server
# ==============================================================================
shinyApp(ui, server)
