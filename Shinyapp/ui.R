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
				src = "https://raw.githubusercontent.com/ocbe-uio/cellmigRation/master/cell_migration_logo.png",
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
				src = "https://raw.githubusercontent.com/ocbe-uio/cellmigRation/master/cell_migration_logo.png",
				width = "50%"
			),
			h4("Welcome to the cellmigRation Shiny app!"), p(),
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
					"Estimating parameters usually takes several minutes.",
					"Please click the 'Submit' button on the left and wait",
					"for the centroid plot to appear below.",
					h1("Matrix image"),
					p(),
					plotOutput("VisualizeImg"),
				),
				tabPanel("Help!",
					img(
						src = "https://raw.githubusercontent.com/ocbe-uio/cellmigRation/master/cell_migration_logo.png",
						width = "30%"
					),
					h4("Welcome to the cellmigRation Shiny app!"), p(),
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