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
			h3("2. Model fit"),
			radioButtons(
				inputId = "who_estimates_parms",
				label = "Model parameters",
				choices = list(
					"Automated parameter estimation" = "auto",
					"Use values below" = "user"
				),
				selected = "user"
			),
			fluidRow(
				column(4, numericInput("lnoise", "L-noise", 1, step = .1)),
				column(4, numericInput("diameter", "Diameter", 10, step = .1)),
				column(4, numericInput("threshold", "Threshold", 5, step = .1))
			),
			numericInput(
				inputId = "num_threads",
				label   = "Number of CPU threads to use",
				value   = parallel::detectCores() - 2,
				min     = 1,
				max     = parallel::detectCores()
			),
			# TODO: add some "please wait" text below. Alternatively, output all
			# terminal messages on a dedicated box/tab
			actionButton("fit_model", "Fit model")
		),
		# ----------------------------------------------------------------------
		# Cell tracking
		# ----------------------------------------------------------------------
		conditionalPanel(
			condition = "output.step == '3'",
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
			condition = "output.step == '4'",
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
			h1("Welcome to the cellmigRation Shiny app!"), p(),
			"Please load a suitable TIFF file using the 'Browse' button on the",
			"left. After the file is loaded, you will be presented with more",
			"options and a help page."
		),
		conditionalPanel(
			condition = "output.slider",
			tabsetPanel(
				tabPanel("1. Original image", imageOutput("image_frame")),
				tabPanel("1. Processed image", plotOutput("processed_image")),
				tabPanel("2. Model fit",
					h1("Instructions"),
					"Estimating parameters usually takes several minutes.",
					"Please click the 'Submit' button on the left and wait",
					"for the centroid plot to appear below.",
					"This operation can take several minutes to complete.",
					p(),
					"When running this application locally, you may want to",
					"adjust the number of CPU threads to match your processor",
					"so operations run in parallel and produce output faster.",
					"We recommend leaving one thread for your computer's",
					"internal operations. In other words, if your CPU has $X$",
					"threads, we recommend setting the number of CPU threads",
					"to use by this app to $X - 1$.",
					h1("Matrix image"),
					p(),
					plotOutput("VisualizeImg"),
				),
				tabPanel("3. Tracking cells",
					plotOutput("VisualizeImgStep3"),
				),
				tabPanel("About/Help",
					img(
						src = "https://raw.githubusercontent.com/ocbe-uio/cellmigRation/master/cell_migration_logo.png",
						width = "30%"
					),
					h1("Welcome to the cellmigRation Shiny app!"), p(),
					"cellmigRation is an R package for tracking cells",
					"and analyzing their trajectories. This Shiny app",
					"provides a graphical version of the package, allowing for",
					"easier access to some of its most important features. To",
					"access all the features of cellmigRation, please use the",
					a(
						href="https://bioconductor.org/packages/cellmigRation/",
						"R package published on Bioconductor"
					), ".",
					h2("Troubleshooting"), p(),
					h3("Frame selection is not working as expected"), p(),
					"Frame selection on the slider and autoplay is disabled",
					"for the processed image for performance purposes.",
					"Use the 'Previous/Next frame' buttons instead.",
					"You can also select a frame on the slider and then press",
					"one of the aforementioned buttons to load a slide.",
					h3("Form fields are too tiny"), p(),
					"Try increasing the screenspace the app is taking,",
					"for example by setting the window to fullscreen.",
					h3("I press 'Fit model' and nothing happens"), p(),
					"This is a very computer-intensive process,",
					"so please allow at lease a few minutes (potentially,",
					"hours, depending on how large your file is and the",
					"hardware it's running on) for the software to run.",
					"If you are running this app locally, you may increase the",
					"number of CPU threads to be used. You can run",
					"'detectCores()' in R to see how high you can go, though",
					"it's a good idea to use a number slightly lower than the",
					"one output by R.",
					h3("Something else is wrong"), p(),
					"Please file a bug report on",
					a(
						href="https://github.com/ocbe-uio/CellMigRation/issues",
						"https://github.com/ocbe-uio/CellMigRation/issues"
					)
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
		),
		h3("Message box"),
		# TODO: improve message box to show server messages. Alternatively, use
		# outputOptions(output, "step", suspendWhenHidden=FALSE)
		# (Source: https://stackoverflow.com/a/39827259/1169233)
		verbatimTextOutput("step"),


	)
)
