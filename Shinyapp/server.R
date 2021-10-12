# Loading necessary packages and options ---------------------------------------
library(remotes)
library(cellmigRation)
options(shiny.maxRequestSize = 1024*1024^2)  # file limit: 1 GB

# Defining the server logic ----------------------------------------------------
server <- function(input, output, session) {
	# Reactive values ----------------------------------------------------------
	frame <- reactiveValues(out = 1)
	x     <- reactiveValues(x1 = NULL, x2 = NULL)
	parms <- reactiveValues(lnoise = NULL, diameter = NULL, threshold = NULL)
	step  <- reactiveValues(current=0)
	# Load imported data -------------------------------------------------------
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
	# Creating image controls --------------------------------------------------
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
	# Determining slide to show ------------------------------------------------
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
	# Render imported data -----------------------------------------------------
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
	# 1. Displaying data -------------------------------------------------------
	output$processed_image <- renderPlot({
		req(input$imported_tiff)
		x$x1 <- cellmigRation::LoadTiff(
			tiff_file  = input$imported_tiff$datapath,
			experiment = input$project_name,
			condition  = input$project_condition,
			replicate  = input$replicate
		)
		# FIXME #68: doesn't reselect if frame is changed on the slider (only buttons work)
		# X1 <- readRDS(...) # ASK: what is this about?
		# Store in variables for now
		time_var <- input$frame_duration
		res_var <- input$pixel_size
		invert_background <- input$invert_background
		VisualizeImg(
			img_mtx = x$x1@images$images[[frame$out]],
			las = 1,
			main = paste("Stack num.", frame$out)
		)
	})
	# Model and tracking -------------------------------------------------------------
	observeEvent(input$fit_model, {
		# Actually fitting model (or using user values) ------------------------
		x$x1 <- LoadTiff(
			tiff_file  = input$imported_tiff$datapath,
			experiment = input$project_name,
			condition  = input$project_condition,
			replicate  = input$replicate
		) # TODO: DRY: move this and L:291 to one reactive function?
		if (input$who_estimates_parms == "auto") {
			message(Sys.time(), " - Optimizing model parameters. Please wait")
			x$x1 <- OptimizeParams(tc_obj = x$x1, threads = input$num_threads)
			parms$lnoise    <- x$x1@optimized$auto_params$lnoise
			parms$diameter  <- x$x1@optimized$auto_params$diameter
			parms$threshold <- x$x1@optimized$auto_params$threshold
		} else {
			message(Sys.time(), " - Using model with given parameters")
			parms$lnoise    <- input$lnoise
			parms$diameter  <- input$diameter
			parms$threshold <- input$threshold
		}
		# Switching active tab -------------------------------------------------
		updateTabsetPanel( # tab changing happens here, even if moved up or down
			session,
			inputId = "post_load",
			selected = "Model and tracking"
		)
		# Rendering plot -------------------------------------------------------
		output$VisualizeImg <- renderPlot({
			message(Sys.time(), " - Performing bandpass")
			b <- cellmigRation:::bpass(
				image_array = x$x1@images$images[[frame$out]],
				lnoise = parms$lnoise,
				lobject = parms$diameter,
				threshold = parms$threshold
			)
			message(Sys.time(), " - Finding signal peaks")
			pk <- cellmigRation:::pkfnd(
				im = b,
				th = parms$threshold,
				sz = cellmigRation:::NextOdd(parms$diameter)
			)
			message(Sys.time(), " - Calculating centroids")
			cnt <- cellmigRation:::cntrd(
				im = b,
				mx = pk,
				sz = cellmigRation:::NextOdd(parms$diameter)
			)
			# Visualize Centroids
			message(Sys.time(), " - Plotting")
			VisualizeImg(
				img_mtx = b,
				las = 1,
				main = paste0("Stack num. ", frame$out)
			)
			cellmigRation:::VisualizeCntr(
				centroids = cnt, width_px = ncol(b), height_px = nrow(b)
			)
			step$current <- 3
			message(Sys.time(), " - Ready for cell tracking")
		})
		output$step <- renderText(step$current)
	})
	# 3. Cell tracking ---------------------------------------------------------
	observeEvent(input$track_cells, {
		message(Sys.time(), " - Tracking cells. Please wait")
		x$x2 <- cellmigRation:::CellTracker(
			tc_obj     = x$x1,
			lnoise     = parms$lnoise,
			diameter   = parms$diameter,
			threshold  = parms$threshold,
			maxDisp    = input$max_disp,
			threads    = input$num_threads,
			show_plots = FALSE,
			verbose    = FALSE
		)
		# Switching active tab -------------------------------------------------
		output$VisualizeImg <- renderPlot({
			message(Sys.time(), " - Plotting tracks")
			VisualizeImg(
				img_mtx = x$x2@proc_images$images[[frame$out]],
				las = 1,
				main = paste0("Stack num. ", frame$out)
			)
			cellmigRation:::VisualizeCntr(
				centroids = x$x2@centroids[[frame$out]],
				width_px = ncol(x$x2@proc_images$images[[frame$out]]),
				height_px = nrow(x$x2@proc_images$images[[frame$out]])
			)
			message(Sys.time(), " - Ready for step 4")
			step$current <- 4
		})
	})
	# 4. Output data ----------------------------------------------------------
	# TODO #64: after step 3 is done, call visualizeCellTracks() to plot
	observeEvent(input$extract_trajectories, {
		message(Sys.time(), " - Extracting trajectories")
		tracks_df <- getTracks(x$x2) # FIXME: NULL
		step$current <- 5
		message(Sys.time(), " - Trajectories extracted")
		print(str(tracks_df)) #TEMP
	})
	observeEvent(input$extract_summary, {
		message(Sys.time(), " - Extracting summary")
		x$x2 <- ComputeTracksStats(x$x2) #FIXME: Warning: Error in seq.default: feil forteikn i «by»-argument
		cell_summary <- getCellsStats(x$x2)
		step$current <- 6
		message(Sys.time(), " - Summary extracted")
		print(str(cell_summary)) #TEMP
	})
}
