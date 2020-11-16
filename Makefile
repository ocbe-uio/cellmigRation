runshiny:
	R -e "shiny::runApp('Shinyapp', port=3029)"
deployshiny:
	R -e "rsconnect::deployApp('Shinyapp', appName = 'cellmigRation')"