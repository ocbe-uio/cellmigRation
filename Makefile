runshiny:
	R -e "shiny::runApp('Shinyapp')"
deployshiny:
	R -e "rsconnect::deployApp('Shinyapp', appName = 'cellmigRation')"