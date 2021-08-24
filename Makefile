README.md: README.Rmd
	R -e "rmarkdown::render('README.Rmd', output_format = 'github_document')"
runshiny:
	R -e "shiny::runApp('Shinyapp', port=3029)"
deployshiny:
	R -e "rsconnect::deployApp('Shinyapp', appName = 'cellmigRation')"