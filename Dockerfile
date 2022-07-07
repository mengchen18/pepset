FROM rocker/shiny:4.2.1

RUN mkdir ~/projects

RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('mengchen18/pepset', upgrade = 'never')"

CMD ["R", "-e", "shiny::runApp(system.file(package = 'pepset', 'app', 'app.R'), host = '0.0.0.0', port = 3838)"]
