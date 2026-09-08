#!/bin/bash

echo "Iniciando VVO apps..."

R -e "setwd('/app/VVO_rnaseq'); shiny::runApp('.', host='0.0.0.0', port=8180, launch.browser=FALSE)" &

sleep 3

R -e "setwd('/app/VVO_genomics'); shiny::runApp('.', host='0.0.0.0', port=8181, launch.browser=FALSE)" &

sleep 3

R -e "
library(httpuv)
landing_html <- paste(readLines('/app/landing/index.html', warn = FALSE), collapse = '\n')
app <- list(call = function(req) {
  list(status = 200L, headers = list('Content-Type' = 'text/html'), body = landing_html)
})
startServer('0.0.0.0', 8000, app)
while (TRUE) { service(); Sys.sleep(0.01) }
" &

wait