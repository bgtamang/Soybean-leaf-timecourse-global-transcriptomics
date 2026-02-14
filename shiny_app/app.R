# GmJAG1 RNA-Seq Shiny Dashboard (Phase3)
# Launch script

# Source global, ui, and server
source("global.R")
source("ui.R")
source("server.R")

# Run the app
shinyApp(ui = ui, server = server)
