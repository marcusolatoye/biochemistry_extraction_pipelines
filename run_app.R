library(shiny)

# Set the working directory
setwd("C:/Users/olatoyeo/OneDrive - J.R. Simplot Company/Documents/2023/Biochemistry")

# Run the Shiny app
shinyApp(ui = source("ui.R")$value, server = source("server.R")$value)

