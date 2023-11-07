library(shiny)

shinyUI(fluidPage(
  titlePanel("Metabolomics Data Reformatting App"),
  sidebarLayout(
    sidebarPanel(
      textInput("directory", "Directory Path", value = ""),
      textInput("runner", "Runner", value = ""),
      textInput("experiment", "Experiment", value = ""),
      textInput("file_names", "File Names (comma-separated, no space)", value = ""),
      textInput("time_value", "Time Point", value = ""),
      actionButton("analyze", "Submit"),
      actionButton("exit", "Exit")
    ),
    mainPanel(
      textOutput("message")
    )
  )
))
