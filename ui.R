library(shiny)

shinyUI(fluidPage(
  tags$style(HTML("
    .sidebar {
      width: 350px; /* Adjust the width to your preference */
      background-color: #DAA520; /* Set your desired color here #f8f9fa */
    }
    .main {
      margin-left: 360px; /* Adjust the margin to make room for the wider sidebar */
    }
  ")),
  titlePanel("Metabolomics Data Reformatting App"),
  sidebarLayout(
    sidebarPanel(
      width = 5,  # Set the width of the sidebarPanel
      textInput("directory", "Directory Path Format \n C:/Users/username/OneDrive - J.R. Simplot Company/Documents", value = ""),
      textInput("runner", "Runner", value = ""),
      textInput("experiment", "Experiment", value = ""),
      fileInput("files", "Upload CSV Files", multiple = TRUE),
      #textInput("file_names", "File Names (comma-separated, no space)", value = ""),
      textInput("time_value", "Time Point", value = ""),
      textInput("replications", "Number of Reps", value = ""),
      actionButton("analyze", "Submit"),
      actionButton("exit", "Exit")
    ),
    mainPanel(
      textOutput("message")
    )
  )
))
