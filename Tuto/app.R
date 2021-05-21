
# load the shiny library and an example app
library(shiny)
# Define server logic required to draw a histogram ----



server <- function(input, output) {
 
}


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Clustering analyis of genes from pairs of C.albicans isolate based on their expression profile"),
  navbarPage("test"),
  sidebarLayout( 
    sidebarPanel("sidebar panel"),
    position = "right",
    mainPanel("main panel")
  )
)
shinyApp(ui = ui, server = server)