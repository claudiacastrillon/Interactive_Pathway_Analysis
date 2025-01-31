library(shiny)
library(DT)  # For DataTables

# Define the user interface for the app
ui <- fluidPage(
    titlePanel("Interactive Pathway Analysis"),
    sidebarLayout(
        sidebarPanel(
            # DataTable to display pathway analysis results
            dataTableOutput("pathwayTable")
        ),
        mainPanel(
            # Network visualization
            visNetworkOutput("visNetwork")
        )
    )
)
