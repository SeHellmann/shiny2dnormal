library(shiny)

# Define UI for Bayesfactor
shinyUI(fluidPage(withMathJax(),
                  theme = "bootstrap.css",

                  title = "Fit-a-normal!",

                  titlePanel("Fit-a-normal!"),

                  fluidRow(
                    column(4,
                           h4("Parameters"),
                           sliderInput("mu1", "Mean of X  (\\(\\mu_{X}\\)):",
                                       min = -2.0, max = 2.0, value = 0.0, step = 0.1, width=NULL),
                           sliderInput("mu2", "Mean of Y  (\\(\\mu_{Y}\\)):",
                                       min = -2.0, max = 2.0, value = 0.0, step = 0.1, width=NULL),
                           sliderInput("sigma1", "SD of X  (\\(\\sigma_{X}\\)):",
                                       min = 0.1, max = 4.0, value = 1.5, step = 0.1, width=NULL),
                           sliderInput("sigma2", "SD of Y  (\\(\\sigma_{Y}\\)):",
                                       min = 0.1, max = 4.0, value = 1.5, step = 0.1, width=NULL),
                           sliderInput("rho", "Correlation  (\\(\\rho\\)):",
                                       min = -1.0, max = 1.0, value = 0.0, step = 0.1, width=NULL),
                           checkboxInput("showTrue", "Show true distribution", value = FALSE),
                           actionButton("reset_it", "Reset & new data"),
                           actionButton("optim", "Let 'optim' find the best fit!")
                    ),
                    column(8,
                           plotOutput("scatterplot")
                    )
                  )# ,
                  #
                  # 	fluidRow(
                  # 		column(4,
                  # 			helpText("RSS-Stats")
                  # 		),
                  # 		column(8,
                  # 			plotOutput("RSS"),
                  # 			helpText("This widgets relies on the RSA package for plotting.")
                  # 		)
                  # 	)
))
