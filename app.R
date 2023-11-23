#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
input <- list(mu1=1, mu2=2, sigma1 = 0.5, sigma2=1, rho=-.2)
library(shiny)

library(mvtnorm)
library(ggplot2)
require(shinyBS)
library(shinyjs)

# Define UI for application that draws a histogram
ui <- fillPage(
  #header=tags$head(tags$style(type='text/css', ".irs-grid-text { font-size: 8pt; }")),
  withMathJax(), useShinyjs(),titlePanel("Two-dimensional normal model"),
  sidebarLayout(position="left",
                mainPanel = mainPanel(plotOutput(outputId = "distPlot", height="350px"),width=8),
                sidebarPanel= sidebarPanel(
                  tags$h4("Model Parameters"),
                         sliderInput("mu1", "Mean of X  (\\(mu_{X}\\)):",
                                     min = -2.0, max = 2.0, value = 1.0, step = 0.1, width=NULL),
                         sliderInput("mu2", "Mean of Y  (\\(mu_{Y}\\)):",
                                     min = -2.0, max = 2.0, value = 1.0, step = 0.1, width=NULL),
                         sliderInput("sigma1", "SD of X  (\\(sigma_{X}\\)):",
                                     min = 0.1, max = 4.0, value = 1.0, step = 0.1, width=NULL),
                         sliderInput("sigma2", "SD of Y  (\\(sigma_{Y}\\)):",
                                     min = 0.1, max = 4.0, value = 1.0, step = 0.1, width=NULL),
                         sliderInput("rho", "Correlation  (\\(rho\\)):",
                                     min = -1.0, max = 1.0, value = 0.0, step = 0.1, width=NULL),
                  width=4
                  #p(actionButton("recalc","Re-run simulation", icon("random")))# HTML(paste("Re-run", "simulation", sep="<br/>"))
                )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
      #observe(input$recalc)
      # empM <- colMeans(Dat)
      # empSigma <- cov(Dat)
      m <- c(input$mu1, input$mu2)
      sigma <- matrix(c(input$sigma1^2,input$sigma1*input$sigma2*input$rho,
                        input$sigma1*input$sigma2*input$rho, input$sigma2^2), nrow=2)
      data.grid <- expand.grid(X1 = seq(-4, 4, length.out=200), X2 = seq(-4, 4, length.out=200))
      q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m, sigma = sigma))
      # loglikelihood <- sum(log(mvtnorm::dmvnorm(Dat, mean = m, sigma = sigma)))
      ggplot() +
        geom_contour(data=q.samp, aes(x=X1, y=X2, z=prob)) +
        #geom_point(data=Dat, aes(x=X1, y=X2))+
        coord_fixed(xlim = c(-4, 4), ylim = c(-4, 4), ratio = 1)+
        ylab("Y")+xlab("X")+
        #annotate(geom="text", x=60, y= 60, hjust=0.5, vjust=0.5, label=as.character(round(loglikelihood, 2)))+
        #annotate(geom="text", x=60, y= 65, hjust=0.5, vjust=0.5, label="log(Likelihood)")+
        theme_minimal()
        # generate bins based on input$bins from ui.R
        # x    <- faithful[, 2]
        # bins <- seq(min(x), max(x), length.out = input$bins + 1)
        #
        # # draw the histogram with the specified number of bins
        # hist(x, breaks = bins, col = 'darkgray', border = 'white',
        #      xlab = 'Waiting time to next eruption (in mins)',
        #      main = 'Histogram of waiting times')
    })
}

# Run the application
shinyApp(ui = ui, server = server)
