library(shiny)
library(mvtnorm)
library(ggplot2)
require(gridExtra)
# construct a new environment for global data sharing
dat <- new.env()
#input=list(sigma1 = 1, sigma2=1, mu1=0,mu2=0, rho=0)
reset_values <- function() {
  dat$i <- 0
  dat$mus <- runif(2, -1.5, 1.5)
  dat$sigmas <- runif(2, 0.3, 1.5)
  dat$rho <- runif(1, -0.9, 0.9)
  dat$Sigma <- matrix(c(dat$sigmas[1]^2, dat$sigmas[1]*dat$sigmas[2]*dat$rho,
                        dat$sigmas[1]*dat$sigmas[2]*dat$rho,dat$sigmas[2]^2), nrow=2)
  sim <- rmvnorm(200, mean=dat$mus, sigma=dat$Sigma)
  sim <- as.data.frame(sim)
  names(sim) <- c("X1", "X2")
  dat$sim <- sim

  # # These variables store indices of the optimal fit
  # dat$l1 <- lm(dat$y~dat$x)
  dat$loglik_opt <- -sum(log(mvtnorm::dmvnorm(dat$sim, mean =dat$mus, sigma = dat$Sigma)))

  # These variables store the sequence of (manual) optimization steps
  dat$loglik.manual <- c()

  dat$optim.pars <- data.frame()
  dat$optim.pars.plot <- data.frame()
}
reset_values()


# The minimization function for optim
minimize_loglik <- function(X) {
    Sigma <- matrix(c(exp(X[3])^2, exp(X[3])*exp(X[4])*((pnorm(X[5])-0.5)*2),
                      exp(X[3])*exp(X[4])*((pnorm(X[5])-0.5)*2),exp(X[4])^2), nrow=2)
    loglik.optim <- -sum(log(mvtnorm::dmvnorm(dat$sim, mean = c(X[1],X[2]), sigma = Sigma)))
  dat$i <- dat$i +1
  if (dat$i %% 5 ==0 ) {
    # save parameter set for intermediate results
    dat$optim.pars <- rbind(dat$optim.pars,
                            c(c(X[1], X[2], exp(X[3]), exp(X[4]), (pnorm(X[5])-0.5)*2),
                              loglik.optim))
  }
  #if (is.infinite(loglik.optim)) loglik.optim <- 1e+12
  return(loglik.optim)
}
minimize_loglik2 <- function(X) {

  Sigma <- matrix(c(X[3]^2, X[3]*X[4]*X[5],
                    X[3]*X[4]*X[5],X[4]^2), nrow=2)
  loglik.optim <- -sum(log(mvtnorm::dmvnorm(dat$sim, mean = c(X[1],X[2]), sigma = Sigma)))

  # save parameter set for intermediate results
  dat$optim.pars <- rbind(dat$optim.pars, c(X, loglik.optim))
  #if (is.infinite(loglik.optim)) loglik.optim <- 1e+12
  return(loglik.optim)
}
negloglik <- function(X) {
  Sigma <- matrix(c(X[3]^2, X[3]*X[4]*X[5],
                    X[3]*X[4]*X[5],X[4]^2), nrow=2)
  return(-sum(log(mvtnorm::dmvnorm(dat$sim, mean = c(X[1],X[2]), sigma = Sigma))))
}

min_function <- minimize_loglik
min_function2 <- minimize_loglik2


shinyServer(function(input, output, session) {

  output$scatterplot <- renderPlot({

    input$reset_it	# this call triggers this function whenever the reset button is clicked.
    input$optim

    # if the optim has been called, we want to show the sequence of optimization steps
    # always increase the to-be-plotted sequence until it reaxhes the end
    if (nrow(dat$optim.pars.plot) < nrow(dat$optim.pars)) {
      dat$optim.pars.plot <- dat$optim.pars[1:(nrow(dat$optim.pars.plot)+1), ]
      invalidateLater(150, session)
    }


    dat$loglik.manual <- c(dat$loglik.manual,
                           negloglik(c(input$mu1, input$mu2, input$sigma1, input$sigma2, input$rho))) # residual sum of squares
    # ---------------------------------------------------------------------
    # Scatterplot

    if (input$showTrue == TRUE) {
      title <- paste0("Raw data + theoretical distribution + true distribution")
    } else {
      title <-  paste0("Raw data + theoretical distribution")
    }
    m <- c(input$mu1, input$mu2)
    sigma <- matrix(c(input$sigma1^2,input$sigma1*input$sigma2*input$rho,
                      input$sigma1*input$sigma2*input$rho, input$sigma2^2), nrow=2)
    data.grid <- expand.grid(X1 = seq(-4, 4, length.out=200), X2 = seq(-4, 4, length.out=200))
    q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m, sigma = sigma))

    p1 <- ggplot() +
      geom_contour(data=q.samp, aes(x=X1, y=X2, z=prob)) +
      geom_point(data=dat$sim, aes(x=X1, y=X2))+
      coord_fixed(xlim = c(min(dat$sim[,1],-4), max(dat$sim[,1],4)),
                  ylim = c(min(dat$sim[,2],-4), max(dat$sim[,2],4)),  ratio = 1)+
      ylab("Y")+xlab("X")+
      #annotate(geom="text", x=60, y= 60, hjust=0.5, vjust=0.5, label=as.character(round(loglikelihood, 2)))+
      #annotate(geom="text", x=60, y= 65, hjust=0.5, vjust=0.5, label="log(Likelihood)")+
      theme_minimal()

    # optim abline
    if (nrow(dat$optim.pars.plot) > 0) {
      X <- c(t(dat$optim.pars.plot[nrow(dat$optim.pars.plot),]))
      # X[3] <- exp(X[3])
      # X[4] <- exp(X[4])
      # X[5] <- (pnorm(X[5])-0.5)*2
      Sigma2 <- matrix(c(X[3]^2, X[3]*X[4]*X[5],
                      X[3]*X[4]*X[5],X[4]^2), nrow=2)
      p.samp <- cbind(data.grid, prob=mvtnorm::dmvnorm(data.grid, mean =c(X[1],X[2]), sigma = Sigma2))
      p1 <- p1 + geom_contour(data=p.samp, aes(x=X1, y=X2, z=prob), color="darkblue", linetype="dashed")
    }
    if (input$showTrue == TRUE) {

      t.samp <- cbind(data.grid, prob=mvtnorm::dmvnorm(data.grid, mean =dat$mus, sigma = dat$Sigma))
      p1 <- p1 + geom_contour(data=t.samp, aes(x=X1, y=X2, z=prob), color="darkgreen", linetype="dashed",
                              linewidth=1.2)
    }

    ## ======================================================================
    ## loglik. plot
    ## ======================================================================

    XLIM <- c(1, max(10, length(dat$loglik.manual), nrow(dat$optim.pars.plot)))
    if (nrow(dat$optim.pars.plot) > 0) {
      YLIM <- c(dat$loglik_opt*0.9,
                max(c(dat$loglik.manual, dat$optim.pars.plot[ ,6])[!is.infinite(c(dat$loglik.manual, dat$optim.pars.plot[ ,6]))])*1.2)
    } else {
      YLIM <- c(dat$loglik_opt*0.9, max(dat$loglik.manual[!is.infinite(dat$loglik.manual)])*1.2)
    }

    # ---------------------------------------------------------------------
    # Manual optimization steps
    p2 <- ggplot(data.frame(x=1:length(dat$loglik.manual), y=dat$loglik.manual), aes(x=x, y=y))+
      geom_point(shape=20, col="red")+
      geom_line(aes(group=1),col="red")+xlim(XLIM)+ylim(YLIM)+
      xlab("Step")+ylab("negative log-likelihood")+ggtitle("Negative log-likelihood\nSmaller values = better fit")+
      theme_minimal()
    # plot(1:length(dat$loglik.manual), dat$loglik.manual, type="o", pch=20, col="red",
    #      xlim=XLIM, ylim=YLIM, xlab="Step", ylab="negative log-likelihood", main="Negative log-likelihood\nSmaller values = better fit")

    # ---------------------------------------------------------------------
    # 'optim' optimization steps
    if (nrow(dat$optim.pars.plot) > 0) {
      p2 <- p2 + geom_line(data=data.frame(x=1:nrow(dat$optim.pars.plot), y=dat$optim.pars.plot[ ,6]),
                           aes(group=2),col="darkblue")+
        geom_point(data=data.frame(x=1:nrow(dat$optim.pars.plot), y=dat$optim.pars.plot[ ,6]), shape=20, col="darkblue")
    }
    p2 <- p2 +geom_hline(yintercept= dat$loglik_opt, col="red", lty="dashed")
    gridExtra::grid.arrange(p1, p2, ncol=2)
  })


  # Observe reset button and do something it has been clicked
  observe({
    if (input$reset_it == 0) return()
    isolate({
      reset_values()
      updateSliderInput(session, "mu1", value = 0)
      updateSliderInput(session, "mu2", value = 0)
      updateSliderInput(session, "sigma1", value = 1.5)
      updateSliderInput(session, "sigma2", value = 1.5)
      updateSliderInput(session, "rho", value = 0)
      invalidateLater(1, session)
    })
  })


  # Observe optim button and do something it has been clicked
  observe({
    if (input$optim == 0) return()

    isolate({
      # reset the optimization history
      dat$optim.pars <- data.frame()
      dat$optim.pars.plot <- data.frame()

      # optimize parameters. Reduce relative tolerance (otherwise the optimizer stays quite long around the final line)
      result <- optim(par = c(mean(dat$sim[,1]), mean(dat$sim[,2]), log(sd(dat$sim[,1])), log(sd(dat$sim[,2])), qnorm(cor(dat$sim[,1], dat$sim[,2])/2+0.5))+
                        rnorm(5, 0, 0.3),
                      min_function, control=list(reltol=.002, maxit=200))
      result$par <- with(result, c(par[1], par[2], exp(par[3]), exp(par[4]), (pnorm(par[5]))-0.5)*2)
      # result <- minqa::bobyqa(par = c(0, 0, 0.8, 0.8, 0), min_function2,
      #                 lower = c(-3, -3, 0, 0 ,-1), upper=c(5, 5, 5, 5, 1),
      #                 control = list(rhobeg=1,rhoend=0.0001, maxfun=300))

      # show the sequence of optimization steps
      invalidateLater(1, session)

    })
  })

})
