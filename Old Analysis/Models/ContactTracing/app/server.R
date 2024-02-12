library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("Cairo"); library("shinyjs")


Server <- function(input, output) {

  #Generation Time
  GenTime <- function(T2, R0) {
    G = T2 * ((R0-1)/log(2))
    return(G)
  }
  
  #Model Betas
  beta <- function(R0) {
    gamma <- 1/(GenTime(3.3,2.8))
    return(R0*gamma)
  } 
  
  #Linear Increase in contact tracing
  lin_cont <- function(x,m,comp) {
    y = m*(x) + 0
    return(y)
  }
  
  SIRS1 <- function(time, init, parms) {
    
    beta <- beta(parms["R0"])
    
    incidence <- beta*init["I"]*init["S"]*parms["N"]
    inctime <- c(time, incidence)
    inc_cont1 <<- rbind(inc_cont1, inctime)
    inc_cont1 <<- inc_cont1[!duplicated(inc_cont1$X0), ]
    
    dS = - beta*init["I"]*init["S"] 
    dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
    dR = parms["gamma"]*init["I"]
    
    return(list(c(dS, dI, dR)))
  }
  SIRS1SIS <- function(time, init, parms) {
    
    beta <- beta(parms["R0"])
    
    incidence <- beta*init["I"]*init["S"]*parms["N"]
    inctime <- c(time, incidence)
    inc_cont1 <<- rbind(inc_cont1, inctime)
    inc_cont1 <<- inc_cont1[!duplicated(inc_cont1$X0), ]
    
    dS = parms["gamma"]*init["I"] - beta*init["I"]*init["S"] 
    dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
    
    return(list(c(dS, dI)))
  }
  
  #SIRS Model To Identify where contact tracing = incidence
  SIRS2 <- function(time, init, parms) {
    
    beta <- beta(parms["R0"])
    linear <<- lin_cont(time, parms["ramp_rate"])
    
    incidence <- beta*init["I"]*init["S"]*parms["N"]
    inctime <- c(time, incidence, linear, beta, linear/ (incidence))
    inc_cont <<- rbind(inc_cont, inctime)
    inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
    
    beta <- beta(parms["R0"]) * (1-(inc_cont[nrow(inc_cont),3]/(inc_cont[nrow(inc_cont),2]))*parms["efficacy"])
    
    dS = - beta*init["I"]*init["S"] 
    dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
    dR = parms["gamma"]*init["I"]
    
    return(list(c(dS, dI, dR)))
  }
  SIRS2SIS <- function(time, init, parms) {
    
    beta <- beta(parms["R0"])
    linear <<- lin_cont(time, parms["ramp_rate"])
    
    incidence <- beta*init["I"]*init["S"]*parms["N"]
    inctime <- c(time, incidence, linear, beta, linear/ (incidence))
    inc_cont <<- rbind(inc_cont, inctime)
    inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
    
    beta <- beta(parms["R0"]) * (1-(inc_cont[nrow(inc_cont),3]/(inc_cont[nrow(inc_cont),2]))*parms["efficacy"])
    
    dS = parms["gamma"]*init["I"] - beta*init["I"]*init["S"] 
    dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
    
    return(list(c(dS, dI)))
  }
  
  #SIRS Model used for Final Model Output
  SIRS3 <- function(time, init, parms) {
    
    timeequal = timeequal
    
    if(time >= timeequal) { # if time is equal
      
      beta <- beta(parms["R0"])*(1-(1*parms["efficacy"]))
      incidence <- beta*init["I"]*init["S"]*parms["N"]
      
      inctime <- c(time, incidence, incidence, beta, 1)
      inc_cont <<- rbind(inc_cont, inctime)
      inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
      
      dS = - beta*init["I"]*init["S"]
      dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
      dR = parms["gamma"]*init["I"] 
      
    } else{ #if inc is over 1000
      
      beta <- beta(parms["R0"])
      linear <<- lin_cont(time, parms["ramp_rate"])
      
      incidence <- beta*init["I"]*init["S"]*parms["N"]
      inctime <- c(time, incidence, linear, beta, 
                   linear/ (incidence))
      inc_cont <<- rbind(inc_cont, inctime)
      inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
      
      beta <- beta(parms["R0"]) * (1-(inc_cont[nrow(inc_cont),3]/(inc_cont[nrow(inc_cont),2]))*parms["efficacy"])
      
      dS = - beta*init["I"]*init["S"] 
      dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
      dR = parms["gamma"]*init["I"]
      
    }
    
    return(list(c(dS, dI, dR)))
  }
  SIRS3SIS <- function(time, init, parms) {
    
    timeequal = timeequal
    
    if(time >= timeequal) { # if time is equal
      
      beta <- beta(parms["R0"])*(1-(1*parms["efficacy"]))
      incidence <- beta*init["I"]*init["S"]*parms["N"]
      
      inctime <- c(time, incidence, incidence, beta, 1)
      inc_cont <<- rbind(inc_cont, inctime)
      inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
      
      dS = parms["gamma"]*init["I"] - beta*init["I"]*init["S"] 
      dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
      
    } else{ #if inc is over 1000
      
      beta <- beta(parms["R0"])
      linear <<- lin_cont(time, parms["ramp_rate"])
      
      incidence <- beta*init["I"]*init["S"]*parms["N"]
      inctime <- c(time, incidence, linear, beta, 
                   linear/ (incidence))
      inc_cont <<- rbind(inc_cont, inctime)
      inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
      
      beta <- beta(parms["R0"]) * (1-(inc_cont[nrow(inc_cont),3]/(inc_cont[nrow(inc_cont),2]))*parms["efficacy"])
      
      dS = parms["gamma"]*init["I"] - beta*init["I"]*init["S"] 
      dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
      
    }
    
    return(list(c(dS, dI)))
  }
  
  
  contact_trac_run <- function(efficacy, ramp_rate, R0, target, iinput, rinput) {
    
    #Parameters
    N <- 5.5*10^6
    init <- c(S = (N-(iinput + rinput))/N, I = iinput/N, R = rinput/N)
    
    parms = c(gamma = 1/(GenTime(3.3,2.8)), 
              N = 5.5*10^6,
              efficacy = efficacy,
              ramp_rate = ramp_rate, 
              R0 = R0,
              target = target)
    
    #Intersect
    inc_cont <<- data.frame(0, beta(parms["R0"])*init["I"]*init["S"]*parms["N"] , 0, 0, 0)
    times <- seq(0, 200, by = 1)
    out1 <- data.frame(rk(y = init, func = SIRS2, times = times, parms = parms, method = "rk4"))
    
    t <- inc_cont[,1][which(inc_cont[,5] < 1 & inc_cont[,5] > 0.5)]
    t1 <- split(t, cumsum(c(1, diff(t) != 0.5)))[[1]]
    timeequal <<- t1[length(t1)]
    
    
    #Actually Run the Model 
    inc_cont <<- data.frame(0, beta(parms["R0"])*init["I"]*init["S"]*parms["N"] , 0, 0, 0)
    times <- seq(0, 200, by = 1)
    out1 <<- data.frame(rk(y = init, func = SIRS3, times = times, parms = parms, method = "rk4"))
    
    #inc_cont[1,1] <<- 0; inc_cont[1,2] <<- 0
    colnames(inc_cont) <<- c("time", "incidence", "contacttraccap", "beta", "con/inc rat")
    
    timesince <- inc_cont[,1][which(inc_cont[,2] < target)][1]
    return(c("Time1000" = timeequal, "RelTimetoTarget" = timesince))
  }
  contact_trac_runSIS <- function(efficacy, ramp_rate, R0, target, iinput, rinput) {
    
    #Parameters
    N <- 5.5*10^6
    init <- c(S = (N-input$iinput)/N, I = input$iinput/N)
    
    parms = c(gamma = 1/(GenTime(3.3,2.8)), 
              N = 5.5*10^6,
              efficacy = efficacy,
              ramp_rate = ramp_rate, 
              R0 = R0,
              target = target)
    
    #Intersect
    inc_cont <<- data.frame(0, beta(parms["R0"])*init["I"]*init["S"]*parms["N"] , 0, 0, 0)
    times <- seq(0, 200, by = 1)
    out1 <- data.frame(rk(y = init, func = SIRS2SIS, times = times, parms = parms, method = "rk4"))
    
    t <- inc_cont[,1][which(inc_cont[,5] < 1 & inc_cont[,5] > 0.5)]
    t1 <- split(t, cumsum(c(1, diff(t) != 0.5)))[[1]]
    timeequal <<- t1[length(t1)]
    
    
    #Actually Run the Model 
    inc_cont <<- data.frame(0, beta(parms["R0"])*init["I"]*init["S"]*parms["N"] , 0, 0, 0)
    times <- seq(0, 200, by = 1)
    out1 <<- data.frame(rk(y = init, func = SIRS3SIS, times = times, parms = parms, method = "rk4"))
    
    #inc_cont[1,1] <<- 0; inc_cont[1,2] <<- 0
    colnames(inc_cont) <<- c("time", "incidence", "contacttraccap", "beta", "con/inc rat")
    
    timesince <- inc_cont[,1][which(inc_cont[,2] < target)][1]
    return(c("Time1000" = timeequal, "RelTimetoTarget" = timesince))
  }
  
  observe({
    toggleState("rinput", input$model == "SIR")
  })
  
  model <- reactive({
    
    N <- 5.5*10^6
    
    if(input$model == "SIS") {
      
      tracdet <- contact_trac_run(input$efficacy,input$cont_ramp, 
                                  input$R0, input$target, input$iinput, input$rinput)
      
    } 
    
    if(input$model == "SIR") {
      tracdet <- contact_trac_runSIS(input$efficacy,input$cont_ramp, 
                                  input$R0, input$target, input$iinput, input$rinput)
      
    }
    
    list(incidence = inc_cont, infect = out1, details = tracdet)
  })
  
  modelbase <- reactive({
    
    times <- seq(0, 200, by = 1)
    N <- 5.5*10^6
    
    if(input$model == "SIS") {
      
      init <- c(S = (N-input$iinput)/N, I = input$iinput/N)
      
      parms = c(gamma = 1/(GenTime(3.3,2.8)), 
                N = 5.5*10^6,
                R0 = input$R0)
      
      inc_cont1 <<- data.frame(0, beta(parms["R0"])*init["I"]*init["S"]*parms["N"])
      out <<- data.frame(rk(y = init, func = SIRS1SIS, times = times, parms = parms, method = "rk4"))
      colnames(inc_cont1) <- c("time", "incidence")
    } 
    
    if(input$model == "SIR") {
      init <- c(S = (N-(input$iinput + input$rinput))/N, I = input$iinput/N, R = input$rinput/N)
      parms = c(gamma = 1/(GenTime(3.3,2.8)), 
                N = 5.5*10^6,
                R0 = input$R0)
      
      inc_cont1 <<- data.frame(0, beta(parms["R0"])*init["I"]*init["S"]*parms["N"])
      out <<- data.frame(rk(y = init, func = SIRS1, times = times, parms = parms, method = "rk4"))
      colnames(inc_cont1) <- c("time", "incidence")
    }
    
    list(incidence = inc_cont1, infect = out1)
    
  })
  
  
  
  output$Plot1 <- renderPlot({
    
    statsinc <- melt(model()$incidence, id.vars =  c("time"), measure.vars = c("incidence", "contacttraccap"))
    
    ggplot(data = statsinc, aes(x = (time), y = value, col = variable))+ theme_bw() + 
      labs(x ="Time (Days since Contact Tracing Initiation)", y = "Daily Incidence", color = "Population") + scale_y_continuous(limits = c(0,1000), expand = c(0,0)) + 
      scale_x_continuous(limits= c(model()$details["TimeContactTracEqualInc"], model()$details["TimeTargetReach"]), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
            axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
      geom_line(size = 1.02, stat = "identity") + 
      geom_line(data = modelbase()$inc_cont1, aes(x = modelbase()$incidence[1], y = modelbase()$incidence[2]))
    
  })
  
  output$Plot2 <- renderPlot({
    N <- 5.5*10^6
    
    ggplot(data = model()$infect, aes(x = (time), y = I*N))+ theme_bw() +
      labs(x ="Time (Days since Contact Tracing Initiation)", y = "Total Infected", color = "Population") + scale_y_continuous(limits = c(0, 20000),expand = c(0,0)) + 
      scale_x_continuous(limits = c(model()$details["TimeContactTracEqualInc"], model()$details["TimeTargetReach"]), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
            axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
      geom_line(size = 1.02, stat = "identity")
    
  })
  
    
  output$progressBoxdeath <- renderValueBox({
    
    valueBox(tags$p(paste0("Time to Target: ", model()$details[2]), style = "font-size: 50%;"),
             paste0("From when Ccntact tracing Capacity = Incidence"),
             color = "red")
  })
  
  
  output$instructions <- renderUI({
    HTML(paste0("<strong> PLACEHOLDER DESCRIPTION.</strong>",
               sep ="<br/>"))
  })
  
}
  
