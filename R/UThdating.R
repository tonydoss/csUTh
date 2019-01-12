#' UThdating
#'
#' UThdating calculates closed-system Th-230/U ages, including detrital correction.
#'
#' The data need to be in a tab-separated text file named 'IoliteExport_All_Integrations.txt'.
#'
#' The following columns need to be present in the data file: U234_U238_CORR, U234_U238_CORR_Int2SE, Th230_U238_CORR, Th230_U238_CORR_Int2SE, Th232_U238_CORR, Th232_U238_CORR_Int2SE.
#'
#' @param sample_name_choice Name of the sample to solve. It needs to be exactly as in the data file. For example, you want to solve for sample MK16, enter 'MK16'. Default: 'MK16'
#' @param nbitchoice Number of iterations in the model. Have at least 100. Default: 100.
#' @param detcorrectionchoice Whether to do a detrital correction. Enter 'Y' for yes, or 'N' for no. Default: 'Y'.
#' @param keepfiltereddata Whether to do save filtered data on which an outlier test was performed. Only recommended if all analyses of a same sample are supposed to give the same age. Enter 'Y' for yes, or 'N' for no. Default: 'N'.

UThdating <- function(sample_name_choice = 'MK16',
                      nbitchoice = 100,
                      detcorrectionchoice = 'Y',
                      keepfiltereddata = 'N',
                      print_summary = TRUE
){

  if (!require("deSolve")) install.packages("deSolve")
  if (!require("ggplot2")) install.packages("ggplot2")

  library(deSolve)
  library(ggplot2)

  l234 <- 2.8262e-6 # 234U decay constant (a-1)
  l230 <- 9.1577e-6 # 230Th decay constant (a-1)

  # name of sample to solve
  sample_name <- sample_name_choice

  # nb of times optimisation is repeated (for each sample)
  nbit <- nbitchoice

  lowerbound <- c(2, 01.0) # lower bound values for age (log10(yr)) and initial (234U/238U)
  upperbound <- c(6, 10.0) # upper bound values for age (log10(yr)) and initial (234U/238U)

  # use detrital correction (Y) or not (N)
  detcorrection <- detcorrectionchoice

  # composition detritus
  R28det <- 0.8
  R28det_err <- 0.08
  R08det <- 1
  R08det_err <- 0.05
  R48det <- 1
  R48det_err <- 0.02

   # import iolite results
  iolite_results <- read.table("input/IoliteExport_All_Integrations.txt",
                               header = TRUE, sep = "\t", comment.char = "")

  # create dataframe with data only for samples to solve
  data <- subset(iolite_results, (grepl(sample_name, X)))
  # data <- subset(ABL_AV_C, MATERIAL == sample_name)
  # data <- data[-c(8,11),]
  # number of samples to solve
  number_sampletosolve <- nrow(data)

  # parameters for detrial correction
  data$B <- (R08det - data$Th230_U238_CORR)/(R28det - data$Th232_U238_CORR)
  data$b <- (R48det - data$U234_U238_CORR)/(R28det - data$Th232_U238_CORR)
  data$r1 <- data$Th232_U238_CORR/(R28det - data$Th232_U238_CORR)
  data$r2 <- R28det/(R28det - data$Th232_U238_CORR)

  # detrital-corrected ratios
  data$U234_U238_DET_CORR <- data$U234_U238_CORR - data$b*data$Th232_U238_CORR
  data$U234_U238_DET_CORR_ERR <- sqrt(data$b^2*(data$r2^2*data$Th232_U238_CORR_Int2SE^2 +
                                                  data$r1^2*R28det_err^2) +
                                        data$r2^2*data$U234_U238_CORR_Int2SE +
                                        data$r1^2*R48det_err^2)
  data$Th230_U238_DET_CORR <- data$Th230_U238_CORR - data$B*data$Th232_U238_CORR
  data$Th230_U238_DET_CORR_ERR <- sqrt(data$B^2*(data$r2^2*data$Th232_U238_CORR_Int2SE^2 +
                                                   data$r1^2*R28det_err^2) +
                                         data$r2^2*data$Th230_U238_CORR_Int2SE +
                                         data$r1^2*R08det_err^2)

  # create vectors
  time_results <- vector(mode="numeric", length=number_sampletosolve)
  err_time_results <- vector(mode="numeric", length=number_sampletosolve)
  R48i_results <- vector(mode="numeric", length=number_sampletosolve)
  err_R48i_results <- vector(mode="numeric", length=number_sampletosolve)
  time2se_results <- vector(mode="numeric", length=number_sampletosolve)
  R48i2se_results <- vector(mode="numeric", length=number_sampletosolve)

  # repeat loop for each sample
  for (count in 1:number_sampletosolve){
    if (detcorrection == 'Y'){
      U48meas <- data$U234_U238_DET_CORR[count]
      Th0U8meas <- data$Th230_U238_DET_CORR[count]
      err_R08 <- data$Th230_U238_DET_CORR_ERR[count]
      err_R48 <- data$U234_U238_DET_CORR_ERR[count]
    } else {
      U48meas <- data$U234_U238_CORR[count]
      Th0U8meas <- data$Th230_U238_CORR[count]
      err_R08 <- data$Th230_U238_CORR_Int2SE[count]
      err_R48 <- data$U234_U238_CORR_Int2SE[count]
    }
    # create list for optimisation results
    sol=list()
    # create vectors
    U48calc <- vector(mode="numeric", length=nbit)
    Th0U8calc <- vector(mode="numeric", length=nbit)
    time <- vector(mode="numeric", length=nbit)
    R48i <- vector(mode="numeric", length=nbit)
    time_2se <- vector(mode="numeric", length=nbit)
    R48i_2se <- vector(mode="numeric", length=nbit)

    # repeat optimisation 'nbit' number of times for a given sample
    for (i in 1:nbit){
      # pick (234U/238U) and (230Th/238U) within range of measured values
      U48target <- runif(1, U48meas - err_R48, U48meas + err_R48)
      Th0U8target <- runif(1, Th0U8meas - err_R08, Th0U8meas + err_R08)

      # start optimisation with random age and initial (234U/23U) taken from the range of values allowed
      init_time <- runif(1, lowerbound[1], upperbound[1])
      init_R48i <- runif(1, lowerbound[2], upperbound[2])
      paraminit <- c(init_time, init_R48i)

      # function to minimise
      funmin <- function(x) {
        t <- 10^x[1] # age in yr
        U48i <- x[2] # intial (234U/238U)

        U48calc <- 1 + (U48i - 1)*exp(-l234*t) # (234U/238U)
        Th0U8calc <- 1 - exp(-l230*t) + (U48calc - 1)*(l230/(l230 - l234))*(1 - exp((l234 - l230)*t))

        fmin <- sum((U48calc - U48target)^2 + (Th0U8calc - Th0U8target)^2 ) # function to minimise
      }

      # optimisation
      sol <- optim(paraminit, funmin, method = "L-BFGS-B",
                   lower = lowerbound, upper = upperbound, control = list(factr = 1e-8))
      # store calculated age, initial (234U/23U) and calculated activity ratios for each optimisation
      time[i] <- 10^sol$par[1]
      R48i[i] <- sol$par[2]
      U48calc[i] <- 1 + (R48i[i] - 1)*exp(-l234*time[i]) # (234U/238U)
      Th0U8calc[i] <- 1 - exp(-l230*time[i]) - (U48calc[i] - 1)*(l230/(l234 - l230))*(1 - exp((l234 - l230)*time[i]))
    }

    # store results from all optimisations for a given sample
    results <- as.data.frame(cbind(time, R48i, U48calc, Th0U8calc))
    # take the median of all ages and initial (234U/23U)
    time_median <- median(results$time)
    time_2se <- 2*sd(results$time)/sqrt(nbit)
    R48i_median <- median(results$R48i)
    R48i_2se <- 2*sd(results$R48i)/sqrt(nbit)

    # library(ggplot2)
    # ggplot(data=results, aes(time)) + geom_histogram()

    # calculate error on age and initial (234U/23U)
    # gamma_0 <- 1 + (U48meas - 1)*exp(l234*median_time)
    # k1 <- l230/(l230 - l234)*(1 - exp(-(l230 - l234)*median_time))
    # k2 <- l234*(U48meas - 1)*exp(l234*median_time)
    # D <- l230*(exp(-l230*median_time) + (U48meas-1)*exp(-(l230 - l234)*median_time))
    # err_time <- sqrt((err_R08^2 + k1^2*err_R48^2)/D^2)#/2
    # # err_time <- sqrt((err_R08^2 + k1^2*err_R48^2 + 2*k1*cov(Th0U8meas, U48meas))/D^2)
    # err_R48i <- sqrt(k2*err_time^2 + (exp(l234*median_time)*err_R48)^2 + 2*k2*exp(l234*median_time)/D*
    #                    (k1*err_R48^2))#/1000

    # store age, error on age and initial (234U/23U) for each sample
    time_results[count] <- time_median
    time2se_results[count] <- time_2se
    R48i_results[count] <- R48i_median
    R48i2se_results[count] <- R48i_2se
  }

  final_results <- as.data.frame(cbind(as.character(data$X),
                                       data[,27:36], data[,41:44],
                                       round(time_results/1000,3), round(time2se_results/1000,3),
                                       round(R48i_results,3), round(R48i2se_results,3)))
  colnames(final_results)[1] <- c("ID")
  colnames(final_results)[16:19] <- c("Age (kyr)", "2se", "(234U/238U)i", "2se")

  remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }

  final_results_filtered <- final_results
  final_results_filtered$`Age (kyr)` <- remove_outliers(final_results$`Age (kyr)`)
  final_results_filtered <- final_results_filtered[!is.na(final_results_filtered$`Age (kyr)`),]

  if (keepfiltereddata == 'N'){
    write.table(final_results, file = paste("output/",sample_name,".csv", sep = ""), sep = ",", row.names = F)
    plotdata <- final_results
  } else if (keepfiltereddata == 'Y'){
    write.table(final_results_filtered, file = paste("output/",sample_name,".csv", sep = ""), sep = ",", row.names = F)
    plotdata <- final_results_filtered
  }

  # change column name of initial (234U/238U) error so it can be used to show error bars
  colnames(plotdata)[16:19] <- c("Age (kyr)", "2se", "(234U/238U)i", "2se#2")

  # plot initial (234U/238U)
  p2 <- ggplot(plotdata, aes(ID, `(234U/238U)i`)) + # plot ages
    geom_errorbar(aes(ymin = (`(234U/238U)i` - `2se#2`),ymax = (`(234U/238U)i` + `2se#2`)), width=0.1) + # plot error bars
    geom_point(size=5) + # plot points
    xlab("Sample ID") + # x axis label
    ylab(expression("Initial ("^234*"U/"^238*"U)")) # y axis label

  p2

  # plot ages
  p1 <- ggplot(plotdata, aes(ID, `Age (kyr)`)) + # plot ages
    geom_errorbar(aes(ymin = (`Age (kyr)` - `2se`),ymax = (`Age (kyr)` + `2se`)), width=0.1) + # plot error bars
    geom_point(size=5) + # plot points
    xlab("Sample ID") + # x axis label
    ylab("Age (ka)") # y axis label

  print(p1)

  ggsave(paste("output/",sample_name," - Age.png",sep = ""))

  # display results
  print(plotdata[,16:19])

  print(paste('Mean age: ',round(mean(plotdata$`Age (kyr)`, na.rm = TRUE),1),
              '+/-', round(2*sd(plotdata$`Age (kyr)`, na.rm = TRUE)/
                             sqrt(length(plotdata$`Age (kyr)`)), 1), ' ka'))
}
