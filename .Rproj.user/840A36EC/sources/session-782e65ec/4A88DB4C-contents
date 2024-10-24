# This example replicates the setup of the shinyapp https://shiny.berryconsultants.com/peter/onearm/

# Load packages
library(ggplot2)
source('R/gjh_Functions.R', local = TRUE)

# Define the values of the sliders on the left sidebar
input <- list()
input$p0 <- 0.10 # HTML("p<sub>0</sub>: Ineffective response rate"), min = 0.00, max = 1.00,
input$p1 <- 0.20 # HTML("p<sub>1</sub>: Desired response rate"), min = 0.00, max = 1.00,
input$alpha <- 0.10 # HTML("&alpha;: Efficacy upper bound on the area below p<sub>0</sub>"), min = 0.00, max = 1.00,
input$beta <- 0.50 # HTML("&beta;: Efficacy upper bound on the area below p<sub>1</sub>"), min = 0.00, max = 1.00,
input$gamma <- 0.50 # HTML("&gamma;: Futility lower bound on the area below p<sub>0</sub>"), min = 0.00, max = 1.00,

# Define the values of the sliders on the Results tab
# Parameters for plot "Varying the Sample Size"
input$p2 <- 0.30 # HTML("True response rate (p):"), min = 0.00, max = 1.00,

# Parameters for plot "Varying the True Response Rate"
input$N <- 24 # HTML("Sample size (N):"), min = 1, max = 60,

# Check value restrictions
# p1 must be >= p0 (see also assignment of glb.p1 and glb.Delta in each plot); move p1 together with p0
if ( input$p1 < input$p0 ) { input$p1 <- input$p0 }

# beta must be >= alpha (see also assignment of glb.beta in each plot); move beta together with alpha
if ( input$beta < input$alpha ) { input$beta <- input$alpha }

# gamma must be >= alpha (see also assignment of glb.gamma in each plot); move gamma together with alpha
if ( input$gamma < input$alpha ) { input$gamma <- input$alpha }

# Create the two plots on the Results tab
output <- list()

# Plot "Varying the True Response Rate"
###################################################################################################
# Set Parameters
###################################################################################################
#
# The meaning of the parameters is explained in the technical report "Inference for Binary Dose 
# Finding Studies Based on Confidence Distributions", version 0.1, by Guenter Heimann
#
glb.p0      <- input$p0
# Ineffective response rate (reference value).
glb.Delta    <- max( input$p1 , input$p0 ) - input$p0 # p1 must be >= p0
# Clinically relevant effect needed to define the trial efficacy criterion.
glb.alpha   <- input$alpha
# Efficacy upper bound on the area below glb.p0; maximum area under the confidence density smaller 
# than glb.Delta=0 (= no treatment effect); this is needed for the trial efficacy criterion.
glb.beta    <- max( input$beta , input$alpha ) # beta must be >= alpha
# Efficacy upper bound on the area below glb.p1; maximum area under the confidence density SMALLER 
# than the clinically relevant effect; this is needed for the success criterion.
glb.epsilon <- 0.8
# desired minimum probability of success.
glb.var_option <- "boundary" # "fixed" or "boundary" or "plus" for the estimator variance
# Choose how to modify the variance, see glb.var_parameter below
glb.var_parameter <- 0 # 1/4 for "fixed", 0 or 1 or 10 for "boundary, 1 or 2 for "plus"
# Choose the parameter for the var_option
# - for "fixed" set this to a non-negative real number, so the estimator variance will be glb.var_parameter / n in all cases; e.g. 1/4 is the maximum possible Bernoulli variance, 0.0 is the minimum possible variance
# - for "boundary", set this to a non-negative integer value, so the estimator variance will be lower bounded by the variance for glb.var_parameter successes; e.g. 1 means only replacing the variance when it is zero (no successes or no failures), 10 is the recommended minimum number of successes and failures for the normal approximation to be acceptable, 0 replaces nothing (allows for 0 variance)
# - for "plus", set this to a non-negative integer value, so the estimator variance will be replaced by the variance of data with glb.var_parameter additional successes and glb.var_parameter additional failures if the number of successes or the number of failures is smaller than 10 (which is required for the normal approximation to be acceptable); e.g. 1 adds one success and one failure, 2 adds two successes and two failures (known as the "plus four" rule), 0 replaces nothing (allows for 0 variance)
#
###################################################################################################


###################################################################################################
# "ByParameter" requires to fix the sample size (n), it will vary the success probability under treatment (p) parameter
#
glb.gamma   <- max( input$gamma, input$alpha ) # gamma must be >= alpha
# Futility lower bound on the area below glb.p0; maximum area under the confidence density larger 
# than glb.Delta=0 (= no treatment effect); this is needed for the failure criterion
glb.n       <- input$N
# Sample size of the treatment group
glb.step    <- 0.01
# defines the possible values for the success probability under treatment; basically all values
# starting with glb.p0 going up to 1 in steps of glb.step will be explored.
glb.xstep   <- 0.1 #0.04
# defines the steps on the x-axis for plotting
glb.primaryendpoint   <- c("True response rate (p)")
# to be used for the x-axis label
#
###################################################################################################

###################################################################################################
# Calculate titles, file names, numbers, etc which may depend upon parameters
###################################################################################################
#
glb.p1    <- max( input$p1 , input$p0 ) # p1 must be >= p0
glb.NumberOfScenario <- paste("One-arm_p0=",glb.p0, "_p1=", glb.p1, "_N=", glb.n, "_a=", glb.alpha, "_b=", glb.beta,
                              "_g=", glb.gamma, "_e=", glb.epsilon, "_var=", glb.var_option, glb.var_parameter, sep="")
glb.NameOfScenario   <- paste("One-arm: p0=",  glb.p0, 
                              " p1=",    glb.p1, 
                              " N=",        glb.n, 
                              #" alpha=",    glb.alpha,
                              #" beta=",     glb.beta,
                              " gamma=",    glb.gamma,
                              #" epsilon=",  glb.epsilon,
                              " Var=", glb.var_option, glb.var_parameter, 
                              sep="") 
glb.filename       <- paste(glb.NumberOfScenario, ".pdf", sep ="")
glb.globaltitle    <- ""#paste("OC for ", glb.NameOfScenario, sep="") 
#      glb.outputpath     <- file.path(gjh.DirectoryOutputs, glb.filename) 
glb.xaxis          <- glb.primaryendpoint#paste(glb.primaryendpoint, "under treatment") 
# 
###################################################################################################

###################################################################################################
# Calculate Success and Failure probabilities
###################################################################################################
#
# The vectors calculated below are obtained following the formulae (13) and (17) as described in 
# the technical report "Inference for Binary Dose Finding Studies Based on Confidence Distributions", 
# version 0.1, by Guenter Heimann
glb.successprobability<-glf.successbyparameter.one(glb.alpha,glb.beta,glb.Delta,glb.p0,glb.step,glb.n,glb.var_option,glb.var_parameter)
glb.failureprobability<-glf.failurebyparameter.one(glb.gamma,glb.p0,glb.step,glb.n,glb.var_option,glb.var_parameter)
glb.neitherprobability<-1-glb.successprobability-glb.failureprobability
glb.probability <- seq(glb.p0,1, by=glb.step)
glb.alphavector <- rep(glb.alpha,times=length(glb.probability))
glb.betavector  <- rep(glb.beta,times=length(glb.probability))
glb.epsilonvector <- rep(glb.epsilon,times=length(glb.probability))
glb.results <- t(rbind(glb.probability, glb.successprobability,glb.failureprobability,glb.neitherprobability,
                       glb.alphavector, glb.betavector,glb.epsilonvector))
#
###################################################################################################

###################################################################################################
# Prepare Results for Plotting
###################################################################################################
#
text1<-paste("\u03B1=", glb.alpha, sep="")
text2<-paste("\u03B2=", glb.beta, sep="")
text3<-paste("\u03B5=", glb.epsilon, sep="")
#
glb.results <- as.data.frame(glb.results)
rownames(glb.results)=NULL
glb.DataToPlot <- rbind(transform(glb.results, variable='Probability of Efficacy',     value=glb.successprobability),
                        transform(glb.results, variable='Probability of Futility',     value=glb.failureprobability),
                        transform(glb.results, variable='Probability of No Decision', value=glb.neitherprobability),
                        transform(glb.results, variable=text1,                  value=glb.alphavector),
                        transform(glb.results, variable=text2,                  value=glb.betavector),
                        transform(glb.results, variable=text3,                  value=glb.epsilonvector))
#
#glb.DataToPlot$value=as.numeric(levels(glb.DataToPlot$value))[glb.DataToPlot$value] # ORIGINAL
glb.DataToPlot$value=as.numeric(glb.DataToPlot$value)
#glb.DataToPlot$glb.probability=as.numeric(levels(glb.DataToPlot$glb.probability))[glb.DataToPlot$glb.probability] # ORIGINAL
glb.DataToPlot$glb.probability=as.numeric(glb.DataToPlot$glb.probability)
# somehow these two variables were factors in the data frame, and they need to be converted to numeric for 
# plotting purposes.
glb.DataToPlot <- subset(glb.DataToPlot,select=c(glb.probability, variable, value))
glb.DataToPlot <- as.data.frame(glb.DataToPlot)
# Reorder the legend
glb.DataToPlot$variable <- factor(glb.DataToPlot$variable, levels = c('Probability of Efficacy', 'Probability of Futility', 'Probability of No Decision', text1 , text2 , text3 ))
#
###################################################################################################

###################################################################################################
# Plot the results 
###################################################################################################
#
#
glb.palette     <- c('green4','red', 'darkorange', 'turquoise4', 'yellow4', 'green4')
glb.linestyles  <- c('solid','solid','solid','longdash','longdash', 'dashed')
glb.linescales  <- c(0.75,0.75,0.75,0.5,0.5,0.5)
theme_set(theme_bw())

my_plot <- ggplot(glb.DataToPlot,aes(glb.probability,value,colour=variable,linetype=variable,linewidth=variable)) + 
  geom_line() + 
  geom_vline(xintercept=glb.p1, colour='gray') +
  scale_colour_manual(values=glb.palette, name=" ") +
  scale_linetype_manual(values=glb.linestyles, name=" ") +
  scale_linewidth_manual(values=glb.linescales, name=" ") +
  ylab('Probability') + 
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,by=.1)) +
  xlab(glb.xaxis) + 
  scale_x_continuous(breaks=seq(glb.p0,1,by=glb.xstep)) +
  theme(legend.position = "inside", legend.position.inside=c(.8,.35)) +
  theme(legend.title=element_blank())+
  labs(title=glb.globaltitle)

# print it in the R Studio
my_plot
#
###################################################################################################
output$plot_by_probability <- my_plot



# Plot "Varying the Sample Size"
###################################################################################################
# Set Parameters
###################################################################################################
#
# The meaning of the parameters is explained in the technical report "Inference for Binary Dose 
# Finding Studies Based on Confidence Distributions", version 0.1, by Guenter Heimann
#
glb.p0      <- input$p0
# Ineffective response rate (reference value).
glb.Delta    <- max( input$p1 , input$p0 ) - input$p0 # p1 must be >= p0
# Clinically relevant effect needed to define the trial efficacy criterion.
glb.alpha   <- input$alpha
# Efficacy upper bound on the area below glb.p0; maximum area under the confidence density smaller 
# than glb.Delta=0 (= no treatment effect); this is needed for the trial efficacy criterion.
glb.beta    <- max( input$beta , input$alpha ) # beta must be >= alpha      # Efficacy upper bound on the area below glb.p1; maximum area under the confidence density SMALLER 
# than the clinically relevant effect; this is needed for the success criterion.
glb.epsilon <- 0.8
# desired minimum probability of success.
glb.var_option <- "boundary" # "fixed" or "boundary" or "plus" for the estimator variance
# Choose how to modify the variance, see glb.var_parameter below
glb.var_parameter <- 0 # 1/4 for "fixed", 0 or 1 or 10 for "boundary, 1 or 2 for "plus"
# Choose the parameter for the var_option
# - for "fixed" set this to a non-negative real number, so the estimator variance will be glb.var_parameter / n in all cases; e.g. 1/4 is the maximum possible Bernoulli variance, 0.0 is the minimum possible variance
# - for "boundary", set this to a non-negative integer value, so the estimator variance will be lower bounded by the variance for glb.var_parameter successes; e.g. 1 means only replacing the variance when it is zero (no successes or no failures), 10 is the recommended minimum number of successes and failures for the normal approximation to be acceptable, 0 replaces nothing (allows for 0 variance)
# - for "plus", set this to a non-negative integer value, so the estimator variance will be replaced by the variance of data with glb.var_parameter additional successes and glb.var_parameter additional failures if the number of successes or the number of failures is smaller than 10 (which is required for the normal approximation to be acceptable); e.g. 1 adds one success and one failure, 2 adds two successes and two failures (known as the "plus four" rule), 0 replaces nothing (allows for 0 variance)
#
###################################################################################################


###################################################################################################
# "BySamplesize" requires to fix the success probability under treatment (p2), it will vary the sample size (n) parameter
#
glb.p2      <- input$p2 #glb.p0 + 2 * glb.Delta #0.3#0.4  
# survival probability under treatment used to evaluate the operating characteristics.
glb.p_vec       <- c(glb.p0,glb.p2)
glb.maxss       <- 60#1000  # ONLY NEEDED FOR "BY SAMPLE SIZE"
# maximum sample size allowed
glb.step    <- 1 
# an integer which defines the possible values for the sample size under placebo; basically all 
# values starting with 1 going up to glb.maxss in steps of glb.step will be explored.
glb.xstep   <- 4#50 
# defines the steps on the x-axis for plotting
#
###################################################################################################

###################################################################################################
# Calculate titles, file names, numbers, etc which may depend upon parameters
###################################################################################################
#
glb.p1    <- max( input$p1 , input$p0 ) # p1 must be >= p0
glb.NumberOfScenario <- paste("One-arm_p0=",glb.p0, "_p1=", glb.p1, "_p2=", glb.p2, "_a=", glb.alpha, "_b=", glb.beta,
                              "_e=", glb.epsilon, "_maxss=", glb.maxss, "_var=", glb.var_option, glb.var_parameter, 
                              sep="") 
glb.NameOfScenario   <- paste("One-arm: p0=",  glb.p0, 
                              " p1=",    glb.p1,
                              " p2=",        glb.p2, 
                              #" alpha=",    glb.alpha,
                              #" beta=",     glb.beta,
                              #" epsilon=",  glb.epsilon,
                              " Var=", glb.var_option, glb.var_parameter, 
                              sep="") 
glb.filename       <- paste(glb.NumberOfScenario, ".pdf", sep ="")
#glb.globaltitle    <- paste(glb.NameOfScenario, " Trt sample size=" , sep="")
glb.globaltitle    <- ""#glb.NameOfScenario
#      glb.outputpath     <- file.path(gjh.DirectoryOutputs, glb.filename) 
glb.xaxis          <- paste("Sample size (N)")# under treatment
# 
###################################################################################################

###################################################################################################
# Calculate Success and Failure probabilities
###################################################################################################
#
# The vectors calculated below are obtained following the formulae (13) and (17) as described in 
# the technical report "Inference for Binary Dose Finding Studies Based on Confidence Distributions", 
# version 0.1, by Guenter Heimann
glb.successprobability<-glf.successbysamplesize.one(glb.alpha,glb.beta,glb.Delta,glb.p_vec,glb.maxss,glb.step,glb.var_option,glb.var_parameter)
glb.samplesize <- seq(1,glb.maxss, by=glb.step)
glb.alphavector <- rep(glb.alpha,times=length(glb.samplesize))
glb.betavector  <- rep(glb.beta,times=length(glb.samplesize))
glb.epsilonvector <- rep(glb.epsilon,times=length(glb.samplesize))
glb.condition  <- ifelse(glb.successprobability<=glb.epsilonvector,glb.samplesize,0)
glb.results <- t(rbind(glb.samplesize, glb.successprobability,glb.alphavector,
                       glb.betavector,glb.epsilonvector))
glb.treatmentsamplesize <- (max(glb.condition)+1)
#
###################################################################################################

###################################################################################################
# Prepare Results for Plotting
###################################################################################################
#
text1<-paste("\u03B1=", glb.alpha, sep="")
text2<-paste("\u03B2=", glb.beta, sep="")
text3<-paste("\u03B5=", glb.epsilon, sep="")

glb.results <- as.data.frame(glb.results)
rownames(glb.results)=NULL
glb.DataToPlot <- rbind(transform(glb.results, variable='Probability of Efficacy', value=glb.successprobability),
                        transform(glb.results, variable=text1,              value=glb.alphavector),
                        transform(glb.results, variable=text2,              value=glb.betavector),
                        transform(glb.results, variable=text3,              value=glb.epsilonvector))
#
#glb.DataToPlot$value=as.numeric(levels(glb.DataToPlot$value))[glb.DataToPlot$value] # ORIGINAL
glb.DataToPlot$value=as.numeric(glb.DataToPlot$value)
#glb.DataToPlot$glb.samplesize=as.numeric(levels(glb.DataToPlot$glb.samplesize))[glb.DataToPlot$glb.samplesize] # ORIGINAL
glb.DataToPlot$glb.samplesize=as.numeric(glb.DataToPlot$glb.samplesize)
# somehow these two variables were factors in the data frame, and they need to be converted to numeric for 
# plotting purposes.
glb.DataToPlot <- subset(glb.DataToPlot,select=c(glb.samplesize, variable, value))
glb.DataToPlot <- as.data.frame(glb.DataToPlot)
#glb.globaltitle<- paste(glb.globaltitle, glb.treatmentsamplesize, sep="")
# Reorder the legend
glb.DataToPlot$variable <- factor(glb.DataToPlot$variable, levels = c('Probability of Efficacy', text1 , text2 , text3 ))
#
###################################################################################################


###################################################################################################
# Plot the results 
###################################################################################################
#
#
glb.palette     <- c('green4', 'turquoise4', 'yellow4', 'green4')
glb.linestyles  <- c('solid','longdash','longdash', 'dashed')
glb.linescales  <- c(0.75,0.5,0.5,0.5)
theme_set(theme_bw())

my_plot <- ggplot(glb.DataToPlot,aes(glb.samplesize,value,colour=variable,linetype=variable,linewidth=variable)) + 
  geom_line() + 
  scale_colour_manual(values=glb.palette, name=" ") +
  scale_linetype_manual(values=glb.linestyles, name=" ") +
  scale_linewidth_manual(values=glb.linescales, name=" ") +
  ylab('Probability') + 
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,by=.1)) +
  xlab(glb.xaxis) + 
  # xlab('day 28 survival rate under treatment') + 
  scale_x_continuous(breaks=seq(0,glb.maxss,by=glb.xstep)) +
  theme(legend.position = "inside", legend.position.inside=c(.8,.35)) +
  theme(legend.title=element_blank())+
  labs(title=glb.globaltitle)

# print it in the R Studio
my_plot
#
###################################################################################################
output$plot_by_samplesize <- my_plot
