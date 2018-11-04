install.packages("spls")
install.packages("foreach")
install.packages("doParallel")

library(spls)
library(plyr)
library(foreach)
library(doParallel)
source("spls_functions.R")

# read files 
pw.exp <- read.delim("wax_exp_input.txt", header=T, stringsAsFactors=F)
pw<-as.data.frame(t(pw.exp[,-1]))
tf.exp <- read.delim("maize.TF.exp324_input.txt", header=T, stringsAsFactors=F)
tf<-as.data.frame(t(tf.exp[,-1]))

################################################################################
#Dimension Reduction Step
################################################################################
#obtain spls coefficients
beta <- spls.coeff(tf, pw)
write.table(beta, "Output_spls_beta.txt", row.names=FALSE, sep="\t")

# calculations for t-like statisic for spls coefficients  
boot.tstat <- coeff.boot(tf,pw)
write.table(boot.tstat, "bootstrap_spls_output.txt", row.names=FALSE, sep="\t")



