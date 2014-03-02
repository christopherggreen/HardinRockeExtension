# Monte Carlo simulation to compute the
# Wishart degrees of freedom m associated
# with the Hardin-Rocke approximation to 
# the MCD covariance estimate
#
# See Hardin and Rocke (2005) for the 
# original simulation in the maximum
# breakdown point case. 
# 
# This script extends their work to 
# sampling fractions greater than the
# maximum breakdown point case 
# floor( ( n+v+1 )/2 )/n.
#
# Results of the simulation are documented
# in the technical report " ", 
# Christopher G. Green and R. Douglas Martin
#
#
# Christopher G. Green
# 2014-02-24
#

# run simulations in parallel

library(parallel)

thecluster <- makePSOCKcluster(4)
# make this reproducible
clusterSetRNGStream(cl = thecluster, 2014)

# initialize each node
tmp.rv <- clusterEvalQ( cl = thecluster, {
  require( CerioliOutlierDetection )
  require( HardinRockeExtension )
  N.SIM <- 10#5000
  B.SIM <- 10#250
 
  my.pid <- Sys.getpid()
  cat("My pid is ", my.pid, "\n")
  logfile <- paste("Simulation_AllAlphas_Parallel_logfile_",my.pid,".txt",sep="")
  cat("Initialized\n\n", file=logfile)

  invisible(NULL)
})

# build the pairs of sample size n and dimension p
hr.cm.params <- expand.grid(
	    list(
		  p=c(3,5,7,10,15,20),
		  n=c(50,100,250,500,750,1000)
		)
	  )
# adding more coverage for small sample sizes
hr.cm.params <- rbind( hr.cm.params, within( 
  expand.grid(list(p=c(3,5,7,10,15,20), ratio=c( 3,5,7,9,11 ) )), 
  {
    n <- p * ratio
    rm(ratio)
  }
))
# remove any duplicates
hr.cm.params <- unique(hr.cm.params)
# want to run most expensive cases first
hr.cm.params <- hr.cm.params[ order( hr.cm.params$n, hr.cm.params$p, decreasing=TRUE ), ]

# add maximum breakdown point case to the params data set
hr.cm.params[,"mbp"] <- apply( hr.cm.params, 1, function(x) floor( (x[2] + x[1] + 1)/2 )/x[2] )

# want each case to be a column so that we can use parLapply
hr.cm.params <- data.frame(t(as.matrix(hr.cm.params)))
      
mcd.alphas <- c(0.55,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.99,0.995) 
clusterExport(cl = thecluster, "hr.cm.params")
clusterExport(cl = thecluster, "mcd.alphas")

#
# using parLapply here to prevent simplification of the
# results (as parApply would attempt to do)
#
cat("Starting run at ", format(Sys.time()), "\n")

hr.cm.results.all.pre <- parLapply(cl = thecluster, 
  X = hr.cm.params[1:5], function(pn) {
    cat("Starting case p = ",pn[1]," n = ",pn[2]," at time ", 
	  format(Sys.time()), " \n",file=logfile,append=TRUE)
    results <- hr.cm(p = pn[1] , n = pn[2], N=N.SIM, B=B.SIM, 
	  mcd.alpha=unique(c(pn[3],mcd.alphas)), logfile=logfile)
    cat("Finished case p = ",pn[1]," n = ",pn[2]," at time ", 
	  format(Sys.time()), " \n",file=logfile,append=TRUE)
	data.frame(p=pn[1],n=pn[2],mbp=pn[3],results)
  }
)
cat("Run completed at ", format(Sys.time()), "\n")

stopCluster(thecluster)

hr.cm.results.all.pre[[1]]

hr.cm.results.all <- do.call("rbind", hr.cm.results.all.pre )

save("hr.cm.results.all", file="hr.cm.results.all.final.rda")
save.image()
