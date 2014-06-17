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

thecluster <- makePSOCKcluster(10)
# make this reproducible
clusterSetRNGStream(cl = thecluster, 2014)

# initialize each node
tmp.rv <- clusterEvalQ( cl = thecluster, {
  require( CerioliOutlierDetection )
  require( HardinRockeExtension )
  N.SIM <- 5000
  B.SIM <- 250
 
  my.pid <- Sys.getpid()
  cat("My pid is ", my.pid, "\n")
  logfile <- paste("Simulation_AlphaNear1LargeN_Parallel_logfile_",my.pid,".txt",sep="")
  cat("Initialized\n\n", file=logfile)

  invisible(NULL)
})

# build the pairs of sample size n and dimension p
hr.cmstarbig.params <- expand.grid(
	    list(
		  p=c(3,5,7,10,15,20),
		  n=c(500,750,1000)
		)
	  )

# remove any duplicates
hr.cmstarbig.params <- unique(hr.cmstarbig.params)

# add maximum breakdown point case to the params data set
hr.cmstarbig.params[,"mbp"] <- apply( hr.cmstarbig.params, 1, function(x) floor( (x[2] + x[1] + 1)/2 )/x[2] )

# run each case 10 times
nnn <- nrow(hr.cmstarbig.params)
hr.cmstarbig.params <- apply( hr.cmstarbig.params, 2, function(x) rep(x,10) )
hr.cmstarbig.params <- cbind(hr.cmstarbig.params, "Case"=rep(1:10, each=nnn))

# want to run most expensive cases first
hr.cmstarbig.params <- hr.cmstarbig.params[ order( hr.cmstarbig.params[,"n"], hr.cmstarbig.params[,"p"], decreasing=TRUE ), ]

# want each case to be a column so that we can use parLapply
hr.cmstarbig.params <- data.frame(t(as.matrix(hr.cmstarbig.params)))
names(hr.cmstarbig.params) <- as.character(seq(ncol(hr.cmstarbig.params)))
      
mcd.alphas <- c(0.90,0.95,0.99,0.995) 
clusterExport(cl = thecluster, "hr.cmstarbig.params")
clusterExport(cl = thecluster, "mcd.alphas")

#
# using parLapply here to prevent simplification of the
# results (as parApply would attempt to do)
#

cat("Starting run at ", format(Sys.time()), "\n")

hr.cmstarbig.results.all.pre <- parLapply(cl = thecluster, 
  X = hr.cmstarbig.params, function(pn) {
    cat("Starting case p = ",pn[1]," n = ",pn[2]," at time ", 
	  format(Sys.time()), " \n",file=logfile,append=TRUE)
    results <- hr.cm(p = pn[1] , n = pn[2], N=N.SIM, B=B.SIM, 
	  mcd.alpha=unique(mcd.alphas), logfile=logfile)
    cat("Finished case p = ",pn[1]," n = ",pn[2]," at time ", 
	  format(Sys.time()), " \n",file=logfile,append=TRUE)
	data.frame(p=pn[1],n=pn[2],mbp=pn[3],case=pn[4],results)
  }
)
cat("Run completed at ", format(Sys.time()), "\n")

stopCluster(thecluster)

hr.cmstarbig.results.all.pre[[1]]

hr.cmstarbig.results.all <- do.call("rbind", hr.cmstarbig.results.all.pre )

save("hr.cmstarbig.results.all", file="hr.cmstarbig.results.all.final.rda")
save.image()
