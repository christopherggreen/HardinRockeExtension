# Replicate Table 1 in Hardin and Rocke (2005),
# page 942
#
# Christopher G. Green
# 2014-02-24
#

# run simulations in parallel

library(parallel)

thecluster <- makePSOCKcluster(4)
# initialize each node
tmp.rv <- clusterEvalQ( cl = thecluster, {
  require( CerioliOutlierDetection )
  require( HardinRockeExtension )
  N.SIM <- 5000
  B.SIM <- 500


  NULL
})

# want each case to be a column so that we can use parLapply
hr.cases <- data.frame(t(as.matrix(expand.grid(list(p=c(5,10,20),n=c(50,100,500,1000))))))
clusterExport(cl = thecluster, "hr.cases")

#
# using parLapply here to prevent simplification of the
# results (as parApply would attempt to do)
#
system.time(hrResults <- parLapply(cl = thecluster, 
  X = hr.cases, function(pn) {
    #cat("Trial p = ",p," n = ",n,"\n")
    hrSim(p = pn[1] , n = pn[2], N=N.SIM, B=B.SIM)
  }
))/60

stopCluster(thecluster)

if ( FALSE ) {
# print column means for each result
lapply( hrResults, function(x) 100*colMeans(x) )

# print column stdevs for each result
lapply( hrResults, function(x) 100*apply(x,2,sd) )

}

allmeans <- as.data.frame(t(rbind( hr.cases, 
	100*sapply( hrResults, function(x) colMeans(x) )
)))
row.names(allmeans) <- NULL

allstds <- as.data.frame(t(rbind( hr.cases, 
	100*sapply( hrResults, function(x) apply(x,2,sd) )
)))
row.names(allstds) <- NULL

if ( FALSE ) {
# format for easier comparision to hardin and rocke paper
reshape(allmeans[,c("p","n","CHI2.CON")  ], direction="wide", idvar="p",timevar="n")
reshape(allmeans[,c("p","n","HRASY.CON") ], direction="wide", idvar="p",timevar="n")
reshape(allmeans[,c("p","n","HRPRED.CON")], direction="wide", idvar="p",timevar="n")

reshape(allstds[,c("p","n","CHI2.CON")   ], direction="wide", idvar="p",timevar="n")
reshape(allstds[,c("p","n","HRASY.CON")  ], direction="wide", idvar="p",timevar="n")
reshape(allstds[,c("p","n","HRPRED.CON") ], direction="wide", idvar="p",timevar="n")
}

save.image()

