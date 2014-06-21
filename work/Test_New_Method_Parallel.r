# Replicate Table 1 in Hardin and Rocke (2005),
# page 942
#
# Christopher G. Green
# 2014-02-24
#

# run simulations in parallel

require( CerioliOutlierDetection )
require( HardinRockeExtension )
library(parallel)

thecluster <- makePSOCKcluster(4)

N.SIM <- 530#0
B.SIM <- 50#0

# initialize each node
tmp.rv <- clusterEvalQ( cl = thecluster, {
  require( CerioliOutlierDetection )
  require( HardinRockeExtension )

  my.pid <- Sys.getpid()
  cat("My pid is ", my.pid, "\n")
  logfile <- paste("Test_New_Method_Parallel_logfile_",my.pid,".txt",sep="")
  cat("Initialized\n\n", file=logfile)

  invisible(NULL)
})

# want each case to be a column so that we can use parLapply
hr.cases <- data.frame(t(as.matrix(expand.grid(list(p=c(5,10,20),n=c(50,100,500,1000),mcd.alpha=c(0.60,0.75,0.95))))))
#clusterExport(cl = thecluster, "hr.cases")

#
# using parLapply here to prevent simplification of the
# results (as parApply would attempt to do)
#
system.time(hrResults <- lapply(hr.cases, function(pn,clst,ns,bs) {
    cat("Trial p = ",pn[1]," n = ",pn[2],"\n")
    hrSimNewParallel(cl=clst, p = pn[1] , n = pn[2], 
	mcd.alpha=pn[3], N=ns, B=bs, lgf=logfile)
  }, clst=thecluster, ns=N.SIM, bs=B.SIM
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
reshape(allmeans[,c("p","n","mcd.alpha","CHI2.CON")  ], direction="wide", idvar=c("mcd.alpha","p"),timevar="n")
reshape(allmeans[,c("p","n","mcd.alpha","CGASY.CON") ], direction="wide", idvar=c("mcd.alpha","p"),timevar="n")
reshape(allmeans[,c("p","n","mcd.alpha","CGPRED.CON")], direction="wide", idvar=c("mcd.alpha","p"),timevar="n")

reshape(allstds[,c("p","n","mcd.alpha","CHI2.CON")   ], direction="wide", idvar=c("mcd.alpha","p"),timevar="n")
reshape(allstds[,c("p","n","mcd.alpha","CGASY.CON")  ], direction="wide", idvar=c("mcd.alpha","p"),timevar="n")
reshape(allstds[,c("p","n","mcd.alpha","CGPRED.CON") ], direction="wide", idvar=c("mcd.alpha","p"),timevar="n")
}

save.image()

