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

thecluster <- makePSOCKcluster(10)

N.SIM <- 5000
B.SIM <- 500

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
hr.cases <- as.matrix(expand.grid(list(p=c(5,10,20),n=c(50,100,500,1000),
	mcd.alpha=c(0.65,0.75,0.85,0.95,0.99))))
hr.cases <- rbind(cbind(hr.cases,FALSE),
	t(apply(as.matrix(expand.grid(list(p=c(5,10,20),n=c(50,100,500,1000)))),
		1,function(x) c(x,mcd.alpha=floor( (x[2] + x[1] + 1)/2 )/x[2],mbp=TRUE) )))
hr.cases <- unique(hr.cases)
hr.cases <- hr.cases[ order(hr.cases[,"mcd.alpha"]),]
hr.cases <- hr.cases[ order(hr.cases[,"n"],decreasing=TRUE),]
dimnames(hr.cases)[[2]] <- c("p","n","mcd.alpha","mbp")
hr.cases <- data.frame(t(hr.cases))
#clusterExport(cl = thecluster, "hr.cases")

#
# using parLapply here to prevent simplification of the
# results (as parApply would attempt to do)
#

cat("Starting run at ", format(Sys.time()), "\n")

hrResults <- lapply(hr.cases, function(pn,clst,ns,bs) {
    cat("Trial p = ",pn[1]," n = ",pn[2],"\n")
    hrSimNewParallel(cl=clst, p = pn[1] , n = pn[2], 
	mcd.alpha=pn[3], N=ns, B=bs, lgf=logfile)
  }, clst=thecluster, ns=N.SIM, bs=B.SIM
)

cat("Run completed at ", format(Sys.time()), "\n")

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


allmeans.hrsimnew05 <- allmeans
allstds.hrsimnew05 <- allstds

save("allmeans.hrsimnew05", file="allmeans.hrsimnew05.rda")
save("allstds.hrsimnew05", file="allstds.hrsimnew05.rda")



save.image()

