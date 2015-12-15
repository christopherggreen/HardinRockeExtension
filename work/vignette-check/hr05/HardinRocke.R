### R code from vignette source '/homes/cggreen/r/lib/HardinRockeExtension/doc/HardinRocke.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: init
###################################################



###################################################
### code chunk number 2: setup (eval = FALSE)
###################################################
# Replicate Table 1 in Hardin and Rocke (2005),
# page 942
#
# Christopher G. Green
# 2014-02-24
#
# run simulations in parallel

require( RhpcBLASctl             )
require( parallel )
require( CerioliOutlierDetection )
require( HardinRockeExtension )

# force single-threaded BLAS if possible
omp_set_num_threads(1)

###################################################
### code chunk number 3: makecluster (eval = FALSE)
###################################################
#thecluster <- makePSOCKcluster(10)

thecluster <- makePSOCKcluster(
	c(
		rep("student1",1),
		rep("student2",1),
		rep("student3",1),
		rep("student4",1),
		rep("student5",1),
		rep("student6",1),
		rep("student7",1),
		rep("student8",1),
		rep("student9",1),
		rep("student10",1),
		rep("student11",0),
		rep("student12",0),
		rep("student13",0),
		rep("student14",0),
		rep("student15",0),
		rep("student16",0),
		rep("laplace1",0),
		rep("laplace2",0),
		rep("laplace3",0),
		rep("laplace4",0),
		rep("laplace5",0),
		rep("laplace6",0),
		rep("laplace7",0),
		rep("laplace8",0),
		rep("laplace9",0),
		rep("laplace10",0),
		rep("laplace11",0),
		rep("laplace12",0),
		rep("laplace13",0),
		rep("laplace14",0),
		rep("laplace15",0),
		rep("laplace16",0),
		rep("newton1",0),
		rep("newton2",0)
	),
	outfile="~/hr05/hr05_replicate.log",
	useXDR=FALSE # these are all little endian machines
)

# make reproducible
clusterSetRNGStream(cl = thecluster, 2015)

###################################################
### code chunk number 4: setsimsize (eval = FALSE)
###################################################
N.SIM <- 5000
B.SIM <- 500


###################################################
### code chunk number 5: clusterinit (eval = FALSE)
###################################################
# initialize each node
tmp.rv <- clusterEvalQ( cl = thecluster, {
  require( CerioliOutlierDetection )
  require( HardinRockeExtension )
  require( mvtnorm )

  my.pid <- Sys.getpid()
  cat("My pid is ", my.pid, "\n")
  my.nodename <- gsub("\\.stat\\.washington\\.edu","",Sys.info()["nodename"])
  cat("My pid is ", my.pid, " on node ", my.nodename, "\n")
  logfile <- paste("~/hr05/Test_Old_Method_Parallel_logfile_",my.nodename,"_",my.pid,".txt",sep="")
  Sys.sleep(30)
  cat("Initialized\n\n", file=logfile)
  cat("My pid is ", my.pid, "on node ", my.nodename, "\n",file=logfile,append=TRUE)

  invisible(NULL)
})


###################################################
### code chunk number 6: buildcases (eval = FALSE)
###################################################
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

cat("Params used:\n", file="hrcases.txt")
sink( "hrcases.txt", append=TRUE)
print( hr.cases )
cat("ncol of params:", ncol(hr.cases),"\n", file="hrcases.txt", append=TRUE)
cat("N.SIM:",N.SIM,"\n", file="hrcases.txt", append=TRUE)
sink( )

###################################################
### code chunk number 7: runsim (eval = FALSE)
###################################################
cat("Starting run at ", format(Sys.time()), "\n")

hrResults <- lapply(hr.cases, function(pn,clst,ns,bs) {
    cat("Trial p = ",pn[1]," n = ",pn[2],"\n")
    hrSimParallel(cl=clst, p = pn[1] , n = pn[2], 
      mcd.alpha=pn[3], alpha=0.05, N=ns, B=bs, lgf=logfile)
  }, clst=thecluster, ns=N.SIM, bs=B.SIM
)

cat("Run completed at ", format(Sys.time()), "\n")
stopCluster(thecluster)

save(hrResults, file="hrResults.20151213.rda")
save(hr.cases, file="hr.cases.20151213.rda")

###################################################
### code chunk number 8: calcstats (eval = FALSE)
###################################################

allmeans <- as.data.frame(t(rbind( hr.cases, 
  100*sapply( hrResults, function(x) colMeans(x) )
)))
row.names(allmeans) <- NULL
allmeans$n <- factor(sprintf("%04d", allmeans$n), ordered=TRUE)
allmeans$mcd.alpha <- sprintf("%0.3f", allmeans$mcd.alpha)
allmeans$mcd.alpha[ allmeans$mbp==1 ] <- "MBP"
allmeans$mcd.alpha <- factor(allmeans$mcd.alpha)

allstds <- as.data.frame(t(rbind( hr.cases, 
  100*sapply( hrResults, function(x) apply(x,2,sd) )
)))
row.names(allstds) <- NULL
allstds$n <- factor(sprintf("%04d", allstds$n), ordered=TRUE)
allstds$mcd.alpha <- sprintf("%0.3f", allstds$mcd.alpha)
allstds$mcd.alpha[ allstds$mbp==1 ] <- "MBP"
allstds$mcd.alpha <- factor(allstds$mcd.alpha)

allmeans.mbp <- subset( allmeans, mbp==1 )
allstds.mbp  <- subset( allstds , mbp==1 )

###################################################
### code chunk number 9: tabularize (eval = FALSE)
###################################################
# format for easier comparision to hardin and rocke paper
# just the mbp cases first

cat("Just MBP values of mcd.alpha.\n")
cat("Means:\n")
print(reshape(allmeans.mbp[,c("p","n","CHI2.CON")  ],
  direction="wide", idvar=c("p"),timevar="n")[,c(1,5:2)])
print(reshape(allmeans.mbp[,c("p","n","HRASY.CON") ], 
  direction="wide", idvar=c("p"),timevar="n")[,c(1,5:2)])
print(reshape(allmeans.mbp[,c("p","n","HRPRED.CON")], 
  direction="wide", idvar=c("p"),timevar="n")[,c(1,5:2)])

cat("Standard deviations:\n")
print(reshape(allstds.mbp[,c("p","n","CHI2.CON")   ],
  direction="wide", idvar=c("p"),timevar="n")[,c(1,5:2)])
print(reshape(allstds.mbp[,c("p","n","HRASY.CON")  ],
  direction="wide", idvar=c("p"),timevar="n")[,c(1,5:2)])
print(reshape(allstds.mbp[,c("p","n","HRPRED.CON") ], 
  direction="wide", idvar=c("p"),timevar="n")[,c(1,5:2)])

# all values of mcd.alpha
cat("All values of mcd.alpha.\n")
cat("Means:\n")
print(reshape(allmeans[,c("p","n","mcd.alpha","CHI2.CON")  ],
  direction="wide", idvar=c("mcd.alpha","p"),timevar="n")[,c(1:2,6:3)])
print(reshape(allmeans[,c("p","n","mcd.alpha","HRASY.CON") ], 
  direction="wide", idvar=c("mcd.alpha","p"),timevar="n")[,c(1:2,6:3)])
print(reshape(allmeans[,c("p","n","mcd.alpha","HRPRED.CON")], 
  direction="wide", idvar=c("mcd.alpha","p"),timevar="n")[,c(1:2,6:3)])

cat("Standard deviations:\n")
print(reshape(allstds[,c("p","n","mcd.alpha","CHI2.CON")   ],
  direction="wide", idvar=c("mcd.alpha","p"),timevar="n")[,c(1:2,6:3)])
print(reshape(allstds[,c("p","n","mcd.alpha","HRASY.CON")  ],
  direction="wide", idvar=c("mcd.alpha","p"),timevar="n")[,c(1:2,6:3)])
print(reshape(allstds[,c("p","n","mcd.alpha","HRPRED.CON") ], 
  direction="wide", idvar=c("mcd.alpha","p"),timevar="n")[,c(1:2,6:3)])


print(reshape(allmeans.mbp[,c("p","n","CHI2.RAW")  ],
  direction="wide", idvar=c("p"),timevar="n")[,c(1,5:2)])
print(reshape(allmeans.mbp[,c("p","n","HRASY.RAW") ], 
  direction="wide", idvar=c("p"),timevar="n")[,c(1,5:2)])
print(reshape(allmeans.mbp[,c("p","n","HRPRED.RAW")], 
  direction="wide", idvar=c("p"),timevar="n")[,c(1,5:2)])


print(reshape(allmeans.mbp[,c("p","n","CHI2.SM")  ],
  direction="wide", idvar=c("p"),timevar="n")[,c(1,5:2)])
print(reshape(allmeans.mbp[,c("p","n","HRASY.SM") ], 
  direction="wide", idvar=c("p"),timevar="n")[,c(1,5:2)])
print(reshape(allmeans.mbp[,c("p","n","HRPRED.SM")], 
  direction="wide", idvar=c("p"),timevar="n")[,c(1,5:2)])


