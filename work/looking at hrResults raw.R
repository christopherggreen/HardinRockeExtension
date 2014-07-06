load("hrResults01raw.rda")
load("hrResults05raw.rda")
load("hrcases.rda")
hr.cases

summary(hrResults01[[1]][,"CGPRED.CON"])
boxplot( hrResults01[[1]][,"CGPRED.CON"] )

hrResults01.cgpredcon <- cbind(
	apply(t(hr.cases),2,function(x) rep(x,each=5000)),
	do.call("c",lapply( hrResults01, function(x) x[,"CGPRED.CON",drop=TRUE] ))
)
hrResults01.cgpredcon <- data.frame( hrResults01.cgpredcon, row.names=NULL )
names( hrResults01.cgpredcon ) <- c("p","n","mcd.alpha","mbp","cgpred.con")

hrResults01.cgpredcon$mcd.alpha.bin <- ifelse( hrResults01.cgpredcon$mbp, "MBP", 
  format( hrResults01.cgpredcon$mcd.alpha, 3 ) )
# make ordered factor with MBP first
thesorter <- function(x) { n <- length(x); x[c(n,1:(n-1))] }
hrResults01.cgpredcon$mcd.alpha.bin <- factor( hrResults01.cgpredcon$mcd.alpha.bin, ordered=TRUE,
    levels=thesorter(sort(unique(hrResults01.cgpredcon$mcd.alpha.bin)))  )


library( lattice )
library( cggRutils )

xxx <- bwplot( mcd.alpha.bin ~ 100*cgpred.con | p*n,
	data=hrResults01.cgpredcon,
	strip=strip.custom(strip.names=FALSE,strip.levels=TRUE),
	scales=list(x = list(relation="free")),
	layout=c(1,3),
	as.table=TRUE,
	xlab="PERCENTAGE OF SIMULATED DATA IDENTIFIED AS OUTLYING",
	panel=function(...){
		panel.bwplot(...)
		panel.abline(v=1, lty="dotted", col="black")
	}
)
trellis.device(windows, theme=splus.theme)
#tss <- trellis.par.get("superpose.symbol")
#tss$pch <- 1:10
#trellis.par.set("superpose.symbol",tss)
print(xxx)
dev.off()

trellis.device(pdf, theme=splus.theme, file="PercentageSimulatedDataOutlier01.pdf")
print(xxx)
dev.off()



hrResults05.cgpredcon <- cbind(
	apply(t(hr.cases),2,function(x) rep(x,each=5000)),
	do.call("c",lapply( hrResults05, function(x) x[,"CGPRED.CON",drop=TRUE] ))
)
hrResults05.cgpredcon <- data.frame( hrResults05.cgpredcon, row.names=NULL )
names( hrResults05.cgpredcon ) <- c("p","n","mcd.alpha","mbp","cgpred.con")

hrResults05.cgpredcon$mcd.alpha.bin <- ifelse( hrResults05.cgpredcon$mbp, "MBP", 
  format( hrResults05.cgpredcon$mcd.alpha, 3 ) )
# make ordered factor with MBP first
thesorter <- function(x) { n <- length(x); x[c(n,1:(n-1))] }
hrResults05.cgpredcon$mcd.alpha.bin <- factor( hrResults05.cgpredcon$mcd.alpha.bin, ordered=TRUE,
    levels=thesorter(sort(unique(hrResults05.cgpredcon$mcd.alpha.bin)))  )


xxx <- bwplot( mcd.alpha.bin ~ 100*cgpred.con | p*n,
	data=hrResults05.cgpredcon,
	strip=strip.custom(strip.names=FALSE,strip.levels=TRUE),
	scales=list(x = list(relation="free")),
	layout=c(1,3),
	as.table=TRUE,
	xlab="PERCENTAGE OF SIMULATED DATA IDENTIFIED AS OUTLYING",
	panel=function(...){
		panel.bwplot(...)
		panel.abline(v=5, lty="dotted", col="black")
	}
)
trellis.device(windows, theme=splus.theme)
print(xxx)
dev.off()

trellis.device(pdf, theme=splus.theme, file="PercentageSimulatedDataOutlier05.pdf")
print(xxx)
dev.off()




