#allmeans.hrsimold01 allmeans.hrsimold05  
#allmeans.hrsimnew01 allmeans.hrsimnew05  

allmeans01 <- merge(
	allmeans.hrsimnew01,
	allmeans.hrsimold01,
	by=c("p","n","mcd.alpha"),
	all=TRUE,
	suffixes=c("NEW","OLD")
)

head(allmeans01)

allmeans01[ ,c("p","n","mcd.alpha","CGPRED.CONNEW","CGPRED.CONOLD")]


allmeans05 <- merge(
	allmeans.hrsimnew05,
	allmeans.hrsimold05,
	by=c("p","n","mcd.alpha"),
	all=TRUE,
	suffixes=c("NEW","OLD")
)

head(allmeans05)
allmeans05$mcd.alpha.bin <- ifelse( allmeans05$mbpNEW, "MBP", format( allmeans05$mcd.alpha, 3 ) )

allmeans05[ ,c("p","n","mcd.alpha","CGPRED.CONNEW","CGPRED.CONOLD")]

library(lattice)
xyplot( CGPRED.CONNEW ~ CGPRED.CONOLD | p*n,
	groups = mcd.alpha.bin,
	data = allmeans05,
	#layout = c( 4, 4 ),
	panel = function(...) {
		panel.superpose(...)
		panel.abline(c(0,1), lty="dashed", color="black")
		panel.abline(h=5,lty="dashed",color="black")
		panel.abline(v=5,lty="dashed",color="black")
	},
	auto.key = list(space="top")
)

