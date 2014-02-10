library( CerioliOutlierDetection )
library( HardinRockeExtension )

N.SIM <- 5000
B.SIM <- 500

system.time(hrResults.p5.n50   <- hrSim(p=5,n=50  ,N=N.SIM,B=B.SIM))/60
system.time(hrResults.p5.n100  <- hrSim(p=5,n=100 ,N=N.SIM,B=B.SIM))/60
system.time(hrResults.p5.n500  <- hrSim(p=5,n=500 ,N=N.SIM,B=B.SIM))/60
system.time(hrResults.p5.n1000 <- hrSim(p=5,n=1000,N=N.SIM,B=B.SIM))/60

save.image()

system.time(hrResults.p10.n50   <- hrSim(p=10,n=50  ,N=N.SIM,B=B.SIM))/60
system.time(hrResults.p10.n100  <- hrSim(p=10,n=100 ,N=N.SIM,B=B.SIM))/60
system.time(hrResults.p10.n500  <- hrSim(p=10,n=500 ,N=N.SIM,B=B.SIM))/60
system.time(hrResults.p10.n1000 <- hrSim(p=10,n=1000,N=N.SIM,B=B.SIM))/60

save.image()

system.time(hrResults.p20.n50   <- hrSim(p=20,n=50  ,N=N.SIM,B=B.SIM))/60
system.time(hrResults.p20.n100  <- hrSim(p=20,n=100 ,N=N.SIM,B=B.SIM))/60
system.time(hrResults.p20.n500  <- hrSim(p=20,n=500 ,N=N.SIM,B=B.SIM))/60
system.time(hrResults.p20.n1000 <- hrSim(p=20,n=1000,N=N.SIM,B=B.SIM))/60

save.image()

if ( FALSE ) {

100*colMeans( hrResults.p5.n50   )
100*colMeans( hrResults.p5.n100  )
100*colMeans( hrResults.p5.n500  )
100*colMeans( hrResults.p5.n1000 )
100*colMeans( hrResults.p10.n50   )
100*colMeans( hrResults.p10.n100  )
100*colMeans( hrResults.p10.n500  )
100*colMeans( hrResults.p10.n1000 )
100*colMeans( hrResults.p20.n50   )
100*colMeans( hrResults.p20.n100  )
100*colMeans( hrResults.p20.n500  )
100*colMeans( hrResults.p20.n1000 )

}

if ( FALSE ) {

100*colStdevs( hrResults.p5.n50   )
100*colStdevs( hrResults.p5.n100  )
100*colStdevs( hrResults.p5.n500  )
100*colStdevs( hrResults.p5.n1000 )
100*colStdevs( hrResults.p10.n50   )
100*colStdevs( hrResults.p10.n100  )
100*colStdevs( hrResults.p10.n500  )
100*colStdevs( hrResults.p10.n1000 )
100*colStdevs( hrResults.p20.n50   )
100*colStdevs( hrResults.p20.n100  )
100*colStdevs( hrResults.p20.n500  )
100*colStdevs( hrResults.p20.n1000 )

}


allmeans <- rbind(
	c(5,50,100*colMeans( hrResults.p5.n50 )),
	c(10,50,100*colMeans( hrResults.p10.n50 )),
	c(20,50,100*colMeans( hrResults.p20.n50 )),
	c(5,100,100*colMeans( hrResults.p5.n100 )),
	c(10,100,100*colMeans( hrResults.p10.n100 )),
	c(20,100,100*colMeans( hrResults.p20.n100 )),
	c(5,500,100*colMeans( hrResults.p5.n500 )),
	c(10,500,100*colMeans( hrResults.p10.n500 )),
	c(20,500,100*colMeans( hrResults.p20.n500 )),
	c(5,1000,100*colMeans( hrResults.p5.n1000 )),
	c(10,1000,100*colMeans( hrResults.p10.n1000 )),
	c(20,1000,100*colMeans( hrResults.p20.n1000 ))
)

colStdevs <- function(x) apply(x, 2, sd)
allstds <- rbind(
	c(5,50,100*colStdevs( hrResults.p5.n50 )),
	c(10,50,100*colStdevs( hrResults.p10.n50 )),
	c(20,50,100*colStdevs( hrResults.p20.n50 )),
	c(5,100,100*colStdevs( hrResults.p5.n100 )),
	c(10,100,100*colStdevs( hrResults.p10.n100 )),
	c(20,100,100*colStdevs( hrResults.p20.n100 )),
	c(5,500,100*colStdevs( hrResults.p5.n500 )),
	c(10,500,100*colStdevs( hrResults.p10.n500 )),
	c(20,500,100*colStdevs( hrResults.p20.n500 )),
	c(5,1000,100*colStdevs( hrResults.p5.n1000 )),
	c(10,1000,100*colStdevs( hrResults.p10.n1000 )),
	c(20,1000,100*colStdevs( hrResults.p20.n1000 ))
)
	
allmeans <- as.data.frame(allmeans)
allstds  <- as.data.frame(allstds)

names(allmeans)[1:2] <- c("p","n")
names(allstds )[1:2] <- c("p","n")

reshape(allmeans[,c("p","n","CHI2.RAW")], direction="wide", idvar="p",timevar="n")
reshape(allmeans[,c("p","n","CHI2.CON")], direction="wide", idvar="p",timevar="n")
reshape(allmeans[,c("p","n","CHI2.SM")], direction="wide", idvar="p",timevar="n")
reshape(allmeans[,c("p","n","HRASY.RAW")], direction="wide", idvar="p",timevar="n")
reshape(allmeans[,c("p","n","HRASY.CON")], direction="wide", idvar="p",timevar="n")
reshape(allmeans[,c("p","n","HRASY.SM")], direction="wide", idvar="p",timevar="n")
reshape(allmeans[,c("p","n","HRPRED.RAW")], direction="wide", idvar="p",timevar="n")
reshape(allmeans[,c("p","n","HRPRED.CON")], direction="wide", idvar="p",timevar="n")
reshape(allmeans[,c("p","n","HRPRED.SM")], direction="wide", idvar="p",timevar="n")


reshape(allstds[,c("p","n","CHI2.RAW")], direction="wide", idvar="p",timevar="n")
reshape(allstds[,c("p","n","CHI2.CON")], direction="wide", idvar="p",timevar="n")
reshape(allstds[,c("p","n","CHI2.SM")], direction="wide", idvar="p",timevar="n")
reshape(allstds[,c("p","n","HRASY.RAW")], direction="wide", idvar="p",timevar="n")
reshape(allstds[,c("p","n","HRASY.CON")], direction="wide", idvar="p",timevar="n")
reshape(allstds[,c("p","n","HRASY.SM")], direction="wide", idvar="p",timevar="n")
reshape(allstds[,c("p","n","HRPRED.RAW")], direction="wide", idvar="p",timevar="n")
reshape(allstds[,c("p","n","HRPRED.CON")], direction="wide", idvar="p",timevar="n")
reshape(allstds[,c("p","n","HRPRED.SM")], direction="wide", idvar="p",timevar="n")

save.image()

