head(hr.cmstar.results.all)
require(lattice)

aggregate( m ~ n + p + mcd.alpha, data=hr.cmstar.results.all,
  function(x) c(median(x), max(x)-min(x)) ) 

mrange <- aggregate( m ~ n + p + mcd.alpha, data=hr.cmstar.results.all,
  function(x) max(x)-min(x) ) 	

xyplot( m ~ n | p,
	groups = mcd.alpha,
	data = mrange,
	ylab = "Max - Min",
	xlab = "Sample Size",
	auto.key = list("top")
)

xyplot( m ~ p | n,
	groups = mcd.alpha,
	data = mrange,
	ylab = "Max - Min",
	xlab = "Sample Size",
	auto.key = list("top")
)



head(hr.cmstarbig.results.all)
require(lattice)

aggregate( m ~ n + p + mcd.alpha, data=hr.cmstarbig.results.all,
  function(x) c(median(x), max(x)-min(x)) ) 

mrangebig <- aggregate( m ~ n + p + mcd.alpha, data=hr.cmstarbig.results.all,
  function(x) max(x)-min(x) ) 	

xyplot( m ~ n | p,
	groups = mcd.alpha,
	data = mrangebig,
	ylab = "Max - Min",
	xlab = "Sample Size",
	auto.key = list("top")
)

xyplot( m ~ p | n,
	groups = mcd.alpha,
	data = mrangebig,
	ylab = "Max - Min",
	xlab = "Sample Size",
	auto.key = list("top")
)


xyplot( m ~ n | mcd.alpha,
	groups = p,
	data = mrangebig,
	ylab = "Max - Min",
	xlab = "Sample Size",
	auto.key = list("top")
)


