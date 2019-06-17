rm( list=ls() )
#library(gxemm)
library(GxEMM)
source( '../code/sim_fxn.R' )
source( '../code/plot_fxns.R' )
load( 'Rdata/setup.Rdata' )

allprev	<- array( NA, dim=c( nx, maxit, 3	), dimnames=list( 1:nx, 1:maxit, c( 'both', '1', '2' ) ) )
for( it in 1:maxit )
	for( xval in sample(nx) )
try({
tryCatch({
	load( paste0( 'Rdata/adj/'	, xs[xval], '_', it, '_', modes[1], '.Rdata' ) )
}, error=function(e) file.remove( paste0( 'Rdata/adj/'	, xs[xval], '_', it, '_', modes[1], '.Rdata' ) ) )
	allprev[xval,it,]	<- prevs
	rm( prevs )
	adjs	<- apply( allprev, 1:3, function(x) adjfxn( x )  )
	save( adjfxn, adjs, allprev, file='Rdata/adj.Rdata' )
})
