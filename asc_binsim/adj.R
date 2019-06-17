rm( list=ls() )
load( 'Rdata/setup.Rdata' )
source( '../code/plot_fxns.R' )

allprev			<- array( NA, dim=c( nx, maxit, 3	), dimnames=list( 1:nx, 1:maxit, c( 'both', '1', '2' ) ) )
allprev_pop	<- array( NA, dim=c( nx, maxit, 3	), dimnames=list( 1:nx, 1:maxit, c( 'both', '1', '2' ) ) )
adjs				<- array( NA, dim=c( nx, maxit, 3	), dimnames=list( 1:nx, 1:maxit, c( 'both', '1', '2' ) ) )
for( it in 1:maxit )
	for( xval in sample(nx) )
try({
	load( paste0( 'Rdata/adj/'	, fs[xval], '_', it, '_', 'ldak', '.Rdata' ) )
	tryCatch({
	allprev_pop	[xval,it,]	<- pop_prevs
	}, error=function(e) file.remove( paste0( 'Rdata/adj/'	, fs[xval], '_', it, '_', 'ldak', '.Rdata' ) ) )
	allprev			[xval,it,]	<- prevs
	adjs				[xval,it,]	<- sapply( 1:3, function(i) adjfxn( P=prevs[i], K=pop_prevs[i] ) )

	rm( prevs, pop_prevs )
	save( adjs, allprev, allprev_pop, file='Rdata/adj.Rdata' )
})
save( adjs, allprev, allprev_pop, file='Rdata/adj.Rdata' )
