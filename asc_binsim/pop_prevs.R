rm( list=ls() )
library(GxEMM)
source( './sim_asc.R' )
load( 'Rdata/setup.Rdata' )

mode	<- 'ldak'
for( it in sample(maxit) )
	for( xval in sample(nx) )
try({
	savefile0	<- paste0( 'Rdata/adj/'	, fs[xval], '_', it, '_', mode, '.Rdata' )
	if( file.exists( savefile0 ) )
		next

	print( c( it, mode, xval ) )

	load( 'Rdata/preprocess.Rdata' )              ## load X, Z, K, Xnames
	simdat	<- simfxn_asc( it, X=cbind(1,X), Z, ncaus, L, sig2hom=.35, sig2het=rep(0,K0), tauhet=NA, bin=TRUE, prev=fs[xval], fixmean=TRUE, alpha1=.4, sd.alpha=NA, N )
	y				<- simdat$y
	X				<- simdat$X[,-1]
	Z				<- simdat$Z
	K				<- simdat$K

	prevs	<- c( mean( y ), mean( y[ Z[,1] == 1 ] ), mean( y[ Z[,2] == 1 ] ) )
	pop_prevs	<- simdat$prevs
	save( prevs, pop_prevs, file=savefile0 )
	rm( K, Z, X )

})
