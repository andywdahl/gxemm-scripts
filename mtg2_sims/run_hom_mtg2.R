rm( list=ls() )
library(GxEMM)
source( '../code/sim_fxn.R' )
source( '../code/MRNM.R' )
load( 'Rdata/setup.Rdata' )

it	<- as.numeric( commandArgs(TRUE)[[1]] )
for( xval in sample(nx) ){

	savefile	<- paste0( 'Rdata/mtg2/', xs[xval], '_', it, '.Rdata' )
	sinkfile	<- paste0( 'Rout/mtg2_'	, xs[xval], '_', it, '.Rout' )
	if( file.exists( savefile ) | file.exists( sinkfile ) )	next
	sink(	sinkfile )

	load( 'Rdata/preprocess.Rdata' )              ## load X, Z, K, Xnames
	y	<- simfxn( it, X=cbind(1,X), Z=Z, G=sample_G( seed=it, ncaus, Xnames, lens ), sig2hom=.5, sig2het=rep(0,K0), tauhet=NA, bin=TRUE, alpha1=xs[xval], fixmean=TRUE )

	if( ! file.exists( savefile ) ){
		mrnmtime_hom	<- system.time({
			out_mrnm_hom		<- Rmtg2( y=y, X=X, Z=Z, suffix=paste0( '.', it, '_', xval ), mc.cores=1, type='hom', IDs=IDs )
		})[3]
		print( out_mrnm_hom )
		save( out_mrnm_hom, mrnmtime_hom, file=savefile )
	}

	rm( y, out_mrnm_hom, mrnmtime_hom, Z, X )
	sink()
	print(warnings())
	print('Done')
}#)
