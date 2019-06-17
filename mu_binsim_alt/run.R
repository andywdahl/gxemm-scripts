rm( list=ls() )
library(GxEMM)
source( '../code/sim_fxn.R' )
load( 'Rdata/setup.Rdata' )

it	<- as.numeric( commandArgs(TRUE)[[1]] )
for( mode in sample(modes) )
	for( xval in sample(nx) )
try({
	savefile0	<- paste0( 'Rdata/adj/'	, xs[xval], '_', it, '_', mode, '.Rdata' )
	savefile1	<- paste0( 'Rdata/hom/'	, xs[xval], '_', it, '_', mode, '.Rdata' )
	savefile2	<- paste0( 'Rdata/het/'	, xs[xval], '_', it, '_', mode, '.Rdata' )
	savefile3	<- paste0( 'Rdata/diag/', xs[xval], '_', it, '_', mode, '.Rdata' )
	savefile4	<- paste0( 'Rdata/diag1/',xs[xval], '_', it, '_', mode, '.Rdata' )
	sinkfile	<- paste0( 'Rout/'			, xs[xval], '_', it, '_', mode, '.Rout' )
	if( file.exists( savefile4 ) | file.exists( sinkfile ) )	next
	sink(	sinkfile )

	load( 'Rdata/preprocess.Rdata' )              ## load X, Z, K, Xnames
	G	<- sample_G( seed=it, ncaus, Xnames, lens ) ## load G caus
	# generate pheno
	y		<- simfxn( it, X=cbind(1,X), Z, G, sig2hom=.15, sig2het=c( rev(xs)[xval],  xs[xval] ), tauhet=NA, bin=TRUE, alpha1=.4, fixmean=TRUE )
	rm( G )

	if( mode == 'pcgc' )
		y	<- y + 1

	prevs	<- c( mean( y ), mean( y[ Z[,1] == 1 ] ), mean( y[ Z[,2] == 1 ] ) )
	save( prevs, file=savefile0 )

	tmpdir  <- paste0( '/wynton/scratch/gxemm/tmpdir_', it, '_', mode, '_', xval, '_mu_binsim_alt' )

	if( ! file.exists( savefile1 ) ){
		out_hom		<- GxEMM( y, X, K, Z, binary=( mode == 'pcgc' ), prev=mean( y == 2 ), gtype='hom', tmpdir=tmpdir )
		save( out_hom, file=savefile1 )
	}
	if( ! file.exists( savefile2 ) ){
		out_het		<- GxEMM( y, X, K, Z, binary=( mode == 'pcgc' ), prev=mean( y == 2 ), gtype='iid', tmpdir=tmpdir )
		save( out_het, file=savefile2 )
	}
	if( ! file.exists( savefile3 ) ){
		out_diag	<- GxEMM( y, X, K, Z, binary=( mode == 'pcgc' ), prev=mean( y == 2 ), gtype='free', tmpdir=tmpdir )
		save( out_diag, file=savefile3 )
	}
	if( ! file.exists( savefile4 ) ){
		out_diag1	<- GxEMM( y, X, K, Z, binary=( mode == 'pcgc' ), prev=mean( y == 2 ), gtype='free', etype='free', tmpdir=tmpdir )
		save( out_diag1, file=savefile4 )
	}
	rm( K, Z, X )
	sink()
	print(warnings())
	print('Done')
})
