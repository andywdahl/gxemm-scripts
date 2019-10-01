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
	savefile3	<- paste0( 'Rdata/free/', xs[xval], '_', it, '_', mode, '.Rdata' )
	savefile5	<- paste0( 'Rdata/iid/', xs[xval], '_', it, '_', mode, '.Rdata' )
	sinkfile	<- paste0( 'Rout/'			, xs[xval], '_', it, '_', mode, '.Rout' )
	if( file.exists( savefile5 ) | file.exists( sinkfile ) )	next
	print(	sinkfile )
	sink(	sinkfile )

	load( 'Rdata/preprocess.Rdata' )              ## load X, Z, K, Xnames
	y	<- simfxn( it, X=cbind(1,X), Z=Z, G=sample_G( seed=it, ncaus, Xnames, lens ), sig2hom=.5, sig2het=rep(0,K0), tauhet=NA, bin=TRUE, alpha1=xs[xval], fixmean=TRUE )

	if( mode == 'pcgc' )
		y	<- y + 1

	prevs	<- mean( y )
	save( prevs, file=savefile0 )

	tmpdir  <- paste0( '/wynton/scratch/gxemm/tmpdir_', mode, '_', it, '_', xval, '_mubin' )

	if( ! file.exists( savefile1 ) ){
		time_hom	<- system.time({
		out_hom		<- GxEMM( y, X, K, Z, binary=( mode == 'pcgc' ), prev=mean( y == 2 ), gtype='hom', tmpdir=tmpdir ) #, noise_K0=TRUE
		})[3]
		save( time_hom, out_hom, file=savefile1 )
	}
	if( ! file.exists( savefile2 ) ){
		time_het	<- system.time({
		out_het		<- GxEMM( y, X, K, Z, binary=( mode == 'pcgc' ), prev=mean( y == 2 ), gtype='iid', tmpdir=tmpdir )
		})[3]
		save( time_het, out_het, file=savefile2 )
	}
	if( ! file.exists( savefile3 ) ){
		time_diag	<- system.time({
		out_diag	<- GxEMM( y, X, K, Z, binary=( mode == 'pcgc' ), prev=mean( y == 2 ), gtype='free', etype='free', tmpdir=tmpdir )
		})[3]
		save( time_diag, out_diag, file=savefile3 )
	}
	if( ! file.exists( savefile5 ) ){
		out_het1	<- GxEMM( y, X, K, Z, binary=( mode == 'pcgc' ), prev=mean( y == 2 ), gtype='iid'	, etype='iid'	, tmpdir=tmpdir )
		save( out_het1, file=savefile5 )
	}

	rm( K, Z, X )
	sink()
	print(warnings())
	print('Done')
})
