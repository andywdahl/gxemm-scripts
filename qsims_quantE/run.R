rm( list=ls() )
library(BEDMatrix)
library(GxEMM)
source( '../code/sim_fxn.R' )
load( 'Rdata/setup.Rdata' )

it	<- as.numeric( commandArgs(TRUE)[[1]] )
set.seed( round(as.numeric(Sys.time())) + it )
for( xval in sample(nx) )
	for( sigtype in sample(5) )
try({
	
	if( sigtype == 4 & xval == 1 ) next

	savefile1	<- paste0( 'Rdata/hom/'	, sigtype, '_', xval, '_', it, '.Rdata' )
	savefile2	<- paste0( 'Rdata/het/'	, sigtype, '_', xval, '_', it, '.Rdata' )
	savefile3	<- paste0( 'Rdata/diag/', sigtype, '_', xval, '_', it, '.Rdata' )
	savefile4	<- paste0( 'Rdata/diag1/', sigtype, '_', xval, '_', it, '.Rdata' )
	savefile5	<- paste0( 'Rdata/hom1/', sigtype, '_', xval, '_', it, '.Rdata' )
	sinkfile	<- paste0( 'Rout/'			, sigtype, '_', xval, '_', it, '.Rout' )
	if( file.exists( savefile5 ) | file.exists( sinkfile ) )	next
	print( sinkfile )
	sink(	sinkfile )

	## load X, Z, K, Xnames
	load( 'Rdata/preprocess.Rdata' )
	
	# generate pheno
	if( sigtype == 4 )
		Z	<- Z/sqrt(2)
	y		<- simfxn( it, X=cbind(1,X), Z,
		G=sample_G( seed=it, ncaus, Xnames, lens ),
		sig2hom=all_sig2homs(xval)[ sigtype],
		sig2het=all_sig2hets(xval)[[sigtype]],
		tauhet=all_tauhets	(xval)[[sigtype]]
	)

	tmpdir  <- paste0( '/wynton/scratch/gxemm/tmpdir_', sigtype, '_', it, '_', xval, '_qsim' )

	if( ! file.exists( savefile1 ) ){
		out_hom		<- GxEMM_dev( y, X, K, Z,gtype='hom', tmpdir=tmpdir, noise_K0=TRUE )
		save( out_hom, file=savefile1 )
	}
	if( ! file.exists( savefile2 ) ){
		out_het		<- GxEMM_dev( y, X, K, Z, gtype='iid', tmpdir=tmpdir, noise_K0=TRUE )
		save( out_het, file=savefile2 )
	}
	if( ! file.exists( savefile3 ) ){
		out_diag	<- GxEMM_dev( y, X, K, Z, gtype='free', etype='free', tmpdir=tmpdir, noise_K0=TRUE )
		save( out_diag, file=savefile3 )
	}
	if( ! file.exists( savefile4 ) ){
		out_diag1	<- GxEMM_dev( y, X, K, Z, gtype='free', etype='hom', tmpdir=tmpdir, noise_K0=TRUE )
		save( out_diag1, file=savefile4 )
	}
	if( ! file.exists( savefile5 ) ){
		out_hom1	<- GxEMM_dev( y, X, K, Z, gtype='hom', etype='free', tmpdir=tmpdir, noise_K0=TRUE )
		save( out_hom1, file=savefile5 )
	}
	rm( K, Z, X )
	print(warnings())
	print('Done')
	sink()
})
