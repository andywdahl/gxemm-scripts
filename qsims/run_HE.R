rm( list=ls() )
library(BEDMatrix)
library(GxEMM)
source( '../code/sim_fxn.R' )
load( 'Rdata/setup.Rdata' )

it	<- as.numeric( commandArgs(TRUE)[[1]] )
set.seed( round(as.numeric(Sys.time())) + it )
for( sigtype in sample(5) )
	for( xval in sample(nx) )
try({
		for( sigtype in sample(5) )
	
	if( sigtype == 4 & xval == 1 ) next

	savefile1	<- paste0( 'Rdata_HE/hom/'	, sigtype, '_', xval, '_', it, '.Rdata' )
	savefile2	<- paste0( 'Rdata_HE/het/'	, sigtype, '_', xval, '_', it, '.Rdata' )
	savefile3	<- paste0( 'Rdata_HE/diag/', sigtype, '_', xval, '_', it, '.Rdata' )
	savefile4	<- paste0( 'Rdata_HE/diag1/', sigtype, '_', xval, '_', it, '.Rdata' )
	sinkfile	<- paste0( 'Rout_HE/'			, sigtype, '_', xval, '_', it, '.Rout' )
	if( file.exists( savefile4 ) | file.exists( sinkfile ) )	next
	print( sinkfile )
	sink(	sinkfile )

	## load X, Z, K, Xnames
	load( 'Rdata/preprocess.Rdata' )
	
	# generate pheno
	y		<- simfxn( it, X=cbind(1,X), Z,
		G=sample_G( seed=it, ncaus, Xnames, lens ),
		sig2hom=all_sig2homs(xval)[ sigtype],
		sig2het=all_sig2hets(xval)[[sigtype]],
		tauhet=all_tauhets	(xval)[[sigtype]]
	)

	if( ! file.exists( savefile1 ) ){
		out_hom		<- GxEMM_HE( y, X, K, Z,gtype='hom', use_diag=TRUE )
		save( out_hom, file=savefile1 )
	}
	if( ! file.exists( savefile2 ) ){
		out_het		<- GxEMM_HE( y, X, K, Z, gtype='iid', use_diag=TRUE )
		save( out_het, file=savefile2 )
	}
	if( ! file.exists( savefile3 ) ){
		out_diag	<- GxEMM_HE( y, X, K, Z, gtype='free', etype='free', use_diag=TRUE )
		save( out_diag, file=savefile3 )
	}
	if( ! file.exists( savefile4 ) ){
		out_diag1	<- GxEMM_HE( y, X, K, Z, gtype='free', etype='hom', use_diag=TRUE )
		save( out_diag1, file=savefile4 )
	}
	rm( K, Z, X )
	print(warnings())
	print('Done')
	sink()
})
