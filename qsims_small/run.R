rm( list=ls() )
library(BEDMatrix)
library(GxEMM)
source( '../code/sim_fxn.R' )
load( '../qsims/Rdata/setup.Rdata' )


it	<- as.numeric( commandArgs(TRUE)[[1]] )
for( xval in sample(nx) )try({

	for( sigtype in sample(5) )
	
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
	y		<- simfxn( it, X=cbind(1,X), Z,
		G=sample_G( seed=it, ncaus, Xnames, lens ),
		sig2hom=all_sig2homs(xval)[ sigtype],
		sig2het=all_sig2hets(xval)[[sigtype]],
		tauhet=all_tauhets	(xval)[[sigtype]]
	)
	set.seed(it*7+xval)
	mysub	<- sample( length(y), 1200 )
	y	<- y[mysub]
	X	<- X[mysub,]
	K	<- K[mysub,mysub]
	Z	<- Z[mysub,]

	tmpdir  <- paste0( '/wynton/scratch/gxemm/tmpdir_', sigtype, '_', it, '_', xval, '_smallN' )

	if( ! file.exists( savefile1 ) ){
		out_hom		<- GxEMM( y, X, K, Z,gtype='hom', tmpdir=tmpdir )
		save( out_hom, file=savefile1 )
	}
	if( ! file.exists( savefile2 ) ){
		out_het		<- GxEMM( y, X, K, Z, gtype='iid', tmpdir=tmpdir )
		save( out_het, file=savefile2 )
	}
	if( ! file.exists( savefile3 ) ){
		out_diag	<- GxEMM( y, X, K, Z, gtype='free', etype='free', tmpdir=tmpdir )
		save( out_diag, file=savefile3 )
	}
	if( ! file.exists( savefile4 ) ){
		out_diag1	<- GxEMM( y, X, K, Z, gtype='free', etype='hom', tmpdir=tmpdir )
		save( out_diag1, file=savefile4 )
	}
	if( ! file.exists( savefile5 ) ){
		out_hom1	<- GxEMM( y, X, K, Z, gtype='hom', etype='free', tmpdir=tmpdir )
		save( out_hom1, file=savefile5 )
	}
	rm( K, Z, X )
	print(warnings())
	print('Done')
	sink()
})
