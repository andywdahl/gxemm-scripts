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
	sinkfile	<- paste0( 'Rout/'			, sigtype, '_', xval, '_', it, '.Rout' )
	if( file.exists( savefile3 ) | file.exists( sinkfile ) )	next
	print( sinkfile )
	sink(	sinkfile )

	## load X, Z, K, Xnames
	load( 'Rdata/preprocess.Rdata' )
	
	# generate pheno
	y		<- simfxn( it, X=cbind(1,X), Z,
		G=sample_G( seed=it, ncaus, Xnames, lens ),
		sig2hom=all_sig2homs(xval)[ sigtype],
		sig2het=all_sig2hets(xval)[[sigtype]],
		tauhet=all_tauhets	(xval)[[sigtype]],
		fixmean=FALSE,
		sd.alpha=.1/(1+ncol(X)),
		bin=TRUE, alpha1=.5
	) + 1 # LDAK expects 1/2 coding

	tmpdir  <- paste0( '/wynton/scratch/gxemm/tmpdir_', sigtype, '_', it, '_', xval, '_binsim' )

	if( ! file.exists( savefile1 ) ){
		out_hom		<- GxEMM( y, X, K, Z, binary=TRUE, prev=mean( y == 2 ), gtype='hom', tmpdir=tmpdir )
		save( out_hom, file=savefile1 )
	}
	if( ! file.exists( savefile2 ) ){
		out_het		<- GxEMM( y, X, K, Z, binary=TRUE, prev=mean( y == 2 ), gtype='iid', tmpdir=tmpdir )
		save( out_het, file=savefile2 )
	}
	if( ! file.exists( savefile3 ) ){
		out_diag	<- GxEMM( y, X, K, Z, binary=TRUE, prev=mean( y == 2 ), gtype='free', etype='free', tmpdir=tmpdir )
		save( out_diag, file=savefile3 )
	}
	rm( K, Z, X )
	print(warnings())
	print('Done')
	sink()
})
