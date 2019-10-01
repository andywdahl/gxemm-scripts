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
	y		<- simfxn( it, X=cbind(1,X), as.matrix(Z),
		G=sample_G( seed=it, ncaus, Xnames, lens ),
		sig2hom=all_sig2homs(xval)[ sigtype],
		sig2het=all_sig2hets(xval)[[sigtype]],
		tauhet=all_tauhets	(xval)[[sigtype]]
	)

	if( ! file.exists( savefile1 ) ){
		runtime_hom	<- system.time({
		out_hom		<- GxEMM( y, X, K, Z,gtype='hom' )
		})[3]
		save( runtime_hom, out_hom, file=savefile1 )
	}
	if( ! file.exists( savefile2 ) ){
		runtime_het	<- system.time({
		out_het		<- GxEMM( y, X, K, Z, gtype='iid' )
		})[3]
		save( runtime_het, out_het, file=savefile2 )
	}
	if( ! file.exists( savefile3 ) ){
		runtime_diag	<- system.time({
		out_diag	<- GxEMM( y, X, K, Z, gtype='free', etype='free' )
		})[3]
		save( runtime_diag, out_diag, file=savefile3 )
	}
	if( ! file.exists( savefile4 ) ){
		runtime_diag1	<- system.time({
		out_diag1	<- GxEMM( y, X, K, Z, gtype='free', etype='hom' )
		})[3]
		save( runtime_diag1, out_diag1, file=savefile4 )
	}
	if( ! file.exists( savefile5 ) ){
		runtime_hom1	<- system.time({
		out_hom1	<- GxEMM( y, X, K, Z, gtype='hom', etype='free' )
		})[3]
		save( runtime_hom1, out_hom1, file=savefile5 )
	}
	rm( K, Z, X )
	print(warnings())
	print('Done')
	sink()
})
