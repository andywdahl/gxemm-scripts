rm( list=ls() )
library(GxEMM)
source( './sim_asc.R' )
load( 'Rdata/setup.Rdata' )

it	<- as.numeric( commandArgs(TRUE)[[1]] )
for( mode in sample(modes) )
	for( xval in sample(nx) )
try({
	savefile0	<- paste0( 'Rdata/adj/'	, fs[xval], '_', it, '_', mode, '.Rdata' )
	savefile1	<- paste0( 'Rdata/hom/'	, fs[xval], '_', it, '_', mode, '.Rdata' )
	savefile2	<- paste0( 'Rdata/het/'	, fs[xval], '_', it, '_', mode, '.Rdata' )
	savefile3	<- paste0( 'Rdata/diag/', fs[xval], '_', it, '_', mode, '.Rdata' )
	savefile4	<- paste0( 'Rdata/diag1/',fs[xval], '_', it, '_', mode, '.Rdata' )
	sinkfile	<- paste0( 'Rout/'			, fs[xval], '_', it, '_', mode, '.Rout' )
	if( file.exists( savefile4 ) | file.exists( sinkfile ) )	next
	sink(	sinkfile )

	load( 'Rdata/preprocess.Rdata' )              ## load X, Z, K, Xnames
	simdat		<- simfxn_asc( it, X=cbind(1,X), Z, ncaus, L, sig2hom=.35, sig2het=rep(0,K0), tauhet=NA, bin=TRUE, prev=fs[xval], fixmean=TRUE, alpha1=.4, sd.alpha=NA, N )
	y					<- simdat$y
	X					<- simdat$X[,-1]
	Z					<- simdat$Z
	K					<- simdat$K
	pop_prevs	<- simdat$prevs
	prevs	<- c( mean( y ), mean( y[ Z[,1] == 1 ] ), mean( y[ Z[,2] == 1 ] ) )

	if( mode != 'ldak' ){
		y	<- y + 1
	} else {
		save( pop_prevs, prevs, file=savefile0 )
	}

	if( ! file.exists( savefile1 ) ){
		out_hom		<- GxEMM( y, X, K, Z, binary=( mode != 'ldak' ), prev=ifelse( mode == 'pcgc', fs[xval], NA ), gtype='hom' )
		save( out_hom, file=savefile1 )                                                                 
	}                                                                                                 
	if( ! file.exists( savefile2 ) ){                                                                 
		out_het		<- GxEMM( y, X, K, Z, binary=( mode != 'ldak' ), prev=ifelse( mode == 'pcgc', fs[xval], NA ), gtype='iid' )
		save( out_het, file=savefile2 )                                                                 
	}                                                                                                 
	if( ! file.exists( savefile3 ) ){                                                                 
		out_diag	<- GxEMM( y, X, K, Z, binary=( mode != 'ldak' ), prev=ifelse( mode == 'pcgc', fs[xval], NA ), gtype='free' )
		save( out_diag, file=savefile3 )                                                                
	}                                                                                                 
	if( ! file.exists( savefile4 ) ){                                                                 
		out_diag1	<- GxEMM( y, X, K, Z, binary=( mode != 'ldak' ), prev=ifelse( mode == 'pcgc', fs[xval], NA ), gtype='free', etype='free' )
		save( out_diag1, file=savefile4 )
	}
	rm( K, Z, X )
	sink()
	print(warnings())
	print('Done')
})
