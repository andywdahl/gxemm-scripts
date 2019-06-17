rm( list=ls() )
library(GxEMM)
source( 'gwas_functions.R' ) # for create_X

load( 'data/gwas_phens.Rdata' )

types	<- 'base'
for( phen in sample(gwas_phens) )
	for( type in types )
try({
	set.seed( which( type == types ) )

	savefile1	<- paste0( 'Rdata/hom/'	, phen, '_', type, '.Rdata' )
	savefile2	<- paste0( 'Rdata/iid/'	, phen, '_', type, '.Rdata' )
	savefile3	<- paste0( 'Rdata/free/', phen, '_', type, '.Rdata' )
	sinkfile	<- paste0( 'Rout/'			, phen, '_', type, '.Rout' )
	if( file.exists( savefile3 ) |	file.exists( sinkfile ) ) next
	sink( sinkfile )

	tmpdir	<- paste0( 'tmpdir_', phen, '_', type )
	print( paste0( phen, '_', type )	)

	load( 'data/Y_QC.Rdata' ) ### Y has some crap pruned out relative to Y_raw
	load( 'data/K.Rdata' )
	Z	<- model.matrix( ~ -1 + Y[,'sex'] )

	y	<- Y[,phen]
	y	<- scale( as.numeric( y ) )
		
	### tweak environments
	if( type == 'dummy' ){
		set.seed(1)
		Z	<-  matrix( rbinom( nrow(X), 1, .5 ), ncol=1 )
		Z	<-  cbind( Z, 1 - Z )
	}

	## add Z fixed effects
	X	<- create_X( y=y, all_Xs=Y, phen=phen, rm_sel=FALSE, noadj=FALSE )[,-1,drop=F]   ### create covariate matrix
	sub	<- which( !is.na( rowSums( X ) ) )
	if( ! min( svd( cbind(X,Z)[sub,] )$d ) < .1 )
		stop()

	nas	<- which( rowSums( is.na( cbind( X, y, Z ) ) ) > 0 )
	if( length(nas)>0 ){
		y	<- y[-nas]
		X	<- X[-nas,,drop=F]
		Z	<- Z[-nas,,drop=F]
		K	<- K[-nas,-nas]
	}

	if( ! file.exists( savefile1 ) ){
		runtime_hom	<- system.time({
			out_hom		<- GxEMM( y, X, K, Z, gtype='hom', etype='hom', tmpdir=tmpdir )
		})[3]
		save( runtime_hom, out_hom, file=savefile1 )
	}

	if( ! file.exists( savefile2 ) ){
		runtime_het	<- system.time({
			out_het		<- GxEMM( y, X, K, Z, gtype='iid', etype='hom', tmpdir=tmpdir )
		})[3]
		save( runtime_het, out_het, file=savefile2 )
	}

	if( ! file.exists( savefile3 ) ){
		runtime_diag	<- system.time({
			out_diag	<- GxEMM( y, X, K, Z, gtype='free', etype='free', tmpdir=tmpdir )
		})[3]
		save( runtime_diag, out_diag, file=savefile3 )
	}

	rm( y, K, Z, X )
	print(warnings())
	print('Done')
})
