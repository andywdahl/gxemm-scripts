rm( list=ls() )
library(GxEMM)
load( file='Rdata/setup.Rdata' )

for( pp in rev(ppsub) ){
	savefile1	<- paste0( 'Rdata/hom/'	, pp, '.Rdata' )
	savefile2	<- paste0( 'Rdata/het/'	, pp, '.Rdata' )
	savefile3	<- paste0( 'Rdata/diag/', pp, '.Rdata' )
	sinkfile	<- paste0( 'Rout/'			, pp, '.Rout' )
	tmpdir		<- paste0( 'tmpdir_', pp )
	if( file.exists( savefile3 ) |	file.exists( sinkfile ) ) next
	sink( sinkfile )

	## load y, X, Z
	load( 'Rdata/data.Rdata' )
	if( pp > ncol(Yb) ){
		Z		<- Y[,pp-ncol(Yb),drop=T]
		Z		<- (Z-min(Z))/(max(Z)-min(Z))
		Z		<- cbind( Z, 1-Z )
	} else {
		Z		<- Yb2Z( Yb[,pp] )
	}
	rm( Y, Yb )
		
	## add Z fixed effects
	X	<- cbind( X, X * ( Z[,1] %o% rep(1,ncol(X)) ) )
	X	<- cbind( X, Z[,1] )

	if( ! file.exists( savefile1 ) ){
		runtime_hom	<- system.time({
			out_hom		<- GxEMM( y, X, K, Z, binary=TRUE, gtype='hom', etype='hom', prev=.088, tmpdir=tmpdir )
		})[3]
		save( runtime_hom, out_hom, file=savefile1 )
	}

	if( ! file.exists( savefile2 ) ){
		runtime_het	<- system.time({
			out_het		<- GxEMM( y, X, K, Z, binary=TRUE, gtype='iid', etype='hom', prev=.088, tmpdir=tmpdir )
		})[3]
		save( runtime_het, out_het, file=savefile2 )
	}

	if( ! file.exists( savefile3 ) ){
		runtime_diag	<- system.time({
			out_diag	<- GxEMM( y, X, K, Z, binary=TRUE, gtype='free', etype='free', prev=.088, tmpdir=tmpdir )
		})[3]
		save( runtime_diag, out_diag, file=savefile3 )
	}

	rm( y, K, Z, X )
	sink()
}
