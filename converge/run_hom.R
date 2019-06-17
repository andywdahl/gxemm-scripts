rm( list=ls() )
library(GxEMM)
load( file='Rdata/setup.Rdata' )

savefile	<- paste0( 'Rdata/hom0.Rdata' )
sinkfile	<- paste0( 'Rout/hom0.Rout' )
if( !file.exists( savefile ) &	!file.exists( sinkfile ) ){
	sink( sinkfile )

	## load y, X
	load( 'Rdata/data.Rdata' )
		
	set.seed(28)
	runtime_hom	<- system.time({
		out_hom		<- GxEMM( y, X, K, Z=rnorm(nrow(X)), binary=TRUE, gtype='hom', etype='hom', prev=.088, tmpdir='tmpdir' )
	})[3]
	save( runtime_hom, out_hom, file=savefile )
	rm( y, K, X )
	sink()
}

load( savefile )
print(out_hom$h2)
print(sqrt(out_hom$h2Covmat[1,1]))
