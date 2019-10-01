rm( list=ls() )
library(BEDMatrix)
library(GxEMM)
load( 'Rdata/setup.Rdata' )

### K names
load( '/wynton/scratch/gxemm/converge_data_nov_2017/GRM_K_from_Na/K.Rdata' )
N		<- nrow(K)
IDs	<- rownames(K)
length( IDs )

### G names
dat			<- BEDMatrix( paste0( '/wynton/scratch/gxemm/converge_data_nov_2017/chr', 1, '.info95maf05.hwe6.mddstress' ) )
Gnames	<- sapply( rownames(dat), function(x) strsplit( x, '_MD' )[[1]][1] )
length( Gnames )
all( IDs %in% Gnames )
rm( dat )

#### sample Z
set.seed( 4040 )
Z	<- as.matrix(rnorm(N))
X	<- Z
nrow(Z)

#K	<- K*N/sum(diag(K))

## load Xnames
load( '../converge/parsed_data/base.Rdata' )
Xnames	<- rownames(G)
rm(G,Y)

all( Xnames == IDs )

write.table( cbind( IDs, IDs, 0, 0, 0, 0 ), file='/wynton/scratch/gxemm/conv/tmp.fam', sep=' ', row.names=F, col.names=F, quote=F )

save( Xnames, X, K, IDs, Z, file='Rdata/preprocess.Rdata' )

write_kin_loc	<- function(tmpdir,K,index,ldak_loc, X0, IDs=1:nrow(K) ){
	N	<- nrow(K)
	w	<- mean(diag(K))
	K	<- K/w
	prefix	<- paste0( tmpdir, '/K.', index )
	keep = which( lower.tri(K,diag=T) , arr.ind=T )
	non_mis = rep( 1  , length(keep[,1]) )
	keep = cbind( keep , non_mis , K[ keep ] )
	rm( K ); gc()
	keep = keep[ order(keep[,1],keep[,2]) , ]
	write.table( keep						, file=gzfile(paste0( prefix, '.grm.gz' )), col.names=F, row.names=F, quote=F )
	keep = which( lower.tri(K,diag=T) , arr.ind=T )
	write.table( cbind(1:N,1:N)	, file=       paste0( prefix, '.grm.id' ) , col.names=F, row.names=F, quote=F )
}
write_kin_loc(tmpdir='/wynton/scratch/gxemm/converge_data_nov_2017/',K			,'for_MTG2'			,"~/GxEMM/code/ldak5.linux ", X, IDs=IDs )
