rm( list=ls() )
library(BEDMatrix)
load( 'Rdata/setup.Rdata' )

## load X
load( '../converge/parsed_data/base.Rdata' )
X				<- scale(G)
Xnames	<- rownames(G)
N				<- nrow(X)
rm(G,Y)

#### sample Z
set.seed( 4040 )
len1	<- round(.25*N)
len2	<- N-len1
perm	<- sort.list(runif(N))
subs	<- list( perm[1:len1], perm[len1+1:len2] )

Z	<- matrix( 0, N, 2 )
for( i in 1:K0 )
	Z[subs[[i]],i]	<- 1
Z		<- Z * sqrt(N/sum(Z^2))

### don't include intercept for GCTA, so remove a rando column of Z
X		<- cbind( X, Z[,-2] )


## load K
dat	<- BEDMatrix( paste0( '/ye/zaitlenlabstore/andy/raw_data/converge_data_nov_2017/chr', 1, '.info95maf05.hwe6.mddstress' ) )
rownames(dat)	<- sapply( rownames(dat), function(x) strsplit( x, '_MD' )[[1]][1] )
sub	<- intersect( Xnames, rownames(dat) )
rm( dat )

load( '/ye/zaitlenlabstore/andy/raw_data/converge_data_nov_2017/GRM_K_from_Na/K.Rdata' )
dim(K)
K	<- K[sub,sub]
K	<- K*N/sum(diag(K))
dim(K)
dim(X)
dim(Z)

save( K, subs, X, Xnames, Z, file='Rdata/preprocess.Rdata' )
