rm( list=ls() )
load( 'Rdata/setup.Rdata' )

X	<- matrix( rnorm( N_univ*Q ), N_univ, Q )
Z	<- rbinom( N_univ, 1, .25 )
Z	<- cbind( Z, 1-Z )
X	<- cbind( X, Z[,-2] )

save( Z, X, file='Rdata/preprocess.Rdata' )
