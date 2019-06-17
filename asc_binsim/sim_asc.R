simfxn_asc	<- function( it, X, Z, S, L, sig2hom, sig2het, tauhet, bin=FALSE, prev=.2, fixmean=TRUE, alpha1=NA, sd.alpha=.1, N ){
	set.seed( it )

	Gcaus	<- scale( matrix( rbinom( N_univ*S, 2, .5 ), N_univ, S ) )
	ws		<- colMeans( Z^2 )

	allbetas	<- sapply( sig2het, function(sig) rnorm( S, sd=sqrt( sig/S ) ) )

	if( is.na( tauhet[1] ) ){
		epsilon	<- sqrt(1-sig2hom-sum(ws*sig2het)) * rnorm(nrow(X))
	} else {
		sig2bg	<- 1-sig2hom-sum(ws*sig2het)-sum(ws*tauhet)
		epsilon	<- sqrt(sig2bg + as.numeric((Z^2) %*% tauhet)) * rnorm(nrow(X))
	}

	if( fixmean & bin ){
		alpha		<- rep(0,ncol(X))
		alpha[length(alpha)]<- alpha1
	} else {
		alpha		<- rnorm( ncol(X), sd=sd.alpha )
	}

	allbZ			<- Z %*% t( allbetas )
	y	<- as.numeric(
		X %*% alpha +
		Gcaus %*% rnorm( S, sd=sqrt( sig2hom/S ) ) +
		#sapply( 1:nrow(Z), function(i) Gcaus[i,] %*% allbZ[i,] ) +
		rowSums( Gcaus * allbZ ) +
		epsilon
	)
	yq		<- quantile( y, 1-prev )
	y			<- as.numeric( y > yq )

	prevs	<- c( mean( y ), mean( y[ Z[,1] == 1 ] ), mean( y[ Z[,2] == 1 ] ) )

	#ascertain
	sub			<- c(
		sample( which(y==0), N/2, replace=FALSE ),
		sample( which(y==1), N/2, replace=FALSE )
	)

	Z	<- Z[sub,]
	X	<- X[sub,]
	Gcaus	<- Gcaus[sub,]
	y	<- y[sub]

	G	<- cbind( Gcaus, scale( matrix( rbinom( N*(L-S), 2, .5 ), N, L-S ) ) )
	K	<- 1/ncol(G) * G %*% t(G)

	list( X=X, Z=Z, K=K, y=y, prevs=prevs )

}
