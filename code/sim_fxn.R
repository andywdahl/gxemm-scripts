library(BEDMatrix)

simfxn	<- function( it, X, Z, G, sig2hom, sig2het, tauhet, bin=FALSE, prev=.2, fixmean=TRUE, alpha1=NA, sd.alpha=.1 ){
	set.seed( it )

	S		<- ncol(G)
	ws	<- colMeans( Z^2 )

	allbetas	<- sapply( sig2het, function(sig) rnorm( S, sd=sqrt( sig/S ) ) )

	if( is.na( tauhet[1] ) ){
		epsilon	<- sqrt(1-sig2hom-sum(ws*sig2het)) * rnorm(nrow(G))
	} else {
		sig2bg	<- 1-sig2hom-sum(ws*sig2het)-sum(ws*tauhet)
		epsilon	<- sqrt(sig2bg + as.numeric((Z^2) %*% tauhet)) * rnorm(nrow(G))
	}

	if( fixmean & bin ){
		alpha		<- rep(0,ncol(X))
		alpha[length(alpha)]<- alpha1
	} else {
		alpha		<- rnorm( ncol(X), sd=sd.alpha )
	}

	y	<- as.numeric(
		X %*% alpha +
		G %*% rnorm( S, sd=sqrt( sig2hom/S ) ) +
		sapply( 1:nrow(Z), function(i) G[i,] %*% ( allbetas %*% Z[i,] ) ) +
		epsilon
	)

	if( bin )
		y				<- sapply( y, function(x) ifelse( x>quantile( y, 1-prev ), 1, 0 ) )
	y
}

sample_G	<- function(seed,ncaus,Xnames,lens){
	set.seed( seed )

	caus.chr	<- sample( 1:22, size=ncaus, rep=T, prob=lens )
	ncaus.chr	<- sapply( 1:22, function(chr) sum( caus.chr == chr ) )
	rm( caus.chr )

	G					<- NULL
	for( chr in 1:22 ){
		if( ncaus.chr[chr] == 0 ) next
		suppressMessages(
			dat	<- BEDMatrix( paste0( '/ye/zaitlenlabstore/andy/raw_data/converge_data_nov_2017/chr', chr, '.info95maf05.hwe6.mddstress' ) )
		)
		rownames(dat)	<- sapply( rownames(dat), function(x) strsplit( x, '_MD' )[[1]][1] )
		sub	<- intersect( Xnames, rownames(dat) )
		stopifnot( all.equal( sub, Xnames ) )
		G		<- cbind( G, dat[sub,sample( ncol(dat), ncaus.chr[chr] )] )
		rm( dat )
	}
	G			<- scale(G)
	G[is.na(G)]	<- 0
	G			<- scale(G)
	badsub	<- which( colMeans( is.na(G) ) == 1 )
	if( length(badsub) > 0 )
		stop( 'badsub:' )
	G
}
