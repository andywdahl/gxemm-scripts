library('msm')
gxemm	<- function(
	y, X, K, Z=NULL, type=c('hom','iid','free','full')[1],
	M.all=1, trunc=FALSE,
	env_spec_noise=(type=='free')
){

	if( type == 'hom' ){
		if( !env_spec_noise ){
			gout	<- Rgcta( y=y, X=X, A=list(K), M.all=M.all, trunc=trunc )
			out		<- list( h2=gout$h2s[1], h2.se=gout$ses[1], ll=gout$ll, pval=gout$pval )
		} else {
			gout	<- Rgcta( y=y, X=X, A=c( list(K), lapply( 1:(K0-1), function(k) diag( Z[,k,drop=T]^2	) ) ), M.all=M.all, trunc=trunc )

			h2	<- c( gout$h2s[1], sapply( 1:K0, function(k)
				if( k == K0 ){
					gout$sig2s[1]/sum( gout$sig2s[c(1,1+k,1+K0)] )
				} else {
					gout$sig2s[1]/sum( gout$sig2s[c(1,    1+K0)] )
				}
			))
			deltform	<- lapply( 1:K0, function(k)
				if( k == K0 ){
					as.formula(paste0( '~(x1)/(x1+x', 1+k, '+x', 1+K0, ')'))
				} else {
					as.formula(paste0( '~(x1)/(x1+x',            1+K0, ')'))
				}
			)
			h2Covmat	<- deltamethod( deltform, gout$sig2s[1:(1+K0)], gout$sig2Var, ses=F )

			h2.se	<- c( NA, sqrt( diag( h2Covmat ) ) )

			names(h2)	<- names(h2.se)	<- c( 'hom', paste0( 'h2_', 1:K0 ) )

			out		<- list( h2=h2, h2.se=h2.se, h2Covmat=h2Covmat, ll=gout$ll, pval=gout$pval, sig2s=gout$sig2s[1:(length(gout$sig2s)-1)], sig2Var=gout$sig2Var )

		}

	} else if( type == 'iid' ){
		gout	<- Rgcta( y=y, X=X, A=list(K, K * (Z %*% t(Z))), M.all=M.all, trunc=trunc )
		h2		<- c( gout$h2s[1:2], sum(gout$h2s[1:2]) )
		h2.se	<- c( gout$ses[1:2], deltamethod(
			list( as.formula("~(x1+x2)/(x1+x2+x3)") ),
			gout$sig2s[1:3], gout$sig2Var, ses=T
		))
		names(h2)	<- names(h2.se)	<- c( 'hom', 'het', 'hom+het' )
		out		<- list( h2=h2, h2.se=h2.se, ll=gout$ll, pval=gout$pval )

	} else if( type == 'free' ){

		A	<- c( list(K), lapply( 1:K0, function(k) K * ( Z[,k,drop=F] %*% t(Z[,k,drop=F])	) ) )
		if( env_spec_noise )
			A	<- c( A, lapply( 1:(K0-1), function(k) diag( Z[,k,drop=T]^2	) ) )
		gout<- Rgcta( y=y, X=X, A=A, M.all=M.all, trunc=trunc )

		if( env_spec_noise ){
			h2	<- c( gout$h2s[1], sapply( 1:K0, function(k)
				if( k != K0 ){
					sum( gout$sig2s[c(1,1+k)] )/sum( gout$sig2s[c(1,1+k,1+K0+k,1+K0+K0)] )
				} else {
					sum( gout$sig2s[c(1,1+k)] )/sum( gout$sig2s[c(1,1+k,       1+K0+K0)] )
				}
			))
			deltform	<- lapply( 1:K0, function(k)
				if( k != K0 ){
					as.formula(paste0( '~(x1+x', 1+k, ')/(x1+x', 1+k, '+x', 1+K0+k, '+x', 1+K0+K0, ')'))
				} else {
					as.formula(paste0( '~(x1+x', 1+k, ')/(x1+x', 1+k, '+x',               1+K0+K0, ')'))
				}
			)
			h2Covmat	<- deltamethod( deltform, gout$sig2s[1:(1+2*K0)], gout$sig2Var, ses=F )
		} else {
			h2	<- c( gout$h2s[1], sapply( 1:K0, function(k)
				sum( gout$sig2s[c(1,1+k)] )/sum( gout$sig2s[c(1,1+k,2+K0)] )
			))
			deltform	<- lapply( 1:K0, function(k)
				as.formula(paste0( '~(x1+x', 1+k, ')/(x1+x', 1+k, '+x', 2+K0, ')'))
			)
			h2Covmat	<- deltamethod( deltform, gout$sig2s[1:(2+K0)], gout$sig2Var, ses=F )
		}
		h2.se	<- c( NA, sqrt( diag( h2Covmat ) ) )

		names(h2)	<- names(h2.se)	<- c( 'hom', paste0( 'h2_', 1:K0 ) )

		out		<- list( h2=h2, h2.se=h2.se, h2Covmat=h2Covmat, ll=gout$ll, pval=gout$pval, sig2s=gout$sig2s[1:(length(gout$sig2s)-1)], sig2Var=gout$sig2Var )

	} else if( type == 'full' ){
		stop()
	}
	out
}
