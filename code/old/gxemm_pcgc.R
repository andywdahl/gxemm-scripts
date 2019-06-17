library('msm')
gxemm_pcgc	<- function(
	y, X, K, Z=NULL, type=c('hom','iid','free','full')[1],
	M.all=1, trunc=FALSE,
	env_spec_noise=(type=='free')
){

	if( type == 'hom' & !env_spec_noise ){
		gout	<- Rgcta( y=y, X=X, A=list(K), M.all=M.all, trunc=trunc )
		h2		<- gout$h2s[1]
		h2.se	<- gout$ses[1]
		h2Covmat	<- h2.se^2

	} else if( type == 'hom' & env_spec_noise ){
		gout	<- Rgcta( y=y, X=X, A=c( list(K), lapply( 1:(K0-1), function(k) diag( Z[,k,drop=T]^2	) ) ), M.all=M.all, trunc=trunc )
		h2	<- c( gout$h2s[1], sapply( 1:K0, function(k)
			ifelse( k == K0,
				gout$sig2s[1]/sum( gout$sig2s[c(1,1+k,1+K0)] ),
				gout$sig2s[1]/sum( gout$sig2s[c(1,    1+K0)] )
			)))
		deltform	<- lapply( 1:K0, function(k)
			as.formula(ifelse( k == K0,
				paste0( '~(x1)/(x1+x', 1+k, '+x', 1+K0, ')'),
				paste0( '~(x1)/(x1+x',            1+K0, ')')
			))
		h2Covmat	<- deltamethod( deltform, gout$sig2s[1:(1+K0)], gout$sig2Var, ses=F )
		h2.se	<- c( NA, sqrt( diag( h2Covmat ) ) )
		names(h2)	<- names(h2.se)	<- c( 'hom', paste0( 'h2_', 1:K0 ) )

	} else if( type == 'iid' ){
		gout	<- Rgcta( y=y, X=X, A=list(K, K * (Z %*% t(Z))), M.all=M.all, trunc=trunc )
		h2		<- c( gout$h2s[1:2], sum(gout$h2s[1:2]) )
		h2.se	<- c( gout$ses[1:2], deltamethod( list( as.formula("~(x1+x2)/(x1+x2+x3)") ), gout$sig2s[1:3], gout$sig2Var, ses=F ))
		names(h2)	<- names(h2.se)	<- c( 'hom', 'het', 'hom+het' )

	} else if( type == 'free' & env_spec_noise ){
		A	<- c( list(K),
			lapply( 1:K0, function(k) K * ( Z[,k,drop=F] %*% t(Z[,k,drop=F])	) ),
			lapply( 1:(K0-1), function(k) diag( Z[,k,drop=T]^2	) ) )
		gout<- Rgcta( y=y, X=X, A=A, M.all=M.all, trunc=trunc )
		h2	<- c( gout$h2s[1], sapply( 1:K0, function(k)
			ifelse( k != K0,
				sum( gout$sig2s[c(1,1+k)] )/sum( gout$sig2s[c(1,1+k,1+K0+k,1+K0+K0)] ),
				sum( gout$sig2s[c(1,1+k)] )/sum( gout$sig2s[c(1,1+k,       1+K0+K0)] )
			)))
		deltform	<- lapply( 1:K0, function(k)
			as.formula(ifelse( k != K0,
				paste0( '~(x1+x', 1+k, ')/(x1+x', 1+k, '+x', 1+K0+k, '+x', 1+K0+K0, ')'),
				paste0( '~(x1+x', 1+k, ')/(x1+x', 1+k, '+x',               1+K0+K0, ')')
			)))
			#ifelse( k != K0,
			#	as.formula(paste0( '~(x1+x', 1+k, ')/(x1+x', 1+k, '+x', 1+K0+k, '+x', 1+K0+K0, ')')),
			#	as.formula(paste0( '~(x1+x', 1+k, ')/(x1+x', 1+k, '+x',               1+K0+K0, ')'))
			#))
		h2Covmat	<- deltamethod( deltform, gout$sig2s[1:(1+2*K0)], gout$sig2Var, ses=F )
		h2.se	<- c( NA, sqrt( diag( h2Covmat ) ) )
		names(h2)	<- names(h2.se)	<- c( 'hom', paste0( 'h2_', 1:K0 ) )

	} else if( type == 'free' & !env_spec_noise ){
		A	<- c( list(K), lapply( 1:K0, function(k) K * ( Z[,k,drop=F] %*% t(Z[,k,drop=F])	) ) )
		gout<- Rgcta( y=y, X=X, A=A, M.all=M.all, trunc=trunc )
		h2	<- c( gout$h2s[1], sapply( 1:K0, function(k)
			sum( gout$sig2s[c(1,1+k)] )/sum( gout$sig2s[c(1,1+k,2+K0)] )
		))
		deltform	<- lapply( 1:K0, function(k)
			as.formula(paste0( '~(x1+x', 1+k, ')/(x1+x', 1+k, '+x', 2+K0, ')'))
		)
		h2Covmat	<- deltamethod( deltform, gout$sig2s[1:(2+K0)], gout$sig2Var, ses=F )

		h2.se	<- c( NA, sqrt( diag( h2Covmat ) ) )
		names(h2)	<- names(h2.se)	<- c( 'hom', paste0( 'h2_', 1:K0 ) )

	} else if( type == 'full' ){
		stop()
	}
	df	<- length(gout$sig2s)-1
	list( h2=h2, h2.se=h2.se, h2Covmat=h2Covmat, df=df,
		ll=gout$ll, pval=gout$pval, sig2s=gout$sig2s[1:df], sig2Var=gout$sig2Var )
}
