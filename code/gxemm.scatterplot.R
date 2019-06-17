gxemm.scatterplot	<- function(
	h2s,
	allps,
	labs,
	minp=.01,
	#nlab=length(labs),
	lims=lapply( 1:4, function(ii) c(-.2,1.2) ),
	xlabs=expression( h[greml]^2	, h[greml]^2, h[greml]^2	, h[1]^2 ),
	ylabs=expression( h[iid]^2		, h[1]^2	, h[1]^2	, h[2]^2 )#,
	#cols=c( 1, 1, 'grey', 'orange', 2, 3, 1, 2, 3 )
){
	nh2	<- 9
	y		<- as.numeric( t(h2s) )
	for( ii in 1:4 ){
		if( ii == 1 ){
			x	<- h2s[,'g']
			y	<- h2s[,'iid']
			pvec	<- allps[,'iid']
		} else if( ii == 2 ){
			x	<- h2s[,'g']
			y	<- h2s[,'d1']
			pvec	<- allps[,'d1']
		} else if( ii == 3 ){
			x	<- h2s[,'g']
			y	<- h2s[,'d2']
			pvec	<- allps[,'d2']
		} else if( ii == 4 ){
			x	<- h2s[,'d1']
			y	<- h2s[,'d2']
			pvec	<- allps[,'d2']*0
		}

		plot( x, y, xlim=lims[[ii]], ylim=lims[[ii]], type='n', xlab=xlabs[ii], ylab=ylabs[ii], main='' )
		abline( a=0, b=1, lty=2, col=2 )

		if( length( sub	<- which( pvec < minp ) ) > 0 ){
			#sub	<- sort.list( abs(y-x), dec=T )[1:nlab]
			points( x[-sub], y[-sub]		, pch=16 )
			text(		x[ sub], y[ sub]		, lab=labs[sub] )
			if( sum(pvec!=0,na.rm=T)!=0 )
			text(		x[ sub], y[ sub]-.03, lab=format( pvec, scientific=T, digits=2 )[sub] )
		} else {
			points( x, y, pch=16 )
		}

		if( ii == 4 )
			if( length( sub	<- which( allps[,'delg'] < .01 ) ) > 0 )
				text(		x[ sub], y[ sub]-.15, lab=format( allps[,'delg'] , scientific=T, digits=2 )[sub], col=2 )

	}
}
