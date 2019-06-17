library('msm')
mylines	<- function(x,y,col,points=F,trunc=T,pt.cex=2.0,ysd,asts){

	if( trunc ){
	pchs	<- sapply( y, function(y.i) ifelse( y.i <= 0 | y.i > 1, 17, 16 ) )
	y			<- sapply( y, function(y.i) ifelse( y.i < 0, 0, y.i ) )
	y			<- sapply( y, function(y.i) ifelse( y.i > 1, 1, y.i ) )
	} else {
	pchs	<- sapply( y, function(y.i) 16 )
	}
	if( points ){
		points( x	 , y, col=col, lwd=1, lty=1, pch=pchs, cex=pt.cex )
	} else {
		lines( x	 , y, col=col, lwd=1, lty=1 )
	}

	if(!missing(ysd)){
		lines( rep(x,2)		,y+c(-2,2)  /2*ysd,col=col)
		lines( x+c(-1,1)/20,y+rep(2,2) /2*ysd,col=col)
		lines( x+c(-1,1)/20,y+rep(-2,2)/2*ysd,col=col)
	}
	if(!missing(asts))
		text( x, y+.2, lab=c( '', '*', '**', '***' )[asts+1], cex=2 )

}

adjfxn	<- function( P, K=P ){
	phi		<- dnorm( qnorm( 1-K ) )
	(K^2*(1-K)^2) / ( P*(1-P)*phi^2 )
}

lines1	<- function( x, y, col=1, cex=2.4, lwd=3.5, y.max, pch=NA, lty=1 ){
	if( !missing( y.max ) ){
		pchs	<- sapply( y, function(y.i) ifelse( y.i > y.max, 17, 16 ) )
		y			<- sapply( y, function(y.i) ifelse( y.i > y.max, y.max, y.i ) )
	} else {
		pchs	<- rep( pch, length(y) )
	}
	points( x, y, col=col, pch=pchs, cex=cex )
	lines(	x, y, col=col, lwd=lwd, lty=lty )
}
pointline	<- function( x, y, col=1, cex=4, lwd=1, pch=16, nolines=F, lty=1 ){
	y0	<- mean( y, na.rm=T )
	#y0	<- median( y, na.rm=T )
	sd0	<- sd(	 y, na.rm=T ) 
	points( x							, y0								, col=col, pch=pch, cex=cex )
	if(nolines) return()
	lines( rep( x, 2 )		, y0 + c( -1,1)*sd0	, col=col, lwd=lwd, lty=lty )
	lines( x + c(-.02,.02), y0 + rep(1,2)*sd0	, col=col, lwd=lwd, lty=lty )
	lines( x + c(-.02,.02), y0 - rep(1,2)*sd0	, col=col, lwd=lwd, lty=lty )
}


#h2_equal_test	<- function( sigmas, sigmaCov, prevs_sam, prevs_pop=prevs_sam, adj=T ){
#
#	K0			<- length( sigmas ) - 2
#
#	if( adj == T ){
#	locadj	<- sapply( 1:K0, function(kk) adjfxn( prevs_sam[kk], prevs_pop[kk] )  )
#	} else {
#	locadj	<- sapply( 1:K0, function(kk) 1 )
#	}
#
#	if( K0 == 2 ){
#		h2Covmat <- deltamethod( list(
#			as.formula("~(x1+x2)/(x1+x2+x4)"),
#			as.formula("~(x1+x3)/(x1+x3+x4)")
#		), sigmas, sigmaCov, ses=F )
#
#		h2s	<- sapply( 1:K0, function(k) sum( out_diag$sig2s[c(1,1+k)] )/sum( out_diag$sig2s[c(1,1+k,2+K0)] ) )
#
#		Tmat		<- cbind( 1, -1 ) %*% diag( locadj )
#		Tmath2	<- Tmat %*% h2s
#		chi2val	<- as.numeric( solve( Tmat %*% h2Covmat %*% t(Tmat) ) %*% Tmath2 ) %*% Tmath2 
#
#		#print( locadj * h2s )
#		#print( h2s )
#		#print( Tmat %*% h2Covmat %*% t(Tmat) )
#		#print( diag( locadj ) %*% h2Covmat %*%diag( locadj ) )
#
#
#	} else {
#		#Tmat		<- cbind( diag(2), c(-1,-1) ) %*% diag( locadj )
#		stop(K0)
#	}
#
#	pchisq( chi2val, df=K0-1, lower.tail=F )
#
#}


plot.gxemm	<- function(
	h2greml, h2iid, sig2free, h2free,
	se.h2greml, se.h2iid, se.h2iidsum, se.sig2free, se.h2free,
	p.h2greml, p.h2hom, p.h2het, p.h2iid, p.sig2free,
	pdfname,
	xlab=expression( h['GREML']^2, h['Hom']^2, h['IID']^2, h['Hom']^2+h['IID']^2, h[stress]^2, h[unstress]^2, sigma[g]^2, v[stress], v[unstress], w[stress], w[unstress] ),
	cols=c( 1, 1, '#FFBB00', '#FB6542', 1, '#3F681C', '#375E97', '#3F681C', '#375E97' ),
	ylim=c(-.3,1.02)
){

	if( !missing(pdfname) )
		pdf( pdfname, width=6.8, height=5 )
	par( mar=c(8,4.5,.5,.5) )

	plot( c( 1, 3.1 ), ylim, type='n', axes=F, xlab='', ylab='Variance Explained', bty='n', cex.lab=1.6 )
	axis(2,cex.axis=.75,at=c(-(1:2)/4,0:4/4),lab=paste0( c( -(1:2)/4, 0:4/4 )*100, '%' ),tick=F)
	abline( h=0, lty=1, col=1, lwd=1.5 )
	abline( h=1, lty=1, col=1, lwd=1.5 )
	for( x in c( -(1:2)/4, 1:3/4 ) )
	abline( h=x, lty=3, col='lightgrey', lwd=1.5 )

	xvals	<- c( 1.0, 1.3, 1.45, 1.6, 1.9, 2.1, 2.5, 2.7, 2.8, 3.0, 3.1 )
	mylines( xvals[1], h2greml			, se.h2greml			,trunc=F, col=cols[1], points=T,pt.cex=3, asts=sum( p.h2greml < c( 1e-3,1e-2,.05) ) )
	mylines( xvals[2], h2iid[1]			, se.h2iid[1]			,trunc=F, col=cols[2], points=T,pt.cex=3, asts=sum( p.h2hom		< c( 1e-3,1e-2,.05) ) )
	mylines( xvals[3], h2iid[2]			, se.h2iid[2]			,trunc=F, col=cols[3], points=T,pt.cex=3, asts=sum( p.h2het		< c( 1e-3,1e-2,.05) ) )
	mylines( xvals[4], sum(h2iid)		, se.h2iidsum			,trunc=F, col=cols[4], points=T,pt.cex=3, asts=sum( p.h2iid		< c( 1e-3,1e-2,.05) ) )

	for( kk in 1:2 )
	mylines( xvals[4+kk], h2free[kk], se.h2free[kk]	,trunc=F, col=cols[5+kk], points=T,pt.cex=3 )

	for( kk in 1:5 )
	mylines( xvals[6+kk], sig2free[kk], se.sig2free[kk]	,trunc=F, col=cols[4+kk], points=T,pt.cex=3 )

	axis( 1, at=xvals, lab=xlab, las=2, tick=F, cex.axis=1.4 )

	if( !missing(pdfname) )
		dev.off()

}
