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
		text( x, y+.2, lab=c( '', '*', '**', '***' )[asts+1], cex=1.5 )

}

gxemm.plot	<- function(
	gxemm.hom, gxemm.iid, gxemm.free,
	pdfname,
	cols=c( 'grey', 1, 'gold', 'chocolate3', 'yellow4', '#FB6542', '#3F681C', '#FB6542', '#3F681C' ),
	Enames=1:P,
	main='',
	ylim=c(-.1,1.1)
){

	p.h2g	<- Waldtest( gxemm.hom$h2, gxemm.hom$h2Covmat[1,1] )
	h2g		<- gxemm.hom$h2
	se.h2g<- sqrt(gxemm.hom$h2Covmat[1,1])

	p.h2hom	<- Waldtest( gxemm.iid$h2[1], gxemm.iid$h2Covmat[1,1] )
	p.h2het	<- Waldtest( gxemm.iid$h2[2], gxemm.iid$h2Covmat[2,2] )
	p.h2iid	<- Waldtest( sum(gxemm.iid$h2), sum(gxemm.iid$h2Covmat) )
	h2iid			<- gxemm.iid$h2
	se.h2iid	<- sqrt(diag(gxemm.iid$h2Covmat))
	se.h2iidsum	<- sqrt(sum(gxemm.iid$h2Covmat))

	h2free		<- gxemm.free$h2
	se.h2free	<- sqrt(diag(gxemm.free$h2Covmat))

	P		<- length( h2free )
	Pe	<- ifelse( FALSE, 1+P+P+1, 1+P+P )
	sig2free		<- c(
		gxemm.free$sig2s[1:(P+1)],
		c(gxemm.free$sig2s[(1+P)+1:(P-1)],0)+gxemm.free$sig2s[Pe]
	)
	se.sig2free	<- sqrt(c(
		diag(gxemm.free$sig2Var)[1:(P+1)],
		sapply( 1:(P-1), function(p) sum(gxemm.free$sig2Var[c(1+P+p,Pe),c(1+P+p,Pe)]) ),
		gxemm.free$sig2Var[Pe,Pe]
	))
	p.h2free	<- sapply( 1:P	, function(p) Waldtest( h2free	[p], se.h2free	[p]^2 ) )
	p.free		<- sapply( 1:Pe	, function(p) Waldtest( sig2free[p], se.sig2free[p]^2 ) )

	if( !missing(pdfname) )
		pdf( pdfname, width=6.8, height=5 )

	layout( matrix( c(1,2,1,3), 2, 2 ), width=c(7,5), height=c(1,12) )
	par( mar=rep(0,4) )
	plot.new()
	mtext( side=1, line=-1, text=main, cex=1.5 )

	par( mar=c(5,4.5,.5,1.5) )

	xvals	<- c( 1, 1.6, 2, 2.4, 3.5+seq(-1,1,len=P)*P/10 )

	plot( range(xvals)+c(-.1,.1), ylim, type='n', axes=F, xlab='', ylab='Heritability Estimates', bty='n', cex.lab=1.3 )
	axis(2,cex.axis=.75,at=c(-(1:2)/4,0:4/4),lab=paste0( c( -(1:2)/4, 0:4/4 )*100, '%' ),tick=F)
	abline( h=0, lty=1, col=1, lwd=1.5 )
	abline( h=1, lty=1, col=1, lwd=1.5 )
	for( x in c( -(1:2)/4, 1:3/4 ) )
	abline( h=x, lty=3, col='lightgrey', lwd=1.5 )

	mylines( xvals[1], h2g				, se.h2g			,trunc=F, col=cols[1], points=T,pt.cex=3, asts=sum( p.h2g		< c( 1e-3,1e-2,.05) ) )
	mylines( xvals[5], h2iid[1]		, se.h2iid[1]	,trunc=F, col=cols[2], points=T,pt.cex=3, asts=sum( p.h2hom	< c( 1e-3,1e-2,.05) ) )
	mylines( xvals[6], h2iid[2]		, se.h2iid[2]	,trunc=F, col=cols[3], points=T,pt.cex=3, asts=sum( p.h2het	< c( 1e-3,1e-2,.05) ) )
	mylines( xvals[2], sum(h2iid)	, se.h2iidsum	,trunc=F, col=cols[4], points=T,pt.cex=3, asts=sum( p.h2iid	< c( 1e-3,1e-2,.05) ) )

	for( p in 1:P )
	mylines( xvals[p+2], h2free[p], se.h2free[p]	,trunc=F, col=cols[5+p], points=T,pt.cex=3, asts=sum( p.h2free[p]	< c( 1e-3,1e-2,.05) ) )

	xlab	<- c( expression( h['g']^2, h['IID']^2 )	, sapply( Enames, function(p) substitute(h[pp]^2, list(pp=p)) ), expression( h['hom']^2, h['het']^2 ) )
	axis( 1, at=xvals, lab=xlab, las=2, tick=F, cex.axis=1.4 )

	par( mar=c(5,4.5,.5,.5) )
	xvals	<- c( 1.2, 2+seq(-1,1,len=P)*P/10, 3+seq(-1,1,len=P)*P/10 )

	plot( range(xvals)+c(-.1,.1), range(c(sig2free,-.1))*1.3, type='n', axes=F, xlab='', ylab='Free Model Variance Components', bty='n', cex.lab=1.3 )
	axis(2,cex.axis=.75,tick=F)
	abline( h=0, lty=1, col=1, lwd=1.5 )

	for( p in 1:(1+2*P) )
	mylines( xvals[p], sig2free[p], se.sig2free[p]	,trunc=F, col=cols[4+p], points=T,pt.cex=3, asts=sum( p.free[p]	< c( 1e-3,1e-2,.05) ) )

	xlab	<- c( expression( sigma['hom']^2 ), sapply( Enames, function(p) substitute(v[pp], list(pp=p)) ), sapply( Enames	, function(p) substitute(w[pp], list(pp=p)) ) )
	axis( 1, at=xvals, lab=xlab, las=2, tick=F, cex.axis=1.4 )

	if( !missing(pdfname) )
		dev.off()

}
