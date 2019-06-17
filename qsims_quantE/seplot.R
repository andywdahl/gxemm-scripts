rm( list=ls() )
load( 'Rdata/setup.Rdata' )
source( '../code/plot_fxns.R' )

h2s			<- array( NA, dim=c( 5, nx, maxit, 6	), dimnames=list( 1:5, 1:nx, 1:maxit, c('g','hom','het','hom0','d1','d2') ) )
h2.ses	<- array( NA, dim=c( 5, nx, maxit, 6	), dimnames=list( 1:5, 1:nx, 1:maxit, c('g','hom','het','hom0','d1','d2') ) )
for( sigtype in 1:5 )
	for( xval in 1:nx )
		for( it in 1:maxit )
try({

	savefile1	<- paste0( 'Rdata/hom/'	, sigtype, '_', xval, '_', it, '.Rdata' )
	savefile2	<- paste0( 'Rdata/het/'	, sigtype, '_', xval, '_', it, '.Rdata' )
	savefile3	<- paste0( 'Rdata/diag/', sigtype, '_', xval, '_', it, '.Rdata' )

	load( savefile1 )
	h2s		[sigtype,xval,it,'g']						<- out_hom$h2
	h2.ses[sigtype,xval,it,'g']						<- sqrt(out_hom$h2Covmat[1,1])
	rm( out_hom )

	load( savefile2 )
	h2s		[sigtype,xval,it,c('hom','het')]<- out_het$h2
	h2.ses[sigtype,xval,it,c('hom','het')]<- sqrt(diag(out_het$h2Covmat))
	rm( out_het )

	load( savefile3 )
	h2s		[sigtype,xval,it,c('d1','d2')]	<- out_diag$h2
	h2.ses[sigtype,xval,it,c('d1','d2')]	<- sqrt(diag(out_diag$h2Covmat))
	rm( out_diag )

},silent=T)

h2.ses0	<- apply( h2s, c(1,2,4), sd, na.rm=T )
rm( h2s )

pdf( '~/figs/gxemm/Fig14c.pdf', width=15.5, height=4.2 )
layout( matrix( c(1, 3:7, 2 ), nrow=1 ), widths=c( 1.5, rep(5,5), 4.3 ) )
par( mar=c(5.5,0,4.2,.5) )

plot.new()
mtext( side=2, 'Heritability S.E.'	, cex=1.39, line=-2.1 )

par( mar=c(1.2,0,1.2,.1) )
plot.new()
legend( 'top', bty='n', fill=c( 'grey', cols[1:4] ), leg=expression(h[g]^2,h[hom]^2,h[het]^2, h[1]^2, h[2]^2 ), cex=1.8, y.intersp=1.4 )
legend( 'bottom', bty='n', lty=c(1,NA), pch=c(NA,16), leg=c( 'SD over sims', 'SEs per sims' ), cex=1.8 )


par( mar=c(5.5,0,4.2,.5) )
let.i	<- 0
for( sigtype in c(1:3,5,4) ){

	if( sigtype == 4 ){
		plot( c(.08,.62), c(0,.22), type='n', axes=F, xlab='', main='', ylab='', cex.lab=2 )
	} else {
		plot( xlim			, c(0,.22), type='n', axes=F, xlab='', main='', ylab='', cex.lab=2 )
	}
	box()

	text(.01,.13, letters[ let.i	<- let.i+1 ],cex=2.2)
	mtext( side=1, xlabs[sigtype]			, cex=1.35, line=4.1 )
	mtext( side=3, mains[sigtype]			, cex=1.35, line=1.7 )

	if( sigtype == 4 ){
	print( xs )
print( h2.ses0[sigtype,,'g'] )
		axis(1,cex.axis=1.1,at=xs[-1],lab=round(xs[-1]/.3,2))
	} else {
		axis(1,cex.axis=1.1,at=xs,lab=xs)
	}
	if( sigtype == 1 )
		axis(2,cex.axis=1.1,at=0:4/4*.2)

	lines( xs					, h2.ses0[sigtype,,'g']		, col='grey'	, lwd=3 )
	lines( xs					, h2.ses0[sigtype,,'hom']	, col=cols[1]	, lwd=3 )
	lines( xs					, h2.ses0[sigtype,,'het']	, col=cols[2]	, lwd=3 )
	lines( xs					, h2.ses0[sigtype,,'d1']	, col=cols[3]	, lwd=3 )
	lines( xs					, h2.ses0[sigtype,,'d2']	, col=cols[4]	, lwd=3 )

	for( j in 1:nx ){
		pointline( xs[j], h2.ses[sigtype,j,,'g']	, col='grey'	, cex=3 )
		pointline( xs[j], h2.ses[sigtype,j,,'hom'], col=cols[1]	, cex=3 )
		pointline( xs[j], h2.ses[sigtype,j,,'het'], col=cols[2]	, cex=3 )
		pointline( xs[j], h2.ses[sigtype,j,,'d1']	, col=cols[3]	, cex=3 )
		pointline( xs[j], h2.ses[sigtype,j,,'d2']	, col=cols[4]	, cex=3 )
	}

}
dev.off()
