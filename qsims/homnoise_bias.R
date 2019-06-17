rm( list=ls() )
library(GxEMM)
load( 'Rdata/setup.Rdata' )
load( 'Rdata/preprocess.Rdata' ) # ahat, bhat
#ahat	<- 0.001080078
#bhat	<- 0.088220372
cnst	<- 1/(1+ahat+bhat)

sigtype	<- 5
ws0		<- sapply( 1:nx, function(ix) all_tauhets(ix)[[sigtype]][1] )
sig2	<- mean( ws0 )

source( '../code/plot_fxns.R' )

allps		<- array( NA, dim=c( 5, nx, maxit, length(testnames)), dimnames=list( 1:5, 1:nx, 1:maxit, testnames ) )
h2s			<- array( NA, dim=c( 5, nx, maxit, 9								), dimnames=list( 1:5, 1:nx, 1:maxit, c('g','hom','het','hom0','d1','d2','hom0b','d1b','d2b') ) )


for( xval in 1:nx )
	for( it in 1:maxit )
try({

	savefile4	<- paste0( 'Rdata_HE/diag1/', sigtype, '_', xval, '_', it, '.Rdata' )
	#savefile4	<- paste0( 'Rdata_HE/diag1/', sigtype, '_', xval, '_', it, '.Rdata' )

	load( savefile4 )
	if( out_diag1 == 'FAILED' ) file.remove( savefile4 )
	h2s		[sigtype,xval,it,c('d1','d2')]		<- out_diag1$sig2s[2:3]
	h2s		[sigtype,xval,it,c('d1b','d2b')]	<- out_diag1$h2

	rm(  out_diag1 )
},silent=T)
print( apply( ! is.na(h2s[,,,c('g','hom','d1','d1b')]), c(1,2,4), sum ) )

pdf( '~/figs/gxemm/Fig7.pdf', width=12, height=6 )
par( mfrow=c(1,2) )
par( mar=c(4.5,7,4.5,.5) )

plot( xlim, c(-.8,.8), type='n', axes=F, xlab=xlabs[sigtype], main='Environment-specific variance components\nignoring noise heterogeneity', ylab=expression( hat(v)[k] ), cex.lab=1.9 )
box()
axis(1,cex.axis=1.1,at=xs,lab=round( xs*2.75,2 ))
axis(2,cex.axis=1.1,at=(0:8-4)/5)

lines( xs, 0*xs, col=cols[3], lwd=3, lty=2 )
lines( xs, 0*xs, col=cols[4], lwd=3, lty=2 )

lines( xs, sapply( ws0			- sig2, function(x) cnst*x ), col=cols[3], lwd=3 )
lines( xs, sapply( rev(ws0)	- sig2, function(x) cnst*x ), col=cols[4], lwd=3 )

for( j in 1:nx ){
	pointline( xs[j], h2s[sigtype,j,,'d1']	, col=cols[3], cex=3 )
	pointline( xs[j], h2s[sigtype,j,,'d2']	, col=cols[4], cex=3 )
}

plot( xlim, c(-1,.8), type='n', axes=F, xlab=xlabs[sigtype], main='Environment-specific heritabilities\nignoring noise heterogeneity', ylab=expression( hat(h)[k]^2 ), cex.lab=2 )
box()
axis(1,cex.axis=1.1,at=xs,lab=round( xs*2.75,2 ))
axis(2,cex.axis=1.1,at=(0:9-5)/5,lab=paste0((0:9-5)/5*100,'%'))


lines( xs, sapply( 2.75*(   xs), function(x) .1 /(.1+x+(1-.1-sum(ws*c(x,(2.75*.6-x))))) ), col=cols[3], lwd=3, lty=2 )
lines( xs, sapply( 2.75*(.6-xs), function(x) .1 /(.1+x+(1-.1-sum(ws*c(x,(2.75*.6-x))))) ), col=cols[4], lwd=3, lty=2 )

lines( xs, sapply( ws0			- sig2, function(x) (.1+cnst*x) /(.1+cnst*x+sig2) ), col=cols[3], lwd=3 )
lines( xs, sapply( rev(ws0)	- sig2, function(x) (.1+cnst*x) /(.1+cnst*x+sig2) ), col=cols[4], lwd=3 )

for( j in 1:nx ){
	pointline( xs[j], h2s[sigtype,j,,'d1b']	, col=cols[3], cex=3 )
	pointline( xs[j], h2s[sigtype,j,,'d2b']	, col=cols[4], cex=3 )
}

dev.off()
