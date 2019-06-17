rm( list=ls() )
library(GxEMM)
load( 'Rdata/setup.Rdata' )
source( '../code/plot_fxns.R' )

testmodes	<- c( 'wald' )
pdfs			<- c( '~/figs/gxemm/Fig12.pdf' )
legsubs		<- list(
	c(1,2,3,4,5,6,7,8),
	c(3,4,5,8)
)

names(pdfs)			<- testmodes
names(legsubs)	<- testmodes

for( testmode in testmodes ){

allps		<- array( NA, dim=c( 5, nx, maxit, length(testnames)), dimnames=list( 1:5, 1:nx, 1:maxit, testnames ) )
h2s			<- array( NA, dim=c( 5, nx, maxit, 9								), dimnames=list( 1:5, 1:nx, 1:maxit, c('g','hom','het','hom0','d1','d2','hom0b','d1b','d2b') ) )
for( sigtype in 1:5 )
	for( xval in 1:nx )
		for( it in 1:maxit )
try({

	savefile1	<- paste0( 'Rdata/hom/'	, sigtype, '_', xval, '_', it, '.Rdata' )
	savefile2	<- paste0( 'Rdata/het/'	, sigtype, '_', xval, '_', it, '.Rdata' )
	savefile3	<- paste0( 'Rdata/diag/', sigtype, '_', xval, '_', it, '.Rdata' )

	load( savefile1 )
	if( out_hom == 'FAILED' ) file.remove( savefile1 )
	h2s		[sigtype,xval,it,'g']		<- out_hom$h2
	allps	[sigtype,xval,it,'hom']	<- Waldtest( out_hom$h2, out_hom$h2Covmat[1,1] )

	load( savefile2 )
	if( out_het == 'FAILED' ) file.remove( savefile2 )
	h2s		[sigtype,xval,it,c('hom','het')]	<- out_het$h2
	allps	[sigtype,xval,it,'iid']	<- Waldtest( out_het$h2[2], out_het$h2Covmat[2,2] )

	load( savefile3 )
	if( out_diag == 'FAILED' ) file.remove( savefile3 )
	h2s		[sigtype,xval,it,c('d1','d2')]	<- out_diag$h2
	#allps	[sigtype,xval,it,'h2eq']	<- h2_equal_test( out_diag$h2, out_diag$h2Covmat )

	#Tmat	<- cbind( c(1,-1,0), c(0,0,1) )
	#allps	[sigtype,xval,it,'diag1']	<- MVWaldtest( t(Tmat) %*% out_diag$sig2s[2:4], t(Tmat) %*% out_diag$sig2Var[2:4,2:4] %*% Tmat )
	allps	[sigtype,xval,it,'diag2']	<- MVWaldtest( out_diag$sig2s[2:3], out_diag$sig2Var[2:3,2:3] )

	rm( out_hom, out_het, out_diag, out_diag1 )
},silent=T)
print( apply( ! is.na(h2s[,,,c('g','hom','d1')]), c(1,2,4), sum ) )

pdf( pdfs[testmode], width=16.7, height=7.4 )
layout( cbind( c(1,2), 4+matrix( c(6:10,1:5), 2, 5, byrow=T ), 3:4 ), widths=c( 2.0, rep(5,5), 4.7 ), heights=c(4.2,5.99) )

par( mar=c(1.2,0,4.2,.5) )
plot.new()
mtext( side=2, 'Positive Rate\n(nominal p<.05)', cex=1.39, line=-4.5 )

par( mar=c(5.5,0,.4,.5) )
plot.new()
mtext( side=2, 'Heritability Estimates'	, cex=1.39, line=-3.9 )

par( mar=c(0,1,1.5,0) )
plot.new()
legend( 'left', bty='n',
	lty=c( ltys, 1								)[legsubs[[testmode]]],
	lwd=rep( 3 , length(cols1)+1	)[legsubs[[testmode]]],
	pch=rep( 16, length(cols1)+1 	)[legsubs[[testmode]]],
	col=c( cols1, 2								)[legsubs[[testmode]]],
	leg=c( testlabs, 'Null'				)[legsubs[[testmode]]],
	text.width=strwidth(testlabs[6]),
	cex=1.7, y.intersp=1.1 )

plot.new()
legend( 'topleft', bty='n',
	lty=c( rep( 1, 5 ), 3, 3 ),
	pch=c( rep( 16, 5 ),5, 5 ),
	col=c( 'grey', cols[c(1:4,3:4)] ),
	leg=expression(h[g]^2,h[hom]^2,h[het]^2, h[1]^2, h[2]^2, paste( h[1]^2, '  | Hom Noise' ), paste( h[2]^2, '  | Hom Noise' ) ),
	cex=1.7, y.intersp=1.5, lwd=3 )

par( mar=c(5.5,0,.4,.5) )
let.i	<- 5
for( sigtype in c(1:3,5,4) ){

	if( sigtype == 4 ){
	plot( c(.08,.62), c(-.1,.99), type='n', axes=F, xlab='', main='', ylab='', cex.lab=2 )
	} else {
	plot( xlim, c(-.1,.99), type='n', axes=F, xlab='', main='', ylab='', cex.lab=2 )
	}
	box()

	text(ifelse(sigtype==4,xs[2]+.01,.01),ifelse(sigtype==5,.76,.94), letters[ let.i	<- let.i+1 ],cex=2.4)
	mtext( side=1, xlabs[sigtype]			, cex=1.35, line=4.1 )

	if( sigtype == 4 ){
		axis(1,cex.axis=1.1,at=xs[-1],lab=round(xs[-1]/.3,2))
	} else if( sigtype == 3 ){
		axis(1,cex.axis=1.1,at=xs,lab=round( xs*10/6,2 ))
	} else if( sigtype == 5 ){                    
		axis(1,cex.axis=1.1,at=xs,lab=round( xs*2.75,2 ))
	} else {                                      
		axis(1,cex.axis=1.1,at=xs,lab=xs)
	}

	if( sigtype == 1 )
		axis(2,cex.axis=1.1,at=0:5/5,lab=paste0(0:5*20,'%'))

	if( sigtype %in% 4:5 ){
		lines1( xs, rowMeans(h2s[sigtype,,,'d1b'],na.rm=T)	, col=cols[3]	, pch=5, cex=3.3, lty=3, lwd=1.6 )
		lines1( xs, rowMeans(h2s[sigtype,,,'d2b'],na.rm=T)	, col=cols[4]	, pch=5, cex=3.3, lty=3, lwd=1.6 )
	}

	if( sigtype %in% 1:2 ){
		lines( xs, xs*(sigtype==1),col=cols[1], lwd=3 )
		lines( xs, xs*(sigtype==2),col=cols[2], lwd=3 )
	} else if( sigtype %in% 3 ){
		lines( xs, sapply(   10/6*xs, function(x) x/(x+(1-sum(ws*c(x,(1-x))))) ), col=cols[4], lwd=3 )
		lines( xs, sapply( 1-10/6*xs, function(x) x/(x+(1-sum(ws*c(x,(1-x))))) ), col=cols[3], lwd=3 )
	} else if( sigtype %in% 5 ){
		lines( xs, sapply( 2.75*(   xs), function(x) .1 /(.1+x+(1-.1-sum(ws*c(x,(2.75*.6-x))))) ), col=cols[3], lwd=3 )
		lines( xs, sapply( 2.75*(.6-xs), function(x) .1 /(.1+x+(1-.1-sum(ws*c(x,(2.75*.6-x))))) ), col=cols[4], lwd=3 )
	} else if( sigtype %in% 4 ){
		lines( xs, xs*0+.1,col=cols[3], lwd=3 )
		lines( xs, xs*0+.1,col=cols[4], lwd=3, lty=2 )
	}

	for( j in 1:nx ){
		if( sigtype %in% 1:2 ){
			pointline( xs[j], h2s[sigtype,j,,'hom']	, col=cols[1]	, cex=3 )
			pointline( xs[j], h2s[sigtype,j,,'het']	, col=cols[2]	, cex=3 )
		}
		if( sigtype %in% 1:5 )
			pointline( xs[j], h2s[sigtype,j,,'g']		, col='grey'	, cex=3 )
		if( sigtype %in% 3:5 ){
			pointline( xs[j], h2s[sigtype,j,,'d1']	, col=cols[3], cex=3 )
			pointline( xs[j], h2s[sigtype,j,,'d2']	, col=cols[4], cex=3 )
		}
	}
}

par( mar=c(1.2,0,4.2,.5) )

let.i	<- 0
for( sigtype in c(1:3,5,4) ){

	if( sigtype == 4 ){
	plot( c(.08,.62), c(0,1), type='n', xlab='', ylab='', axes=F )
	} else {
	plot( xlim, c(0,1), type='n', xlab='', ylab='', axes=F )
	}

	box()
	mtext( side=3, mains[sigtype]			, cex=1.45, line=1.7 )
	if( sigtype == 1 )
	axis(2,cex.axis=1.2,at=0:4/4)
	text(ifelse(sigtype==4,xs[2]+.01,.01),ifelse(sigtype==5,.5,.90),letters[ let.i	<- let.i+1 ],cex=2.4)

	lines(	xs, rep(.05,length(xs))				, col=2, lty=1, lwd=3	)
	for( ii in c( 'iid', 'hom', 'diag2' ) )
		lines1( xs, apply(allps[sigtype,,,ii] < .05,1,mean,na.rm=T), col=cols1[ii], y.max=1.01, lty=ltys[ii] )
}
dev.off()
}
source( 'seplot.R' )
