rm( list=ls() )
load( 'Rdata/setup.Rdata' )
source( '../code/plot_fxns.R' )

#for( suffix in c( '', '_bin', '_bina', '_gcta', '_pcgcsmart' ) ){
for( suffix in c( '_bin', '_bina', '_pcgcsmart' ) ){

xlabs			<- expression( h[hom]^2, h[het]^2, Delta, omega, tau )

h2s			<- array( NA, dim=c( 5, nx, maxit, 6	), dimnames=list( 1:5, 1:nx, 1:maxit, c('g','hom','het','hom0','d1','d2') ) )
h2.ses	<- array( NA, dim=c( 5, nx, maxit, 6	), dimnames=list( 1:5, 1:nx, 1:maxit, c('g','hom','het','hom0','d1','d2') ) )
for( sigtype in 1:5 )
	for( xval in 1:nx )
		for( it in 1:maxit )
try({

	if( suffix == '' ){
		savefile1	<- paste0( 'Rdata/hom/'	, sigtype, '_', xval, '_', it, '.Rdata' )
		savefile2	<- paste0( 'Rdata/het/'	, sigtype, '_', xval, '_', it, '.Rdata' )
		savefile4	<- paste0( 'Rdata/diag1/',sigtype, '_', xval, '_', it, '.Rdata' )
	} else if( suffix == '_gcta' ){
		savefile1	<- paste0( '../final_sim/Rdata/hom/'	, sigtype, '_', xval, '_', it, '.Rdata' )
		savefile2	<- paste0( '../final_sim/Rdata/het/'	, sigtype, '_', xval, '_', it, '.Rdata' )
		savefile4	<- paste0( '../final_sim/Rdata/diag1/',sigtype, '_', xval, '_', it, '.Rdata' )
	} else if( suffix %in% c( '_bin', '_bina' ) ){
		savefile1	<- paste0( '../final_sim_bin/Rdata/hom/'	, xs[xval], '_', it, '_', sigtype, '_pcgcsmart', '.Rdata' )
		savefile2	<- paste0( '../final_sim_bin/Rdata/het/'	, xs[xval], '_', it, '_', sigtype, '_pcgcsmart', '.Rdata' )
		#savefile4	<- paste0( '../final_sim_bin/Rdata/diag1/', xs[xval], '_', it, '_', sigtype, '_pcgcsmart', '.Rdata' )
		savefile4	<- paste0( '../final_sim_bin/Rdata/diag/', xs[xval], '_', it, '_', sigtype, '_pcgcsmart', '.Rdata' )
	} else {
		savefile1	<- paste0( 'Rdata/hom/'	, xs[xval], '_', it, '_', sigtype, '_pcgcsmart', '.Rdata' )
		savefile2	<- paste0( 'Rdata/het/'	, xs[xval], '_', it, '_', sigtype, '_pcgcsmart', '.Rdata' )
		savefile4	<- paste0( 'Rdata/diag1/',xs[xval], '_', it, '_', sigtype, '_pcgcsmart', '.Rdata' )
	}

	load( savefile1 )
	if( suffix == '_bina' ){
	h2s		[sigtype,xval,it,'g']		<- out_hom$sig2s[1]
	h2.ses[sigtype,xval,it,'g']		<- sqrt( diag( out_hom$sig2Var )[1] )
	} else {
	h2s		[sigtype,xval,it,'g']		<- out_hom$h2
	h2.ses[sigtype,xval,it,'g']		<- out_hom$h2.se
	}

	load( savefile2 )
	if( suffix == '_bina' ){
	h2s		[sigtype,xval,it,c('hom','het')]	<- out_het$sig2s[1:2]
	h2.ses[sigtype,xval,it,c('hom','het')]	<- sqrt( diag( out_het$sig2Var )[1:2] )
	} else {
	h2s		[sigtype,xval,it,c('hom','het')]	<- out_het$h2[c('hom','het')]
	h2.ses[sigtype,xval,it,c('hom','het')]	<- out_het$h2.se[c('hom','het')]
	}

	load( savefile4 )
	if( suffix == '_bin' ){
		h2s		[sigtype,xval,it,c('hom0','d1','d2')]	<- out_diag$h2
		h2.ses[sigtype,xval,it,c('hom0','d1','d2')] <- out_diag$h2.se
		rm( out_diag )
	} else if( suffix == '_bina' ){
		h2s		[sigtype,xval,it,c('hom0','d1','d2')]	<- out_diag$sig2s[1:3]
		h2.ses[sigtype,xval,it,c('hom0','d1','d2')] <- sqrt( diag( out_diag$sig2Var )[1:3] )
		rm( out_diag )

	} else {
		h2s		[sigtype,xval,it,c('hom0','d1','d2')]	<- out_diag1$h2
		h2.ses[sigtype,xval,it,c('hom0','d1','d2')]<- out_diag1$h2.se
	}

	rm( out_hom, out_het, out_diag1 )
},silent=T)

h2.ses0	<- apply( h2s, c(1,2,4), sd, na.rm=T )
#h2.ses0[,,'het']	<- apply( h2s[,,,'hom'] + h2s[,,,'het'], c(1,2), sd, na.rm=T )

pdf( paste0( '~/figs/gxemm/final_se', suffix, '.pdf' ), width=13.5, height=4.0 )
layout( matrix( c(1, 3:7, 2 ), nrow=1 ), widths=c( 2.0, rep(5,5), 4.8 ) )
par( mar=c(5.5,0,4.2,.5) )

plot.new()
mtext( side=2, 'Heritability S.E.'	, cex=1.39, line=-3.2 )

par( mar=c(2.5,0,4.2,.5) )
plot.new()
legend( 'top', bty='n', fill=c( 'grey', cols[1:4] ), leg=expression(h[g]^2,h[hom]^2,h[het]^2, h[1]^2, h[2]^2 ), cex=1.8, y.intersp=1.4 )
legend( 'bottom', bty='n', lty=c(1,NA), pch=c(NA,16), leg=c( 'SD over sims', 'SEs per sims' ), cex=1.8 )


par( mar=c(5.5,0,4.2,.5) )
let.i	<- 5
for( sigtype in c(1:3,5,4) ){

	plot( xlim, c(0,ifelse( suffix %in% c( '_bin', '_bina' ), .7 ,.2 )), type='n', axes=F, xlab='', main='', ylab='', cex.lab=2 )
	box()

	text(.01,.684, letters[ let.i	<- let.i+1 ],cex=2.2)
	mtext( side=1, xlabs[sigtype]			, cex=1.35, line=4.1 )
	mtext( side=3, mains[sigtype]			, cex=1.65, line=1.7 )

	axis(1,cex.axis=1.1,at=xs,lab=xs)

	x1s	<- 0:5/5
	if( sigtype == 1 )
		axis(2,cex.axis=1.1,at=0:4/4*.2)

	if( sigtype == 1 ){
		lines( xs, h2.ses0[sigtype,,'g'],col='grey', lwd=3 )
	} else if( sigtype == 2 ){
		lines( xs, h2.ses0[sigtype,,'hom'],col=cols[1], lwd=3 )
		lines( xs, h2.ses0[sigtype,,'het'],col=cols[2], lwd=3 )
	} else if( sigtype %in% 3:5 ){
		lines( xs, h2.ses0[sigtype,,'d1'],col=cols[3], lwd=3 )
		lines( xs, h2.ses0[sigtype,,'d2'],col=cols[4], lwd=3 )
	}

	for( j in 1:nx )
		if( sigtype %in% 1 ){
			pointline( xs[j], h2.ses[sigtype,j,,'g']		, col='grey'	, cex=3 )
		} else if( sigtype %in% 2 ){
			pointline( xs[j], h2.ses[sigtype,j,,'hom']	, col=cols[1]	, cex=3 )
			pointline( xs[j], h2.ses[sigtype,j,,'het']	, col=cols[2]	, cex=3 )
		} else if( sigtype %in% 3:5 ){
			pointline( xs[j], h2.ses[sigtype,j,,'d1']	, col=cols[3], cex=3 )
			pointline( xs[j], h2.ses[sigtype,j,,'d2']	, col=cols[4], cex=3 )
		}

}
dev.off()
}
