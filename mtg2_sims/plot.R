rm( list=ls() )
library(GxEMM)
load( 'Rdata/setup.Rdata' )
load( 'Rdata/adj.Rdata' )
source( '../code/plot_fxns.R' )

modes	<- c( modes, 'mtg2' )
mains	<- c( 'REML (w/ LDAK)', 'PCGC (w/ LDAK)', ' REML (w/ MTG2)' )
names(mains)	<- modes


allps		<- array( NA, dim=c( 3, nx, maxit, 3), dimnames=list( modes, 1:nx, 1:maxit, testnames ) )
h2s			<- array( NA, dim=c( 3, nx, maxit, 5), dimnames=list( modes, 1:nx, 1:maxit, c('g','hom','het','hom0','d') ) )
for( it in 1:maxit )
	for( mode in modes[1:2] )
		for( xval in sample(nx) )
try({

	savefile1	<- paste0( 'Rdata/hom/'	, xs[xval], '_', it, '_', mode, '.Rdata' )
	savefile2	<- paste0( 'Rdata/het/'	, xs[xval], '_', it, '_', mode, '.Rdata' )
	savefile3	<- paste0( 'Rdata/free/'	, xs[xval], '_', it, '_', mode, '.Rdata' )
	savefile3	<- paste0( 'Rdata/iid/'	, xs[xval], '_', it, '_', mode, '.Rdata' )

	load( savefile1 )
	if( out_hom == 'FAILED' ) file.remove( savefile1 )
	h2s		[mode,xval,it,'g']		<- out_hom$h2
	allps	[mode,xval,it,'hom']	<- Waldtest( out_hom$h2, out_hom$h2Covmat[1,1] )

	load( savefile2 )
	if( out_het == 'FAILED' ) file.remove( savefile2 )
	h2s		[mode,xval,it,c('hom','het')]	<- out_het$h2[c('hom','het')]
	allps	[mode,xval,it,'iid']	<- Waldtest( out_het$h2[2], out_het$h2Covmat[2,2] )

	load( savefile3 )
	out_diag	<- out_het1
	if( out_diag == 'FAILED' ) file.remove( savefile3 )
	h2s		[mode,xval,it,c('hom0','d')]	<- out_diag$sig2s[1:2] #/sum( out_diag$sig2s )
	allps	[mode,xval,it,'diag3']	<- MVWaldtest( out_diag$sig2s[2], out_diag$sig2Var[2,2] )

	rm( out_hom, out_het, out_diag )
},silent=T)
print( apply( ! is.na(h2s), c(1,2,4), sum )[1,,c('g','hom','d')] )
print( apply( ! is.na(h2s), c(1,2,4), sum )[2,,c('g','hom','d')] )

mode	<- 'mtg2'
for( it in 1:maxit )
	for( xval in sample(nx) )
try({

	savefile1	<- paste0( 'Rdata/mtg2/', xs[xval], '_', it, '.Rdata' )
	savefile2	<- paste0( 'Rdata/rnm/'	, xs[xval], '_', it, '.Rdata' )
	savefile3	<- paste0( 'Rdata/mrnm/', xs[xval], '_', it, '.Rdata' )

	load( savefile1 )
	h2s		[mode,xval,it,'g']		<- out_mrnm_hom$Vg[1]/( out_mrnm_hom$Vg[1] +  out_mrnm_hom$Ve[1] )
	allps	[mode,xval,it,'hom']	<- Waldtest( out_mrnm_hom$Vg[1], out_mrnm_hom$Vg[2]^2 )

	load( savefile2 )
	h2s		[mode,xval,it,'hom']	<- out_rnm$Vg  [1]/( out_mrnm_hom$Vg[1] + out_rnm$Vgxe[1] + out_mrnm_hom$Ve[1] )
	h2s		[mode,xval,it,'het']	<- out_rnm$Vgxe[1]/( out_mrnm_hom$Vg[1] + out_rnm$Vgxe[1] + out_mrnm_hom$Ve[1] )
	allps	[mode,xval,it,'iid']	<- Waldtest( out_rnm$Vgxe[1], out_rnm$Vgxe[2]^2 )

	#load( savefile3 )
	#print( out_mrnm )
	#h2s		[mode,xval,it,c('d1b','d2b')]	<- out_mrnm$h2
	#allps	[mode,xval,it,'diag3']	<- MVWaldtest( out_diag$sig2s[1+1:2], out_diag$sig2Var[1+1:2,1+1:2] )
},silent=T)
print( apply( ! is.na(h2s), c(1,2,4), sum )[3,,c('g','hom','d')] )

pdf( paste0( '~/figs/gxemm/MTG2.pdf' ), width=16.5, height=6 )

layout( cbind( 1:2, 4+matrix( c(4:6,1:3), 2, 3, byrow=T ), c(3,4) ), widths=c( 1.0, rep(6.3,3), 4.0 ), heights=c(4.5,6.99) )

par( mar=c(.2,0,5.2,.5) )
plot.new()
mtext( side=2, 'Positive Rate', cex=1.2, line=-2.0 )

par( mar=c(.9,0,.9,.5) )
plot.new()
mtext( side=2, 'Heritability'	, cex=1.2, line=-2.0 )

par( mar=c(0,0,0,0) )
plot.new()
legend( 'center', bty='n', title='Genetic Tests\n(nominal p<.05)',
	lwd=c( 3, 3, 5, 3, 3, 3 ),
	col=c( cols1, 2 ),
	lty=c( ltys, 1 ),
	leg=c( testlabs, 'Null' ), cex=1.7, y.intersp=1.1 )

par( mar=c(1.5,0,0,0) )
plot.new()
legend( 'top', bty='n', title='Heritability Estimates',
	lty=c( 1, 1 ,1, 1, 1, 1	),  #  3, 3, 2, 2
	lwd=c( 4, 2 ,3, 3, 3, 3	),  #  3, 3, 2, 2
	pch=c( 16,16,16,16,16,16),  # ,5, 5, 1, 1
	col=c( 'grey', 'yellow4', 6, cols[c(1,2,3)] ),
	leg=expression(h[g]^2, paste( h[hom]^2, ' | Hom E' ), paste( h[het]^2, ' | Hom E' ), h[hom]^2, h[het]^2 ), #paste( h[1], ' Adjusted' ), paste( h[2]^2, ' Global Adjust' ) ),
	cex=1.6, y.intersp=1.5 )
legend( 'bottom', bty='n',
	lty=1:2,
	leg=c('Liability Scale','REML Scale'),
	cex=1.6, y.intersp=1.35, lwd=2.5 )

xlim	<- range(xs)

par( mar=c(5.3,0,.6,.5) )
for( mode in modes ){
	plot( xlim, c(-.1,1.6), type='n', axes=F, xlab='', main='', ylab='', cex.lab=2 )
	box()
	abline( h=0	, col='grey', lty=1, lwd=1 )
	abline( h=.5, col='grey', lty=1, lwd=1 )
	abline( h=1	, col='grey', lty=1, lwd=1 )
	x1s	<- c(0:4/4,1.5)
	if( mode == modes[1]  )
		axis(2,cex.axis=1.1,at=x1s,lab=paste0(x1s*100,'%'))

	if( mode %in% c( 'ldak', 'mtg2' ) ){
		lines1( xs, rowMeans( h2s[mode,,,'g']		* adjs, na.rm=T ), col='grey'		, lty=1, lwd=5, pch=16)
		lines1( xs, rowMeans( h2s[mode,,,'hom0']* adjs, na.rm=T ), col=cols[1]	, lty=1, lwd=3, pch=16)
		lines1( xs, rowMeans( h2s[mode,,,'d']		* adjs, na.rm=T ), col=cols[2]	, lty=1, lwd=3, pch=16)
		lines1( xs, rowMeans( h2s[mode,,,'hom']	* adjs, na.rm=T ), col='yellow4', lty=1, lwd=2, pch=16)
		lines1( xs, rowMeans( h2s[mode,,,'het']	* adjs, na.rm=T ), col=6				, lty=1, lwd=3, pch=1 )

		lines1( xs, rowMeans( h2s[mode,,,'g']					, na.rm=T ), col='grey'		, lty=2, lwd=5, pch=16)
		lines1( xs, rowMeans( h2s[mode,,,'hom0']			, na.rm=T ), col=cols[1]	, lty=2, lwd=3, pch=16)
		lines1( xs, rowMeans( h2s[mode,,,'d']					, na.rm=T ), col=cols[2]	, lty=2, lwd=3, pch=16)
		lines1( xs, rowMeans( h2s[mode,,,'hom']				, na.rm=T ), col='yellow4', lty=2, lwd=2, pch=16)
		lines1( xs, rowMeans( h2s[mode,,,'het']				, na.rm=T ), col=6				, lty=2, lwd=3, pch=1 )

	} else if( mode %in% c( 'pcgc' ) ){
		lines1( xs, rowMeans( h2s[mode,,,'g']					, na.rm=T ), col='grey'		, lty=1, lwd=5, pch=16)
		lines1( xs, rowMeans( h2s[mode,,,'hom0']			, na.rm=T ), col=cols[1]	, lty=1, lwd=3, pch=16)
		lines1( xs, rowMeans( h2s[mode,,,'d']	        , na.rm=T ), col=cols[2]	, lty=1, lwd=5, pch=16, cex=3.2)
		lines1( xs, rowMeans( h2s[mode,,,'hom']	      , na.rm=T ), col='yellow4', lty=1, lwd=2, pch=16)
		lines1( xs, rowMeans( h2s[mode,,,'het']				, na.rm=T ), col=6				, lty=1, lwd=3, pch=16)
	}

	mtext( side=1, expression( beta[z] ), cex=1.35, line=3.8 )
	axis(1,cex.axis=1.1,at=xs,lab=xs)
}

par( mar=c(.6,0,3.3,.5) )
for( mode in modes ){
	plot( xlim, c(0,1), type='n', xlab='', ylab='', axes=F )
	box()
	mtext( side=3, mains[mode]			, cex=1.65, line=1.0 )
	x1s	<- 0:2/2
	if( mode == modes[1] )
		axis(2,cex.axis=1.1,at=x1s,lab=paste0(x1s*100,'%'))

	abline( h=.05, col=2, lty=1, lwd=2 )
	for( ii in c( 'iid', 'diag3', 'hom' ) )
		lines1( xs, apply(allps[mode,,,ii]< .05,1,mean,na.rm=T), lty=ltys[ii], col=cols1[ii] )#, y.max=10
}

dev.off()
