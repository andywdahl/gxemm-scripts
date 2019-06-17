rm( list=ls() )
library(GxEMM)
load( 'Rdata/setup.Rdata' )
source( '../code/plot_fxns.R' )

dirs	<- c( '../mu_binsim/', '../asc_binsim/' )

h2s			<- array( NA, dim=c( 2, 2, nx, maxit, 7), dimnames=list( dirs, modes, 1:nx, 1:maxit, c('g','hom','het','d1','d2','d1b','d2b') ) )
h2.ses	<- array( NA, dim=c( 2, 2, nx, maxit, 7), dimnames=list( dirs, modes, 1:nx, 1:maxit, c('g','hom','het','d1','d2','d1b','d2b') ) )
for( it in 1:maxit )
	for( mode in modes )
		for( xval in sample(nx) )
			for( dir in dirs )
try({

if( dir == dirs[1] ){
	fs	<- xs
} else {
	fs		<- round( seq( .01, .5, len=nx ), 2 )
}

	savefile1	<- paste0( dir, 'Rdata/hom/'	, fs[xval], '_', it, '_', mode, '.Rdata' )
	savefile2	<- paste0( dir, 'Rdata/het/'	, fs[xval], '_', it, '_', mode, '.Rdata' )
	savefile3	<- paste0( dir, 'Rdata/diag/'	, fs[xval], '_', it, '_', mode, '.Rdata' )
	savefile4	<- paste0( dir, 'Rdata/diag1/', fs[xval], '_', it, '_', mode, '.Rdata' )
	savefile5	<- paste0( dir, 'Rdata/diag2/', fs[xval], '_', it, '_', mode, '.Rdata' )

	load( savefile1 )
	if( out_hom == 'FAILED' ) file.remove( savefile1 )
	h2s		[dir,mode,xval,it,'g']		<- out_hom$h2
	h2.ses[dir,mode,xval,it,'g']		<- sqrt(diag(out_hom$h2Covmat))

	load( savefile2 )
	if( out_het == 'FAILED' ) file.remove( savefile2 )
	h2s		[dir,mode,xval,it,c('hom','het')]	<- out_het$h2[c('hom','het')]
	h2.ses[dir,mode,xval,it,c('hom','het')]	<-sqrt(diag(out_het$h2Covmat)[1:2])

	load( savefile3 )
	if( out_diag == 'FAILED' ) file.remove( savefile3 )
	h2s		[dir,mode,xval,it,c('d1b','d2b')]	<- out_diag$h2
	h2.ses[dir,mode,xval,it,c('d1b','d2b')]	<- sqrt(diag(out_diag$h2Covmat))

	load( savefile4 )
	if( out_diag1 == 'FAILED' ) file.remove( savefile4 )
	h2s		[dir,mode,xval,it,c('d1','d2')]	<- out_diag1$h2
	h2.ses[dir,mode,xval,it,c('d1','d2')]	<- sqrt(diag(out_diag1$h2Covmat))
	rm( out_hom, out_het, out_diag, out_diag1 )
},silent=T)

print( apply( ! is.na(h2s), c(1,2,3,5), sum )[1,,,c('g','hom','d1','d1b')] )
print( apply( ! is.na(h2s), c(1,2,3,5), sum )[2,,,c('g','hom','d1','d1b')] )
print( apply( ! is.na(h2.ses), c(1,2,3,5), sum )[1,,,c('g','hom','d1','d1b')] )
print( apply( ! is.na(h2.ses), c(1,2,3,5), sum )[2,,,c('g','hom','d1','d1b')] )

pdf( paste0( '~/figs/gxemm/Fig3_se.pdf' ), width=16.5, height=8.2 )

layout( cbind( 1:2, 4+matrix( c(1:2,3:4), 2, 2, byrow=T ), rep(4,2), 8+matrix( c(1:2,3:4), 2, 2, byrow=T ), c(3,3) ), widths=c( 1.3, rep(4.3,2), .4, rep(4.3,2), 4.0 ), heights=c(14,5) )

par( mar=c(.9,0,.9,.5) )
plot.new()
mtext( side=2, 'Standard Errors'		, cex=1.2, line=-4.0 )

par( mar=c(5.5,0,.4,.5) )
plot.new()
mtext( side=2, 'Disease\nPrevalence', cex=1.2, line=-4.0 )

par( mar=c(5,0,5,0) )
plot.new()
legend( 'top', bty='n',
	lty=c( 1, 1 ,1, 1, 3, 3, 2, 2 ),
	pch=c( 16,16,16,16,5, 5, 1, 1 ),
	col=c( cols[c(1,2,3,4,3:4)], 6, 3 ),
	leg=expression(h[g]^2, h[het]^2, h[1]^2, h[2]^2, paste( h[1]^2, ' | Hom E' ), paste( h[2]^2, ' | Hom E' ), paste( h[1]^2, ' Global Adjust' ), paste( h[2]^2, ' Global Adjust' ) ),
	cex=1.6, y.intersp=1.35, lwd=2.5 )
legend( 'bottom', bty='n', lty=c(1,NA), pch=c(NA,16), leg=c( 'SD over sims', 'SEs per sims' ), cex=1.8 )

plot.new()

for( dir in dirs ){

load( paste0( dir, 'Rdata/adj.Rdata' ) )

if( dir == dirs[1] ){
	fs	<- xs
} else {
	fs		<- round( seq( .01, .5, len=nx ), 2 )
}
xlim	<- range(fs)


par( mar=c(0.5,0,4.2,.5) )
for( mode in modes ){

	plot( xlim, c(0,1), type='n', axes=F, xlab='', main='', ylab='', cex.lab=2 )
	box()

	mtext( side=3, mains[mode]			, cex=1.65, line=1.7 )

	#axis(1,cex.axis=1.1,at=xs,lab=xs)
	#x1s	<- 0:5/5
	#if( sigtype == 1 )
	if( mode == modes[1] & dir == dirs[1] )
		axis(2,cex.axis=1.4)


	if( mode %in% c( 'ldak' ) ){
		h2s[dir,mode,,,'d1']	<- h2s[dir,mode,,,'d1']	* adjs[,,'1'   ]
		h2s[dir,mode,,,'d2']	<- h2s[dir,mode,,,'d2']	* adjs[,,'2'   ]
		h2s[dir,mode,,,'g']		<- h2s[dir,mode,,,'g']	* adjs[,,'both']
		h2s[dir,mode,,,'hom']	<- h2s[dir,mode,,,'hom']* adjs[,,'both']
		h2s[dir,mode,,,'het']	<- h2s[dir,mode,,,'het']* adjs[,,'both']

		h2.ses[dir,mode,,,'d1']	<- h2.ses[dir,mode,,,'d1']	* adjs[,,'1'   ]
		h2.ses[dir,mode,,,'d2']	<- h2.ses[dir,mode,,,'d2']	* adjs[,,'2'   ]
		h2.ses[dir,mode,,,'g']		<- h2.ses[dir,mode,,,'g']	* adjs[,,'both']
		h2.ses[dir,mode,,,'hom']	<- h2.ses[dir,mode,,,'hom']* adjs[,,'both']
		h2.ses[dir,mode,,,'het']	<- h2.ses[dir,mode,,,'het']* adjs[,,'both']
	}


	h2.ses0	<- apply( h2s[dir,mode,,,], c(1,3), sd, na.rm=T )
	lines( fs, h2.ses0[,'g'		],col='grey'	, lwd=3 )
	lines( fs, h2.ses0[,'hom'	],col=cols[1]	, lwd=3 )
	lines( fs, h2.ses0[,'het'	],col=cols[2]	, lwd=3 )
	lines( fs, h2.ses0[,'d1'	],col=cols[3]	, lwd=3 )
	lines( fs, h2.ses0[,'d2'	],col=cols[4]	, lwd=3 )

	for( j in 1:nx ){
		pointline( fs[j], h2.ses[dir,mode,j,,'g']		, col='grey'	, cex=3 )
		pointline( fs[j], h2.ses[dir,mode,j,,'hom']	, col=cols[1]	, cex=3 )
		pointline( fs[j], h2.ses[dir,mode,j,,'het']	, col=cols[2]	, cex=3 )
		pointline( fs[j], h2.ses[dir,mode,j,,'d1']	, col=cols[3], cex=3 )
		pointline( fs[j], h2.ses[dir,mode,j,,'d2']	, col=cols[4], cex=3 )
	}

}

par( mar=c(5.5,0,.3,.5) )
for( mode in modes ){
	plot( xlim, 0:1, type='n', xlab='', ylab='', axes=F )
	box()

	if( mode == modes[1] & dir == dirs[1] )
		axis(2,cex.axis=1.4)
	axis(1,cex.axis=1.1,at=fs,lab=fs)
	mtext( side=1, ifelse( dir == dirs[1], expression( mu[1] ), 'Population Prevalence' )	, cex=1.35, line=4.1 )

	p1	<- rowMeans( allprev[,,'1']		, na.rm=T )
	p2	<- rowMeans( allprev[,,'2']		, na.rm=T )
	
	lines1( fs, p1, col=cols[3], lty=1, pch=16, cex=1 )
	lines1( fs, p2, col=cols[4], lty=1, pch=16, cex=1 )
	abline( h=.2, lwd=4, col=1 )
	abline( h=0 , lwd=2, col=2, lty=2 )
	abline( h=.4, lwd=2, col=2, lty=2 )

}

rm( allprev, adjs )

}
dev.off()
