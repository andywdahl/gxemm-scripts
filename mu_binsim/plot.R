rm( list=ls() )
library(GxEMM)
load( 'Rdata/setup.Rdata' )
source( '../code/plot_fxns.R' )

dirs	<- c( '../mu_binsim/', '../asc_binsim/' )

allps		<- array( NA, dim=c( 2, 2, nx, maxit, 5), dimnames=list( dirs, modes, 1:nx, 1:maxit, testnames ) )
h2s			<- array( NA, dim=c( 2, 2, nx, maxit, 7), dimnames=list( dirs, modes, 1:nx, 1:maxit, c('g','hom','het','d1','d2','d1b','d2b') ) )
for( it in 1:maxit )
	for( mode in modes )
		for( xval in sample(nx) )
			for( dir in dirs )
try({

	load( paste0( dir, 'Rdata/adj.Rdata' ) )
	if( dir == dirs[1] ){
		fs	<- xs
	} else {
		fs		<- round( seq( .01, .5, len=nx ), 2 )
	}

	savefile1	<- paste0( dir, 'Rdata/hom/'	, fs[xval], '_', it, '_', mode, '.Rdata' )
	savefile2	<- paste0( dir, 'Rdata/het/'	, fs[xval], '_', it, '_', mode, '.Rdata' )
	savefile3	<- paste0( dir, 'Rdata/diag/'	, fs[xval], '_', it, '_', mode, '.Rdata' )
	savefile4	<- paste0( dir, 'Rdata/diag1/', fs[xval], '_', it, '_', mode, '.Rdata' )

	load( savefile1 )
	if( out_hom == 'FAILED' ) file.remove( savefile1 )
	h2s		[dir,mode,xval,it,'g']		<- out_hom$h2
	#allps	[mode,xval,it,'hom']	<- Waldtest( out_hom$h2, out_hom$h2Covmat[1,1] )

	load( savefile2 )
	if( out_het == 'FAILED' ) file.remove( savefile2 )
	h2s		[dir,mode,xval,it,c('hom','het')]	<- out_het$h2[c('hom','het')]
	allps	[dir,mode,xval,it,'iid']	<- Waldtest( out_het$h2[2], out_het$h2Covmat[2,2] )

	load( savefile3 )
	if( out_diag == 'FAILED' ) file.remove( savefile3 )
	h2s		[dir,mode,xval,it,c('d1b','d2b')]	<- out_diag$h2
	allps	[dir,mode,xval,it,'diag1']	<- MVWaldtest( out_diag$sig2s[1+1:2], out_diag$sig2Var[1+1:2,1+1:2] )

	load( savefile4 )
	if( out_diag1 == 'FAILED' ) file.remove( savefile4 )
	h2s		[dir,mode,xval,it,c('d1','d2')]	<- out_diag1$h2
	#allps	[mode,xval,it,'diag2']	<- MVWaldtest( out_diag1$sig2s[1+1:4], out_diag1$sig2Var[1+1:4,1+1:4] )
	allps	[dir,mode,xval,it,'diag3']	<- MVWaldtest( out_diag1$sig2s[1+1:2], out_diag1$sig2Var[1+1:2,1+1:2] )
	if( mode == 'ldak' ){
		allps	[dir,mode,xval,it,'h2eq*']<- h2_equal_test( out_diag1$h2			,									out_diag1$h2Covmat )
		aadj	<- adjs[xval,it,c('1','2')]
		allps	[dir,mode,xval,it,'h2eq']	<- h2_equal_test( aadj*out_diag1$h2	, diag(aadj) %*%	out_diag1$h2Covmat %*% diag(aadj) )
		rm( aadj )
	} else {
		allps	[dir,mode,xval,it,'h2eq']	<- h2_equal_test( out_diag1$h2      , out_diag1$h2Covmat )
	}

	rm( out_hom, out_het, out_diag, out_diag1 )
	rm( fs, adjs )
},silent=T)

print( apply( ! is.na(h2s), c(1,2,3,5), sum )[1,,,c('g','hom','d1','d1b')] )
print( apply( ! is.na(h2s), c(1,2,3,5), sum )[2,,,c('g','hom','d1','d1b')] )

pdf( paste0( '~/figs/gxemm/Fig3.pdf' ), width=16.5, height=6.9 )

layout( cbind( 1:3, 6+matrix( c(3:4,1:2,5:6), 3, 2, byrow=T ), rep(6,3), 12+matrix( c(3:4,1:2,5:6), 3, 2, byrow=T ), c(4, 5, 5) ), widths=c( 1.5, rep(4.3,2), .2, rep(4.3,2), 4.0 ), heights=c(4.5,6.99,3.5) )

par( mar=c(.2,0,5.2,.5) )
plot.new()
mtext( side=2, 'False Positive Rate\n(nominal p<.05)', cex=1.2, line=-4.0 )

par( mar=c(.9,0,.9,.5) )
plot.new()
mtext( side=2, 'Heritability\nEstimates'	, cex=1.2, line=-4.0 )

par( mar=c(5.5,0,.4,.5) )
plot.new()
mtext( side=2, 'Disease\nPrevalence'	, cex=1.2, line=-4.0 )

par( mar=c(0,0,0,0) )
plot.new()
legend( 'bottom', bty='n',
	lwd=c( 3, 3, 5, 3, 3, 5 ),
	col=c( cols1, 'lightgrey' ),
	lty=c( ltys, 3 ),
	leg=c( testlabs, 'Null' ), cex=1.7, y.intersp=1.1 )

par( mar=c(5,0,1,0) )
plot.new()
legend( 'top', bty='n',
	lty=c( 1, 1 ,1, 1, 3, 3, 2, 2 ),
	pch=c( 16,16,16,16,5, 5, 1, 1 ),
	col=c( cols[c(1,2,3,4,3:4)], 6, 3 ),
	leg=expression(h[g]^2, h[het]^2, h[1]^2, h[2]^2, paste( h[1]^2, ' | Hom E' ), paste( h[2]^2, ' | Hom E' ), paste( h[1]^2, ' Global Adjust' ), paste( h[2]^2, ' Global Adjust' ) ),
	cex=1.6, y.intersp=1.35, lwd=2.5 )
legend( 'bottom', bty='n',
	lty=2,
	col=cols[3:4],
	leg=c('Prevalence in Env. 1','Prevalence in Env. 1'),
	cex=1.6, y.intersp=1.35, lwd=2.5 )

plot.new()


for( dir in dirs ){

load( paste0( dir, 'Rdata/adj.Rdata' ) )
if( dir == dirs[1] ){
	fs	<- xs
} else {
	fs		<- round( seq( .01, .5, len=nx ), 2 )
}
xlim	<- range(fs)


par( mar=c(.3,0,.3,.5) )
for( mode in modes ){

	plot( xlim, c(-.05,1.05), type='n', axes=F, xlab='', main='', ylab='', cex.lab=2 )
	box()
	#text(.01,.684, letters[2+which(modes==mode)],cex=2.2)

	x1s	<- 0:5/5
	if( mode == modes[1] & dir == dirs[1] )
		axis(2,cex.axis=1.1,at=x1s,lab=paste0(x1s*100,'%'))

	lines( fs, fs*0+.35,col='grey', lwd=3, lty=3 )
	lines( fs, fs*0    ,col='grey', lwd=3, lty=3 )

	if( mode %in% c( 'ldak' ) ){
	lines1( fs, rowMeans( h2s[dir,mode,,,'d1']	* adjs[,,'both'], na.rm=T ), col=6      , lty=2, lwd=3, pch=1 )
	lines1( fs, rowMeans( h2s[dir,mode,,,'d2']	* adjs[,,'both'], na.rm=T ), col=3      , lty=2, lwd=3, pch=1 )
	lines1( fs, rowMeans( h2s[dir,mode,,,'d1']	* adjs[,,'1'   ], na.rm=T ), col=cols[3], lty=1, lwd=3, pch=16)
	lines1( fs, rowMeans( h2s[dir,mode,,,'d2']	* adjs[,,'2'   ], na.rm=T ), col=cols[4], lty=1, lwd=3, pch=16)
	lines1( fs, rowMeans( h2s[dir,mode,,,'g']		* adjs[,,'both'], na.rm=T ), col=cols[1], lty=1, lwd=3, pch=16)
	lines1( fs, rowMeans( h2s[dir,mode,,,'het']	* adjs[,,'both'], na.rm=T ), col=cols[2], lty=1, lwd=3, pch=16)

	lines1( fs, rowMeans( h2s[dir,mode,,,'d1b']	* adjs[,,'1'   ], na.rm=T ), col=cols[3], lty=3, lwd=3, pch=5, cex=3 )
	lines1( fs, rowMeans( h2s[dir,mode,,,'d2b']	* adjs[,,'2'   ], na.rm=T ), col=cols[4], lty=3, lwd=3, pch=5, cex=3 )

	} else {
	lines1( fs, rowMeans( h2s[dir,mode,,,'d1']	                , na.rm=T ), col=cols[3], lty=1, lwd=3, pch=16)
	lines1( fs, rowMeans( h2s[dir,mode,,,'d2']	                , na.rm=T ), col=cols[4], lty=1, lwd=3, pch=16)
	lines1( fs, rowMeans( h2s[dir,mode,,,'d1b']	                , na.rm=T ), col=cols[3], lty=3, lwd=3, pch=5, cex=3 )
	lines1( fs, rowMeans( h2s[dir,mode,,,'d2b']	                , na.rm=T ), col=cols[4], lty=3, lwd=3, pch=5, cex=3 )
	lines1( fs, rowMeans( h2s[dir,mode,,,'g']										, na.rm=T ), col=cols[1], lty=1, lwd=3, pch=16 )
	lines1( fs, rowMeans( h2s[dir,mode,,,'het']									, na.rm=T ), col=cols[2], lty=1, lwd=3, pch=16 )
	}
}

par( mar=c(.3,0,3.3,.5) )
for( mode in modes ){
	plot( xlim, c(0,1), type='n', xlab='', ylab='', axes=F )
	box()
	mtext( side=3, mains[mode]			, cex=1.65, line=1.0 )
	x1s	<- 0:2/2
	if( mode == modes[1] & dir == dirs[1] )
		axis(2,cex.axis=1.1,at=x1s,lab=paste0(x1s*100,'%'))

	lines(	fs, .05+0*fs, col='lightgrey', lty=3, lwd=3	)
	for( ii in c( 'iid', 'diag3', 'h2eq', 'h2eq*', 'diag1' ) )
		lines1( fs, apply(allps[dir,mode,,,ii]< .05,1,mean,na.rm=T), lty=ltys[ii], col=cols1[ii] )#, y.max=10
}

par( mar=c(4.5,0,.3,.5) )
for( ii in modes ){
	plot( xlim, 0:1, type='n', xlab='', ylab='', axes=F )
	box()

	mtext( side=1, ifelse( dir == dirs[1], expression( mu[1] ), 'Population Prevalence' ), cex=1.35, line=3.2 )
	axis(1,cex.axis=1.1,at=fs,lab=fs)
	if( ii == modes[1] & dir == dirs[1] )
	axis(2,cex.axis=1.1,at=0:4/10,lab=paste0(0:4*10,'%'))

	lines1( fs, rowMeans( allprev[,,'1'], na.rm=T ), col=cols[3], lty=2, pch=16, cex=1 )
	lines1( fs, rowMeans( allprev[,,'2'], na.rm=T ), col=cols[4], lty=2, pch=16, cex=1 )
	abline( h=.2, lwd=4, col=1 )
	abline( h=0 , lwd=2, col=2, lty=2 )
	abline( h=.4, lwd=2, col=2, lty=2 )

}

rm( allprev, adjs )

}
dev.off()
source( 'seplot.R' )
