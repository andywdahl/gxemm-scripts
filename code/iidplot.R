rm( list=ls() )
library(GxEMM)
source( '../code/plot_fxns_new.R' )
load( file='Rdata/setup.Rdata' )
#load( file='Rdata/data.Rdata' )
#load( 'data/gwas_phens.Rdata' )

#library(ggplot2)
#library(scales)
#
#load( file='Rdata/mysub.Rdata' )
#gwas_phens	<- gwas_phens[mysub]
#gwas_phens	<- setdiff( gwas_phens, 'Area_bc' )
#P	<- length( gwas_phens )

P	<- P+B

phen_names	<- c(Ybnames,Ynames)
type	<- 'base'
nh2		<- 9
h2s		<- array( NA, dim=c(P,P, nh2	), dimnames=list( phen_names, phen_names, c( 'g', 'hom', 'het', 'iid', 'd1', 'd2', 'sig20', 'sig2a', 'sig2b' ) ) )
allps	<- array( NA, dim=c(P,P, 4		), dimnames=list( phen_names, phen_names, c('hom','iid','d1','d2') ) )
for( pp in sample(ppsub) )
	for( qq in sample(qqsub) )
try({

	savefile1	<- paste0( 'Rdata/hom/'	, qq, '_', pp, '_', type, '_allY.Rdata' )
	savefile2	<- paste0( 'Rdata/het/'	, qq, '_', pp, '_', type, '_allY.Rdata' )
	savefile3	<- paste0( 'Rdata/diag/', qq, '_', pp, '_', type, '_allY.Rdata' )

	load( savefile1 )
	h2s		[pp,qq,'g']	<- out_hom$h2
	allps	[pp,qq,'hom']<- Waldtest( out_hom$h2, out_hom$h2Covmat[1,1] )

	load( savefile2 )
	h2s		[pp,qq,c('hom','het')]	<- out_het$h2   
	h2s		[pp,qq,'iid']	<- sum(out_het$h2)
	allps [pp,qq,'iid']	<- Waldtest( out_het$h2[2], out_het$h2Covmat[2,2] )

	load( savefile3 )
	#h2s		[qq,c('hom1')]		<- out_diag$sig2s[1]
	h2s		[pp,qq,c('sig20','sig2a','sig2b')]	<- out_diag$sig2s[1:3]
	h2s		[pp,qq,c('d1','d2')]	<- out_diag$h2[1:2]
	#h2s		[qq,c('d1','d2')]	<- out_diag$sig2s[1+1:2]
	allps		[pp,qq,c('d1','d2')]	<- sapply( 1:2, function(ii) Waldtest( out_diag$sig2s[1+ii], out_diag$sig2Var[1+ii,1+ii] ) )

})

plot.gxemm	<- function(
	h2s,
	allps,
	labs,
	#nlab=length(labs),
	lims=c(-.2,1.2),
	xlabs=expression( h[greml]^2	, h[greml]^2, h[greml]^2	, h[1]^2 ),
	ylabs=expression( h[iid]^2		, h[1]^2	, h[1]^2	, h[2]^2 ),
	cols=c( 1, 1, 'grey', 'orange', 2, 3, 1, 2, 3 )
){
	nh2	<- 9
	is	<- 1:nh2
	x		<- as.factor( matrix( rep( is 	, each=P ), nh2, P, byrow=T )[,qqsub] )
	cols<- c(         matrix( rep( cols	, each=P ), nh2, P, byrow=T )[,qqsub] )
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

		plot( x, y, xlim=lims, ylim=lims, type='n', xlab=xlabs[ii], ylab=ylabs[ii], main='' )
		abline( a=0, b=1, lty=2, col=2 )

		if( length( sub	<- which( pvec < .01 ) ) > 0 ){
			#sub	<- sort.list( abs(y-x), dec=T )[1:nlab]
			points( x[-sub], y[-sub]		, pch=16 )
			text(		x[ sub], y[ sub]		, lab=labs[sub] )
			if( sum(pvec!=0,na.rm=T)!=0 )
			text(		x[ sub], y[ sub]-.10, lab=format( pvec, scientific=T, digits=2 )[sub] )
		} else {
			points( x, y, pch=16 )
		}


	}
}

envs[qqsub]
stop()
for( phen in c( 'Telomere Length', 'Mitochondrial DNA' ) ){
pdf( paste0( '~/figs/gxemm/converge/iid_', phen, '.pdf' ), width=12, height=12 )
par( mar=c(5,5,1,1), mfrow=c(2,2) )
pp	<-  B+which(Ynames==phen)
cat( pp, phen, paste0( '~/figs/gxemm/converge/iid_', phen, '.pdf' ), '\n' )
plot.gxemm(
	h2s=h2s			[pp,qqsub,],
	allps=allps	[pp,qqsub,],
	labs=phen_names[qqsub],
	lims=c(-.5,1.0),
	xlabs=expression( h[greml]^2	, h[greml]^2, h[greml]^2	, h[1]^2 ),
	ylabs=expression( h[iid]^2		, h[1]^2	, h[1]^2	, h[2]^2 )
)
}
dev.off()
