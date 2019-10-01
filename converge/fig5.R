rm( list=ls() )
load( file='Rdata/setup.Rdata' )
library(GxEMM)
source( '../code/gxemm.violin.R' )

h2s		<- array( NA, dim=c(length(envs), 14), dimnames=list( envs, c( 'g', 'iid', 'd1', 'd2', 'dum', 'hom', 'het', 'dum1', 'sig20', 'va', 'vb', 'wa-wb', 'wa', 'wb' ) ) )
allps	<- array( NA, dim=c(length(envs), 8	), dimnames=list( envs, c('hom','iid','d1','d2','het','va','vb','wa-wb') ) )
for( phen in ppsub )
try({

	savefile1	<- paste0( 'Rdata/hom/'	, phen, '.Rdata' )
	savefile2	<- paste0( 'Rdata/het/'	, phen, '.Rdata' )
	savefile3	<- paste0( 'Rdata/diag/', phen, '.Rdata' )

	load( savefile1 )
	if( out_hom == 'FAILED' ) file.remove( savefile1 )
	h2s		[phen,'g']	<- out_hom$h2
	allps	[phen,'hom']<- Waldtest( out_hom$h2, out_hom$h2Covmat[1,1] )

	load( savefile2 )
	if( out_het == 'FAILED' ) file.remove( savefile2 )
	h2s		[phen,c('hom','het')]	<- out_het$h2   
	h2s		[phen,'iid']	<- sum(out_het$h2)
	allps [phen,'het']	<- Waldtest( out_het$h2[2]	, out_het$h2Covmat[2,2] )
	allps [phen,'iid']	<- Waldtest( sum(out_het$h2), sum(out_het$h2Covmat) )

	load( savefile3 )
	if( out_diag == 'FAILED' ) file.remove( savefile3 )
	h2s		[phen,c('sig20','va','vb')]	<- out_diag$sig2s[1:3]
	h2s		[phen,c('d1','d2')]		<- out_diag$h2[1:2]
	allps	[phen,c('d1','d2')]	<- sapply( 1:2, function(ii) Waldtest( out_diag$h2[ii], out_diag$h2Covmat[ii,ii] ) )
	allps	[phen,c('va','vb')]	<- sapply( 1:2, function(ii) Waldtest( out_diag$sig2s[1+ii], out_diag$sig2Var[1+ii,1+ii] ) )

	Tmat	<- c(1,-1)
	h2s		[phen,c('wa-wb')]		<- Tmat %*% out_diag$sig2s[4:5]
	h2s		[phen,c('wa','wb')]	<-          out_diag$sig2s[4:5]
	allps	[phen,c('wa-wb')]		<- 2*Waldtest( abs(Tmat%*%out_diag$sig2s[4:5]), Tmat %*% (out_diag$sig2Var[4:5,4:5] %*% Tmat) )

},silent=T)
suff	<- c( 'Fig5', 'Fig14' )
for( ii in 1:2 ){
	if( ii == 1 ){
		ppsub1	<- intersect( ppsub, 1:B )
		ypos1		<- c(-2.32,-2.77,-3.17,2.3)
		limits1	<- c(-1.2, 1.8 )
		ypos2		<- c(-5.0,-6.0,-7.0,5.2)
		limits2	<- c(-2.5, 4.1)
		xlabs1	<- expression( h[g]^2, h[iid]^2, h[stress]^2, h[not]^2, '', h[hom]^2, h[het]^2 )
		xlabs2	<- expression( sigma[g]^2 , v[stress], v[not], w[stress]-w[not] )
		Bp			<- length(ppsub1)
		goodphens	<- ppsub1
	} else if( ii == 2 ){
		ppsub1	<- intersect( ppsub, B+1:P )
		ppsub1	<- ppsub1[c(1:7,10:12)] ### removes height, stress aggregates, pre.menstural.symptoms
		ypos1		<- c(-1.2,-1.5,-1.8,1.7)
		limits1	<- c(-0.5, 1.4 )
		ypos2		<- c(-5.1,-6.2,-7.3,6.5)
		limits2	<- c(-2.2, 5.3)
		xlabs1	<- expression( h[GREML]^2, h[iid]^2, h[hi]^2, h[lo]^2, '', h[hom]^2, h[het]^2 )
		xlabs2	<- expression( sigma[g]^2 , v[hi], v[lo], w[hi]-w[lo] )
		Qp			<- length(ppsub1)
		goodphens	<- union( ppsub1, goodphens )
	}

	pdf( paste0( '~/figs/gxemm/', suff[ii], 'a.pdf' ), width=6.7, height=3.6 )
	gxemm.h2violin(
		h2s[ppsub1,c( 'g', 'iid', 'd1', 'd2', 'dum', 'hom', 'het' )],
		allps[ppsub1,c('hom','iid','d1','d2','het')],
		ypos=ypos1, limits=limits1,
		xlabs	= xlabs1
	)
	dev.off()

	pdf( paste0( '~/figs/gxemm/', suff[ii], 'b.pdf' ), width=3.7, height=3.6 )
	gxemm.sig2violin(
		h2s[ppsub1,c('sig20', 'va', 'vb', 'wa-wb' )],
		allps[ppsub1,c('va', 'vb', 'wa-wb' )],
		ypos=ypos2, limits=limits2,
		xlabs	= xlabs2
	)
	dev.off()
}

out_tab	<- cbind(
	c(Ybnames,Ynames),
	c( rep( 'Binary', B ), rep( 'Quantitative', P ) ),
	h2s[,c(1,2,3,4,6,7,9,10,11,13,14)],
	allps
)[goodphens,]

colnames( out_tab )	<- c(
	'Environment names', 'Environment type',
	'h2g', 'h2iid', 'h2hi', 'h2lo', 'h2hom', 'h2het', 'sig2hom', 'sig2g_hi', 'sig2g_lo', 'sig2e_hi', 'sige2_lo',
	'h2g_p', 'h2iid_p', 'h2hi_p', 'h2lo_p', 'h2het_p', 'sig2g_hi_p', 'sig2g_lo_p', 'delta_sig2e_p'
)
row.names(out_tab)	<- NULL

out_tab	

write.csv( out_tab, file=paste0( '~/figs/gxemm/STab2.csv' ), row.names=FALSE, quote=FALSE )
