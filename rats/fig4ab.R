rm( list=ls() )
library(GxEMM)
source( '../code/gxemm.violin.R' )


type	<- 'base'
load( 'data/gwas_phens.Rdata' )
load( file='Rdata/mysub.Rdata' )
gwas_phens	<- gwas_phens[mysub]
gwas_phens	<- setdiff( gwas_phens, 'Area_bc' )
P	<- length( gwas_phens )


nh2		<- 14
h2s		<- array( NA, dim=c(P, nh2), dimnames=list( gwas_phens, c( 'g', 'iid', 'd1', 'd2', 'dum', 'hom', 'het', 'dum1', 'sig20', 'va', 'vb', 'wa-wb', 'wa', 'wb' ) ) )
allps	<- array( NA, dim=c(P, 10	), dimnames=list( gwas_phens, c('hom','iid','d1','d2','het','va','vb','diffg','wa-wb', 'free-iid') ) )
for( phen in gwas_phens )
try({

	savefile1	<- paste0( 'Rdata/hom/'	, phen, '_', type, '.Rdata' )
	savefile2	<- paste0( 'Rdata/iid/'	, phen, '_', type, '.Rdata' )
	savefile3	<- paste0( 'Rdata/free/', phen, '_', type, '.Rdata' )

	load( savefile1 )
	h2s		[phen,'g']	<- out_hom$h2
	allps	[phen,'hom']<- Waldtest( out_hom$h2, out_hom$h2Covmat[1,1] )

	load( savefile2 )
	h2s		[phen,c('hom','het')]	<- out_het$h2   
	h2s		[phen,'iid']	<- sum(out_het$h2)
	allps [phen,'het']	<- Waldtest( out_het$h2[2]	, out_het$h2Covmat[2,2] )
	allps [phen,'iid']	<- Waldtest( sum(out_het$h2), sum(out_het$h2Covmat) )

	load( savefile3 )
	h2s		[phen,c('sig20','va','vb')]	<- out_diag$sig2s[1:3]
	h2s		[phen,c('wa-wb')]			<- out_diag$sig2s[4]
	h2s		[phen,c('wa','wb')]		<- c( out_diag$sig2s[4], 0 ) + out_diag$sig2s[5]
	h2s		[phen,c('d1','d2')]		<- out_diag$h2[1:2]

	allps	[phen,c('d1','d2')]	<- sapply( 1:2, function(ii) Waldtest( sum(out_diag$sig2s[c(1,1+ii)]), sum(out_diag$sig2Var[c(1,1+ii),c(1,1+ii)]) ) )
	allps	[phen,c('va','vb')]	<- sapply( 1:2, function(ii) Waldtest( out_diag$sig2s[1+ii], out_diag$sig2Var[1+ii,1+ii] ) )
	allps	[phen,c('wa-wb')]		<- 2*min(
		Waldtest( out_diag$sig2s[4], out_diag$sig2Var[4,4] ),
		Waldtest( -out_diag$sig2s[4], out_diag$sig2Var[4,4] )
	)

	Tmat	<- cbind( c(1,-1,0), c(0,0,1) )
	allps	[phen,'free-iid']	<- MVWaldtest( t(Tmat) %*% out_diag$sig2s[2:4], t(Tmat) %*% out_diag$sig2Var[2:4,2:4] %*% Tmat )

	Tmat	<- cbind( 1,-1 )
	allps	[phen,'diffg']	<- Waldtest( abs( Tmat %*% out_diag$sig2s[2:3] ), Tmat %*% out_diag$sig2Var[2:3,2:3] %*% t(Tmat) )

})

pdf( paste0( '~/figs/gxemm/Fig4a.pdf' ), width=6.7, height=3.6 )
gxemm.h2violin(
	h2s[,c('g', 'iid', 'd1', 'd2', 'dum', 'hom', 'het')],
	allps[,c('hom','iid','d1','d2','het')],
	xlabs	= expression( h[g]^2, h[iid]^2, h[female]^2, h[male]^2, '', h[hom]^2, h[het]^2, '', sigma[g]^2 , v[fe], v[ma], w[fe]-w[ma] ),
	ypos=c( -.75,-.89,-1.04,1.38,1.15),
	limits=c(-0.2,1.1)
)
dev.off()

pdf( paste0( '~/figs/gxemm/Fig4b.pdf' ), width=3.7, height=3.6 )
gxemm.sig2violin(
	h2s[,c('sig20', 'va', 'vb', 'wa-wb' )],
	allps[,c('va', 'vb', 'wa-wb' )],
	xlabs	= expression( sigma[g]^2 , v[F], v[M], w[F]-w[M] ),
	limits=c(-0.5,2.8)
)
dev.off()


info		<- read.csv( 'raw_data/supp_tab_1.csv' )
infosub	<- which( info[,'Measure.in.GSCAN'] %in% gwas_phens )
info		<- info[infosub,]

out_tab	<- cbind(
	as.character( info[,'Measure'] ), as.character( info[,'Phenotype'] ),
	row.names(h2s),
	h2s[,c(1,2,3,4,6,7,9,10,11,13,14)],
	allps
)

row.names(out_tab)	<- NULL
colnames( out_tab )	<- c(
	'Trait names', 'Trait Category', 'Short trait names',
	'h2g', 'h2iid', 'h2female', 'h2male', 'h2hom', 'h2het', 'sig2hom', 'sig2g_female', 'sig2g_male', 'sig2e_female', 'sige2_male',
	'h2g_p', 'h2iid_p', 'h2female_p', 'h2male_p', 'h2het_p', 'sig2g_female_p', 'sig2g_male_p', 'delta_sig2g_p', 'delta_sig2e_p', 'free_vs_iid_p'
)

write.csv( out_tab, file=paste0( '~/figs/gxemm/STab1.csv' ), row.names=FALSE, quote=FALSE )

sort(allps[,'het']			)[1:13]
sort(allps[,'va']				)[1:5]
sort(allps[,'vb']				)[1:5]
sort(allps[,'wa-wb']		)[1:6]
sort(allps[,'free-iid'] )[1:10]
sort(allps[,'diffg']		)[1:10]
