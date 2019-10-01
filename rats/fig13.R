rm( list=ls() )
library(ggplot2)
library(scales)
library(GxEMM)
source( '../code/plot_fxns_new.R' )
load( 'data/gwas_phens.Rdata' )


load( file='Rdata/mysub.Rdata' )
gwas_phens	<- gwas_phens[mysub]
gwas_phens	<- setdiff( gwas_phens, 'Area_bc' )
P	<- length( gwas_phens )

type	<- 'base'

h2s		<- array( NA, dim=c(P,9), dimnames=list( gwas_phens, c( 'g', 'hom', 'het', 'iid', 'd1', 'd2', 'sig20', 'sig2a', 'sig2b' ) ) )
for( phen in gwas_phens )
try({

	savefile1	<- paste0( 'Rdata/hom/'	, phen, '_', type, '.Rdata' )
	savefile2	<- paste0( 'Rdata/iid/'	, phen, '_', type, '.Rdata' )
	savefile3	<- paste0( 'Rdata/free/', phen, '_', type, '.Rdata' )

	load( savefile1 )
	h2s		[phen,'g']	<- out_hom$h2

	load( savefile2 )
	h2s		[phen,c('hom','het')]	<- out_het$h2   
	h2s		[phen,'iid']	<- sum(out_het$h2)

	load( savefile3 )
	h2s		[phen,c('sig20','sig2a','sig2b')]	<- out_diag$sig2s[1:3]
	h2s		[phen,c('d1','d2')]	<- out_diag$h2

})

pdf( paste0( '~/figs/gxemm/Fig13.pdf' ), width=12, height=12 )
par( mar=c(5,5,1,1), mfrow=c(2,2) )

xlabs	<- expression( h[greml]^2	, h[greml]^2, h[greml]^2	, h[female]^2 ) ### first column of Z is female
ylabs	<- expression( h[iid]^2		, h[female]^2	, h[male]^2	, h[male]^2 )
for( ii in 1:4 ){

	x	<- h2s[,c( 'g'	, 'g'	, 'g'	, 'd1' )[ii]]
	y	<- h2s[,c( 'iid', 'd1', 'd2', 'd2' )[ii]]

	plot( x, y, xlim=c(.2,1), ylim=c(.2,1), type='n', xlab=xlabs[ii], ylab=ylabs[ii], main='' )
	abline( a=0, b=1, lty=2, col=2 )
	sub	<- sort.list( abs(y-x), dec=T )[1:10]
	points( x[-sub], y[-sub], pch=16 )
	text(		x[ sub], y[ sub], pch=16, lab=gwas_phens[sub] )

}
dev.off()
