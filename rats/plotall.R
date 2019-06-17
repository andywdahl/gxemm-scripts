rm( list=ls() )
library(GxEMM)
source( '../code/gxemm.plot.R' )

type	<- 'base'
load( 'data/gwas_phens.Rdata' )
load( file='Rdata/mysub.Rdata' )

for( phen in gwas_phens )
try({

	savefile1	<- paste0( 'Rdata/hom/'	, phen, '_', type, '.Rdata' )
	savefile2	<- paste0( 'Rdata/iid/'	, phen, '_', type, '.Rdata' )
	savefile3	<- paste0( 'Rdata/free/', phen, '_', type, '.Rdata' )

	load( savefile1 )
	load( savefile2 )
	load( savefile3 )
	gxemm.plot( out_hom, out_het, out_diag, Enames=c('F','M'), pdfname=paste0( '~/figs/gxemm/rats/', phen, '.pdf' ) )
	rm( out_hom, out_het, out_diag )

})
