rm( list=ls() )

load( 'data/gwas_phens.Rdata' )
load( 'data/Y_QC.Rdata' ) ### Y has some garbage pruned out relative to Y_raw

P	<- length( gwas_phens )

obs		<- colSums( !is.na( Y )[,gwas_phens]  )

pdf( '~/figs/gxemm/rats/missing.pdf', width=24, height=8 )
par(mar=c(13,4,0,0))

plot( 1:P, obs, type='n', axes=F, ylab='# Samples Observed', xlab='' )
axis(2)
points( 1:P, obs, pch=16 )
abline( h=500, lwd=2, col=2, lty=2 )
axis( 1:P, obs, lab=gwas_phens, at=1:P, las=2 )

dev.off()
mysub	<- which( obs > 1000 )
length( mysub )

nrow(Y)
mean( obs )
mean( obs[mysub] )

save( mysub, file='Rdata/mysub.Rdata' )
