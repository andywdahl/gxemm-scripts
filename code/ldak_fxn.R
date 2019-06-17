Rldak <- function( A, y, X, verbose=TRUE, trunc=T, M.all,
	ldak_loc='~/GxEMM/code/ldak5.linux '
	IDs=1:length(y),
	tmpdir=paste0( 'gcta_tmp_', round( abs(rnorm(1) + X[sample(10,1)] + y[sample(10,1)] + length(A) )*1e5, 0 ) ),
	outdir=tmpdir,
	#lrt_inds=NULL,
	kill_tmpdir=TRUE
){
	r <- length(A)

	#### check not already ran
	if( file.exists( paste0( tmpdir, '/y.phen' ) ) ) stop()
	system( paste0( 'mkdir ', tmpdir ) )
	if( tmpdir != outdir )
	system( paste0( 'mkdir ', outdir ) )

	### write pheno + covars
	write_gcta_pheno( y, file=paste0( tmpdir, '/y.phen')	, IDs=IDs )
	write_gcta_pheno( X, file=paste0( tmpdir, '/X.qcovar'), IDs=IDs )

	### write each of the A to a kinship file
	for( i in 1:r )
		write_gcta_kinship( A[[i]], prefix=paste0( tmpdir, '/K.', i ), IDs=IDs, M.all=M.all )
	write.table( paste0( tmpdir, '/K.', 1:r ), file=paste0( tmpdir, '/multi_grm.txt' ), sep='\t', row.names=F, col.names=F, quote=F )

	failed	<- TRUE
	try({
		### pass kinship files and phenotypes/covariates to GCTA
		system( paste0(
			gcta_loc,
			ifelse( r==1, 
				paste0(	'--grm-gz ', tmpdir, '/K.1 '  ),
				paste0( '--mgrm-gz ', tmpdir, '/multi_grm.txt ' )
			), 
			' --pheno ', tmpdir, '/y.phen --reml',
			' --qcovar ', tmpdir, '/X.qcovar ',
			ifelse( trunc, '', ' --reml-no-constrain ' ),
			#ifelse( is.null(lrt_inds), '', paste0( ' --reml-lrt ', lrt_inds, collapse=' ' ) ),
			' --out ', outdir, '/tmp'
		))

		### read GCTA
		out		<- read.table( paste0( outdir, '/tmp.hsq' ), fill=T, head=T, row.names=1 ) 
		sig2s	<- as.numeric( as.character( out[1:(r+2),'Variance'] ) )
		sig2ses	<- as.numeric( as.character( out[1:(r+2),'SE'] ) )
		h2s		<- as.numeric( as.character( out[r+2+1:r,'Variance'] ) )
		ses		<- as.numeric( as.character( out[r+2+1:r,'SE'] ) )
		pval	<- as.numeric( as.character( out['Pval',1] ) )
		ll		<- as.numeric( as.character( out['logL'	,1] ) )
		ll0		<- as.numeric( as.character( out['logL0',1] ) )

		nrows	<- length( scan( paste0( outdir, '/tmp.log' ), what='character', sep='\n', blank.lines.skip=F ) )
		vcov	<- read.table( paste0( outdir, '/tmp.log' ), fill=T, skip=nrows-8-(r+1) )
		toprow	<- which( as.character( vcov[,1] ) == 'Sampling' )

		sig2Var	<- matrix( as.numeric( as.character( unlist( vcov[toprow+1:(r+1),1:(r+1)] ) ) ), r+1, r+1 )

		failed	<- FALSE
	})
	if( kill_tmpdir )
	system( paste0( 'rm -rf ', tmpdir ) )
	if( failed ) return('FAILED')
	return( list( h2s=h2s, sig2s=sig2s, sig2ses=sig2ses, ses=ses, ll=ll, ll0=ll0, pval=pval, sig2Var=sig2Var ) )
}

write_gcta_pheno	<- function(y,file,IDs)
	write.table( cbind( IDs, IDs, y ), file=file, sep='\t', row.names=F, col.names=F, quote=F )

write_gcta_kinship	<- function(K,prefix,IDs,M.all){
	keep = which( lower.tri(K,diag=T) , arr.ind=T )
	non_mis = rep( M.all  , length(keep[,1]) )
	keep = cbind( keep , non_mis , K[ keep ] )
	keep = keep[ order(keep[,1],keep[,2]) , ]
	write.table( keep , file=gzfile(paste0( prefix, '.grm.gz' )) , col.names=F , row.names=F , quote=F )
	write.table( cbind( IDs, IDs ), file=paste0( prefix, '.grm.id' ), col.names=F , row.names=F , quote=F )
}
