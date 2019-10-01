Rmtg2	<- function( y, X, Z, mc.cores=1, suffix='', resX=TRUE, type=c( 'mrnm', 'rnm', 'rnm+', 'hom' ), IDs ){
	if( ncol(Z) != 1 ) stop()
	Z	<- scale(Z[,1])
	#### change below files for your own application
  	famfile   <-        '/wynton/scratch/gxemm/conv/tmp.fam'
	Kfile			<-				'/wynton/scratch/gxemm/conv/K.for_MTG2.grm.gz'
	phenfile	<- paste0('/wynton/scratch/gxemm/conv/tst.dat'	, suffix, '_', type )
	Zfile			<- paste0('/wynton/scratch/gxemm/conv/mtg2.par'	, suffix, '_', type )

	##### write y
	if( resX ){
		X		<- cbind( 1, X )
		PX	<- X %*% solve( t(X) %*% X ) %*% t(X)
		y		<- y - PX %*% y
	}
	y	<- scale(y)
	if( type %in% c( 'mrnm', 'rnm+' ) ){
		stop('I could never get these working')
	} else {
		write.table( cbind(IDs,IDs,y)		, file=phenfile, sep=' ', row.names=F, col.names=F, quote=F )
	}

	##### write Z, and for some reason y again
	N <- length(y)
	if( type == 'rnm' ){
	  write( N  , ncolumns=1, file=Zfile )
		write( 2  , ncolumns=1, file=Zfile, append=T )
		write( Z  , ncolumns=1, file=Zfile, append=T )
	}

	##### wrapper for MTG2
	if( type == 'rnm' ){
	syscall	<- paste0(
		' /netapp/home/andyd/GxEMM/code/mtg2 ',
		' -p ', famfile,
		' -d ', phenfile,
		' -zg ', Kfile,
		' -rnm ', Zfile,
		' -mrnm 1 ', # ``where	-mrnm	specifies	the	number	of	traits''
		' -mod 1 ', # ``--mod	{number	of	traits}''
		' -out Rout/mtg2out/out', suffix, '_', type
	)
	} else if( type == 'hom' ){
	syscall	<- paste0(
		' /netapp/home/andyd/GxEMM/code/mtg2 ',
		' -p ', famfile,
		' -d ', phenfile,
		' -zg ', Kfile,
		' -mod 1 ',
		' -out Rout/mtg2out/out', suffix, '_', type
	)
	}
	print( syscall )
	system( syscall )
	system( paste0( 'rm ', phenfile ) )

	readMTG2( paste0( 'Rout/mtg2out/out', suffix, '_', type ), type )
}

readMTG2	<- function( loadfile, type=c( 'mrnm', 'rnm', 'hom' ) ){
	xx	<- as.matrix( read.table( loadfile, fill=T, header=FALSE ) )
	if( type == 'hom' ){
		stopifnot( all( xx[1:4,1] == c( "Ve" , "Va" , "h2" , "LKH" ) ) )
		Ve	<- as.numeric( xx[1,2:3] )
		Vg	<- as.numeric( xx[2,2:3] )
		Vgxe<- NA
		ll	<- as.numeric( xx[4,2] )
	} else if( type == 'rnm' ){
		stopifnot( all( xx[1:5,1] == c( "Ve" , "Vk" , "Vk" , "cov" , "LKH" ) ) )
		Ve	<- as.numeric( xx[1,2:3] )
		Vg	<- as.numeric( xx[2,2:3] )
		Vgxe<- as.numeric( xx[3,2:3] )
		ll	<- as.numeric( xx[5,2] )
	}
	list( Ve=Ve, Vg=Vg, Vgxe=Vgxe, ll=ll )
}
