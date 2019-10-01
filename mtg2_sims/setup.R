rm( list=ls() )
K0		<- 1
ncaus	<- 1e3

nx		<- 6
xs		<- round( seq( 0, 2.5, len=nx ), 2 )

xlim	<- range(xs)

#maxit	<- 400
maxit	<- 200
cols	<- c( 1, 'orange', '#FB6542', '#3F681C' )

lens	<- c( 333334, 358700, 316584, 325672, 281756, 291988, 257399, 241441, 185975, 225748, 218513, 207067, 160985, 143190, 121926, 125923, 105171, 122970, 84436, 90558, 61487, 52978)
# from `bash count_snps_per_chr.sh > count_snps_per_chr.out' in /ye/zaitlenlabstore/andy/converge_data_nov_2017

modes	<- c( 'ldak', 'pcgc' )

ltys			<- c(  1			, 1											, 1											)
cols1			<- c(  'grey'	, 6											, 'orange'							)
testnames	<- c( 'hom'		, 'iid'									, 'diag3'								)
testlabs	<- c( 'Hom'		, 'IID vs Hom | Hom E', 'IID vs Hom | IID E')
names( ltys )	<- testnames
names( cols1 )<- testnames

mains	<- c( 'REML', 'PCGC' )
names(mains)	<- modes

save.image( 'Rdata/setup.Rdata' )
