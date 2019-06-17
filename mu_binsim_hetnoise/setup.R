rm( list=ls() )
K0		<- 2
ncaus	<- 1e3

nx		<- 7
xs		<- round( seq( 0, .6, len=nx ), 2 )

xlim	<- range(xs)

maxit	<- 200
cols	<- c( 1, 'orange', '#FB6542', '#3F681C' )

lens	<- c( 333334, 358700, 316584, 325672, 281756, 291988, 257399, 241441, 185975, 225748, 218513, 207067, 160985, 143190, 121926, 125923, 105171, 122970, 84436, 90558, 61487, 52978)
# from `bash count_snps_per_chr.sh > count_snps_per_chr.out' in /ye/zaitlenlabstore/andy/converge_data_nov_2017

modes	<- c( 'ldak', 'pcgc' )

ltys			<- c(  1					, 1												, 3												, 1				 , 2												) # 1				,2
cols1			<- c(  'orange'		, 4												, 4												, 5				 , 5													) # 'grey'	,3
testnames	<- c(  'iid'			, 'diag3'									, 'diag1'									, 'h2eq'	 , 'h2eq*'									) #''hom'		,diag2'
testlabs	<- c( 'IID vs Hom', 'Free vs Hom G | Free E', 'Free vs Hom G | Hom E'	,'Equal h2','Equal h2, Global Adjust'	) #'Hom'		,'Free vs Free - Free Noise'
names( ltys )	<- testnames
names( cols1 )<- testnames

mains	<- c( 'REML', 'PCGC' )
names(mains)	<- modes

save.image( 'Rdata/setup.Rdata' )
