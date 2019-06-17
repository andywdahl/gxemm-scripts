rm( list=ls() )
K0		<- 2
ncaus	<- 1e3

ws	<- c( .5, .5 )

nx		<- 7
xs		<- seq( 0, .6, len=7 )
xlim	<- c(-.02,.62)

M.all	<- 666

maxit	<- 200
cols	<- c( 1, 'orange', '#FB6542', '#3F681C', '#375E97' )

ltys			<- c(  1						, 3										,1						, 1								, 3												, 2							,1					) #, 1
cols1			<- c( 'grey'				,'grey'								,'orange'			, 4 							, 4												, 6							,5					) #, 4
testnames	<- c( 'hom'					,'hom1'								,'iid'				, 'diag4' 				, 'diag3' 								, 'diag1'				,'h2eq'			) #, 'diag2'
testlabs	<- c( 'Hom vs Null' ,'Hom vs Null | Free'	,'IID vs Hom'	, 'Free vs Hom G'	, 'Free vs Hom G | Hom E'	, 'Free vs IID'	, 'Equal h2') #, 'Free vs Hom'
names( ltys )	<- testnames
names( cols1 )<- testnames

all_sig2homs	<- function( xval )
	c( xs[xval], 0, 0, 0, .1 )
all_sig2hets	<- function( xval )
	list(
		rep( 0, K0 ),
		rep( xs[xval], K0 ),
		10/6*c( .6-xs[xval], xs[xval] ),
		c( .3, xs[xval] )/sum(ws*c( .3, xs[xval] ))*.1,
		rep( 0, K0 )
	)
all_tauhets	<- function( xval )
	list( NA, NA, NA,
		c( .3, xs[xval] )/sum(ws*c( .3, xs[xval] ))*.9,
		2.75*c( xs[xval], .6-xs[xval] )
	)

xlabs	<- expression( h[hom]^2, h[het]^2, v[2], omega, w[2] )
mains	<- c( 'Homogeneity', 'IID GxE, Hom Noise', 'Free GxE, Hom Noise', 'Free with equal h2', 'Hom G, Free Noise' )

lens	<- c( 333334, 358700, 316584, 325672, 281756, 291988, 257399, 241441, 185975, 225748, 218513, 207067, 160985, 143190, 121926, 125923, 105171, 122970, 84436, 90558, 61487, 52978)
# from `bash count_snps_per_chr.sh > count_snps_per_chr.out' in /ye/zaitlenlabstore/andy/converge_data_nov_2017

save.image( 'Rdata/setup.Rdata' )
