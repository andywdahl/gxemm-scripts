rm( list=ls() )
K0		<- 2
ncaus	<- 1e2

nx		<- 7
fs		<- round( seq( .01, .5, len=nx ), 2 )
xlim	<- range(fs)

maxit	<- 200
cols	<- c( 1, 'orange', '#FB6542', '#3F681C' )


ltys			<- c( 1						, 1												, 3												, 1				 , 2												) # 1				,2
cols1			<- c( 'orange'		, 4												, 4												, 5				 , 5												) # 'grey'	,3
testnames	<- c( 'iid'				, 'diag3'									, 'diag1'									, 'h2eq'	 , 'h2eq*'									) #''hom'		,diag2'
testlabs	<- c( 'IID vs Hom', 'Free vs Hom G | Free E', 'Free vs Hom G | Hom E'	,'Equal h2','Equal h2, Global Adjust'	) #'Hom'		,'Free vs Free - Free Noise'
names( ltys )	<- testnames
names( cols1 )<- testnames

modes	<- c( 'ldak', 'pcgc' )
mains	<- c( 'REML', 'PCGC' )
names(mains)	<- modes

N_univ	<- 5e5
N				<- 1e4
L				<- 1e4
Q				<- 10

save.image( 'Rdata/setup.Rdata' )
