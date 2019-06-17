rm( list=ls() )
sink( 'Rout/setup.Rout' )

load( paste0( 'parsed_data/', 'base', '.Rdata' )  )
X			<- scale(G)

Yb2Z	<- function(yb){
	Z		<- t(sapply( yb, function(x){
		if( is.na(x) ){
			return( c( NA, NA ) )
		}	else if( x == 0 ){
			return( c( 0, 1 ) )
		}	else if( x == 1 ){
			return( c( 1, 0 ) )
		}
	}))
	Z * sqrt(nrow(Z)/sum(Z^2))
}

## load K
load( 'parsed_data/K.Rdata' )
K	<- K[rownames(X),rownames(X)]
K	<- K*nrow(K)/sum(diag(K))

y		<- 1+as.numeric( Yb[,1] == max(Yb[,1],na.rm=T) )

stressall		<- rowSums( Yb[,paste0( 'LS.', 1:16 )] )
stresspc		<- svd( Yb[,paste0( 'LS.', 1:16 )] )$u[,1]
stresspc2		<- svd( scale(Yb[,paste0( 'LS.', (1:16)[-c(11,15)] )]) )$u[,1]
Y						<- cbind( Y, stressall, stresspc, stresspc2 )

save( Yb2Z, y, X, K, Y, Yb, file='Rdata/data.Rdata' )

ppsub	<- c(7,8,16:ncol(Yb),ncol(Yb)+2, ncol(Yb)+4:(ncol(Y)))
qqsub	<- c(1:6						,ncol(Yb)+c(0,1:3,10,11,12,13,14,15 ))

ppsub	<- ppsub[-c(4,13,17)] # LS.2 is redundant with sepdivwid; LS.11 has only 3 in ctrls; LS.15 has only 17
# LS.6/13 have only 52/43 in ctrls

envs	<- c(colnames(Yb),colnames(Y))
envs[ppsub]
envs[qqsub]


envs[qqsub]
qqsub	<- qqsub[-c(9,15)]
envs[qqsub]



rbind(
round( colSums( Yb[Yb[,1]==1,intersect( ppsub, 1:ncol(Yb) )], na.rm=T ), 2 ),
round( colSums( Yb[Yb[,1]==0,intersect( ppsub, 1:ncol(Yb) )], na.rm=T ), 2 )
)

Ybnames	<- c( 'Major Depression', 'Melancholia', "Panic Disorder", 'Anxiety Disorder', "Dysthymia", "Postnatal Depression",
	"CSA", "Stress", "Father MD", "Mother MD", 'Agoraphobia', 'Social Phobia', 'Animal Phobia', 'Situational Phobia', 'Blood Phobia',
	'Close Family Death', 'Divorced_Separated', 'Unemployed', 'Fired', 'Finanical Crisis', 'Legal Probems', 'Serious Illness',
	'Serious Accident', 'Natural Disaster', 'Witnessed Violence', 'Raped', 'Physically Attacked', 'Abused as Child ', 'Neglected as Child', 'Violently Threatened', 'Other Terrible Events', 'Widowed' )
Ynames	<- c( 'Total MD Episodes', 'Neuroticism', 'Family History', 'Cold Mother', 'Authoritarian Mother', 'Protective Mother', 'Cold Father', 'Authoritarian Father', 'Protective Father','Pre-Menstrual MD', 'Height', 'BMI', 'Mitochondrial DNA', 'Telomere Length', 'StressAll', 'StressPC', 'StressPC2' )

rbind(
Ybnames,
round( colSums( Yb[Yb[,1]==1,], na.rm=T ), 2 ),
round( colSums( Yb[Yb[,1]==0,], na.rm=T ), 2 )
)

B	<- length(Ybnames)
P	<- length(Ynames)
save( B, P, Ybnames, Ynames, envs, ppsub, qqsub, file='Rdata/setup.Rdata' )
