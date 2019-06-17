library(grid)
library(ggplot2)
library(scales)

gxemm_violin	<- function(
	h2s,
	allps,
	cols=c( 1, 'orange', 2, 3, 1, 'grey', 1, 2, 3, 5, NA, NA ),
	limits=range(h2s),
	ypos=c( -1.25,-1.60,-1.95,-2.3,1.6)+.05,
	xlabs	= expression( h[GREML]^2, h[iid]^2, h[male]^2, h[female]^2, '', h[hom]^2, h[het]^2, '', sigma[g]^2 , v[ma], v[fe], w[ma]-w[fe] ),
	plot.margin=unit(c(1.0, 0.3, 3.6, 0.6), "cm"),
	ybreaks=c(-.25,0,.25,.5,.75,1,1.25)
){

	if( length( bads	<- which( h2s < limits[1] ) ) > 0 ){
		cat( 'Warning! ', colSums( h2s < limits[1], na.rm=T ), ' h2s above ymax; truncating\n' )
		h2s[ bads ]	<- limits[1]+1e-5
	}
	if( length( bads	<- which( h2s > limits[2] ) ) > 0 ){
		cat( 'Warning! ', colSums( h2s > limits[2], na.rm=T ), ' h2s above ymax; truncating\n' )
		h2s[ bads ]	<- limits[2]-1e-5
	}

	P		<- nrow(h2s)
	y		<- as.numeric( t(h2s) )
	x		<- as.factor( matrix( rep( 1:ncol(h2s), each=P ), ncol(h2s), P, byrow=T ) )
	cols<- c(         matrix( rep( cols				, each=P ), ncol(h2s), P, byrow=T ) )

	p	<- ggplot(data.frame( x=x, y=y ), aes(x=x, y=y,fill=x)) +
		scale_fill_manual(values=cols)+
		geom_violin(trim=TRUE,scale='width') +
		theme(
			legend.position="none",
			plot.margin=plot.margin,
			panel.background=element_rect(fill = "transparent",colour = NA),
			plot.background	=element_rect(fill = "transparent",colour = NA),
			axis.title.y = element_text(size=15),
			axis.ticks = element_blank()
		) +
		scale_y_continuous(labels = percent, name='Variance Explained  ',breaks=ybreaks,limits=limits) +
		scale_x_discrete( name='', labels=xlabs ) +
		theme( axis.text.y  = element_text(size=8) ) +
		theme( axis.text.x  = element_text(size=12) ) 

	for( xx in (0:6-1)/4 )
		p	<- p + geom_hline(yintercept=xx, linetype="dotted", color = "grey", size=0.5)
	for( xx in c(0,1) )
		p	<- p + geom_hline(yintercept=xx, linetype="solid" , color = 1     , size=0.7)
		
	gp	<- gpar(cex = 1.5)
	mains		<- c( 'Heritabilities', 'IID', 'Free' )
	mainpos	<- c(2.5,6.5,10.5)
	for( xx in 1:3 )
	p		<- p+annotation_custom( grob = textGrob(label=mains[xx], gp=gp), ymin=ypos[5], ymax=ypos[5], xmin=mainpos[xx], xmax=mainpos[xx] )

	gp	<- gpar(cex = .8)
	labs	<- c( 'Averages', 'p<.05', 'p<.05/P', paste( 'P= # traits =',P) )
	for( kk in 1:4 )
	p		<- p+ annotation_custom( grob = textGrob(label=labs[kk], gp=gp), ymin = ypos[kk], ymax = ypos[kk], xmin=-.3, xmax=ifelse(kk!=4,-.3,2) )

	xvals	<- c( 1:4, 6:7, 9:12 )
	pvec	<- format( colMeans(h2s[,xvals]), digit=2 ) #, nsmall=3
	for( kk in 1:length(xvals) )
	p		<- p+annotation_custom( grob = textGrob(label=pvec[kk], gp=gp), ymin=ypos[1], ymax=ypos[1], xmin=xvals[kk], xmax=xvals[kk]	)

	xvals	<- c( 1:4, 7, 10, 11, 12 )
	pvec	<- colSums( allps < .05, na.rm=T )
	#pvec	<- round( apply( h2s, 2, function(x) Waldtest( mean(x), var(x)/length(x) ) ), 4 )
	for( kk in 1:length(xvals) )
	p		<- p+annotation_custom( grob = textGrob(label=pvec[kk], gp=gp), ymin=ypos[2], ymax=ypos[2], xmin=xvals[kk], xmax=xvals[kk]	)

	pvec	<- colSums( allps < .05/P, na.rm=T )
	for( kk in 1:length(xvals) )
	p		<- p+annotation_custom( grob = textGrob(label=pvec[kk], gp=gp), ymin=ypos[3], ymax=ypos[3], xmin=xvals[kk], xmax=xvals[kk]	)

	gt <- ggplot_gtable(ggplot_build(p))
	gt$layout$clip[gt$layout$name == "panel"] <- "off"
	grid.draw(gt)

}
