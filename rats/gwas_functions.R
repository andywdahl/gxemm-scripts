#library(emma)
#library(phenix)

################ gwas-running functions
gwas  <- function( y, all_Xs, X, K, phen, chrs=1:21, all_fits=FALSE, rm_sel=TRUE, noadj=FALSE ){

  ##### create covariates, permute
  if( missing( X ) )
    X     <- create_X( y=y, all_Xs=all_Xs, phen=phen, rm_sel=rm_sel, noadj=noadj )   ### create covariate matrix

  sub   <- which( apply( !is.na( cbind( y, X ) ), 1, all ) )
  y     <- y[sub]
  X     <- X[sub,]
  K     <- K[sub,sub]

  pars  <- fit.pars( y, X, K )  #### transform to use OLS

  ##### sanity checking
  if( ! all(rownames(X)==names(y)) | ! all(rownames(X)==rownames(K)) )
    stop('bad names')

  # now test for association with all happy intervals (haplotypes)
  pvals     <- c()
  if( all_fits )
    fit_list  <- list()
  for( chr in paste0('chr',c(1:20,'X'))[ chrs ] ){
    print( chr )
    chr_dir     <- paste('/data/eider/dahl/rat_gwas/raw_data/sdir/additive/',chr,sep='')
    files       <- list.files(chr_dir)
    files       <- files[grep('^@',files)]
    markers     <- unlist(lapply(paste(chr_dir,'/',files,sep=''),load,.GlobalEnv))
    outs        <- lapply(X=mget(markers,envir=.GlobalEnv),FUN=my_anova,trans=pars$trans,y_rot=pars$y_rot,x_rot=pars$x_rot,chr=chr,sub=sub)
    pvals_loc   <- sapply( outs, function(out) out$p )
    names(pvals_loc)  <- markers
    pvals             <- c( pvals, pvals_loc )
    if( all_fits ){
      fits_loc          <- lapply( outs, function(out) out$fits )
      names(fits_loc)   <- markers
      fit_list          <- c( fit_list, fits_loc )
    }
    rm(list=markers,envir=.GlobalEnv)
  }
  out <- list( pvals=pvals, h2=pars$h2 )
  if( all_fits )
    out <- c( out, list( fit_list=fit_list ) )
  return( out )

}

create_X  <- function( y, phen, all_Xs, model_menu_file='raw_data/model_menu_final.txt', rm_sel=TRUE, noadj ){

  ##### load ingredients
  model_menu  <- read.delim(model_menu_file,header=T,colClasses='character',na.string='') #### p2covariate info

  ##### remove unreliable measurements
  p       <- which( model_menu[,2] == phen )
  if( length(p) != 1 ){
    warning( 'no phen found; using all' )
    covs  <- rev( sort( unique( unlist( sapply( 1:nrow(model_menu), function(p){
      strsplit(model_menu[p,'covariates'],',')[[1]] 
    }) ) ) ) )
   add.BW <- TRUE 
  } else {
    covs    <- strsplit(model_menu[p,'covariates'],',')[[1]]
    add.BW <- FALSE
  }

  if( any( covs %in% c( 'test_worked','reliable' ) ) ){
    sel         <- which( covs %in% c( 'test_worked','reliable' ) )
    if( length(sel) == 1 ){
      good.rows   <- which( as.logical(all_Xs[ ,covs[sel] ]) )
    } else if( length(sel) == 2 ){
      if( rm_sel ){
        warning( 'Forcing rm_sel=FALSE' )
        rm_sel  <- FALSE
      }
    } else {
      stop( 'sel is weird length' )
    }
    covs        <- covs[-sel]
    if( rm_sel ){
      if( length( good.rows ) == 0 )
        stop( 'No good rows')
      all_Xs[ -good.rows, ]  <- NA
    }
  }

  X           <- matrix(1,nrow=nrow(all_Xs),ncol=1,dimnames=list( rownames(all_Xs), 'Intercept' ) )
  if (covs[1]!='sex')
    stop('first cov isnt sex')

  ##### encode as model matrix
  for( cov in covs ){
    if( cov %in% c( 'sex' ,'batch','is.albino','facs_day' ) ){
      x             <- as.factor( all_Xs[,cov] )
      contrasts(x)  <- contrasts(x,contrasts=FALSE)[,-length(levels(x))]
      X             <- cbind( X, model.matrix(~x+1, model.frame(~x+1,data.frame(x=x), na.action=function(x)x))[,-1] )
    } else if( cov %in% c( 'Haemalysis','BW_at_day28_pi','age','BW_at_IPGTT','BW_at_day9_pi' ) ){
      X             <- cbind( X, matrix( as.numeric(all_Xs[,cov]), nrow(all_Xs), 1, dimnames=list(rownames(all_Xs),cov) ) )
    } else {
      stop( c( 'unexpected covariate:', cov ) )
    }
  }
  if( add.BW ){
    if( noadj ){
      warning( 'Not adding BW' )
    } else {
      warning( 'Tacking BW measurement on to covariates' )
      cov <- 'BW_at_day_immunization_bc'
      X   <- cbind( X, matrix( as.numeric(all_Xs[,cov]), nrow(all_Xs), 1, dimnames=list(rownames(all_Xs),cov) ) )
    }
  }

  ###### check for singularities
  while( ( nsing <- sum( svd(X[ which( rowSums( is.na( cbind( y, X ) ) ) == 0 ), ] )$d < 1e-7 ) ) > 0 ){
    print( 'singular X; dropping a level' )
    for( i in sort.list( colMeans( is.na(X) ), dec=TRUE ) ){
      X.tmp     <- X[ , -i ]
      nsing.tmp <- sum( svd(X.tmp[ which( rowSums( is.na( cbind( y, X.tmp ) ) ) == 0 ), ] )$d < 1e-6 )
      if( nsing.tmp < nsing ){
        X <- X.tmp
        break
      }
      if( i == sort.list( colMeans( is.na(X) ), dec=TRUE )[ ncol(X) ] )
        stop( 'Crap covariate stuff' )
    }
  }
  X
}

fit.pars  <- function( y, X, K ){
  emma        <- emma.REMLE( y=y, X=X, K=K )    ##### null VCs
  trans       <- mat.sqrt( emma$vg*K+emma$ve*diag(nrow(X)), inv=TRUE )
  y_rot       <- trans %*% y
  x_rot       <- trans %*% X
  beta        <- coef( lm(y_rot   ~ -1+x_rot) )
  h2          <- ( emma$vg*tr(K) )/( emma$vg*tr(K) + emma$ve*nrow(K) )
  return(list( sigma_g=emma$vg, sigma_e=emma$ve, beta=beta, y_rot=y_rot, x_rot=x_rot, h2=h2, trans=trans ))
}

my_anova  <- function(element,trans,y_rot,x_rot,chr,sub) {
  element   <- trans %*% (element[sub,])  ### sub accounts for NA covariates
  fit0=lm(y_rot   ~ -1+x_rot)
  if( chr!='X' ){
    fit1=lm(y_rot ~ -1+x_rot+element)
  } else {
    fit1=lm(y_rot ~ -1+x_rot+element+x_rot[,2]:element)
  }
  p <- anova(fit0,fit1)[,'Pr(>F)'][2]
  return( list( p=p, fits=list( fit0=fit0, fit1=fit1 ) ) )
}

################################################
################ post-processing functions
################################################
dist  <- function( x, y ){                ### computes genomic distance between x and y, lists containing p=phenotype, bp and chr
  if( is.na( x$phen ) | is.na( y$phen ) ) ### phenotypes undefined
    stop( 'NA phenotype' )
  if( as.character( x$chr ) != as.character( y$chr ) | as.character( x$phen ) != as.character( y$phen ) ) ### phenotype or chr disagree
    return( Inf )
  ifelse(
    ( x$bpi >= y$bpi & x$bpi <= y$bpf ) |   ### x0 in y
    ( x$bpf >= y$bpi & x$bpf <= y$bpf ) |   ### xf in y
    ( y$bpi >= x$bpi & y$bpi <= x$bpf ) |   ### y0 in x
    ( y$bpf >= x$bpi & y$bpf <= x$bpf ),    ### yf in x
    0,                                    ### windows overlap
    min(c( abs( x$bpi - y$bpf ), abs( x$bpf - y$bpi ) ))
  )
}

df_dist <- function( hit, df )
  sapply( 1:nrow(df), function(i) dist( x=hit, y=df[i,] ) )

df_trim <- function( df_in, width ){
  df_out  <- df_in[ 1, ]
  for( i in 2:nrow(df_in) ){                ### if new, append to df_out; otherwise update relevant row(s) of df_out
    overlap  <- which( df_dist( df_in[i,], df_out ) <= width )
    if( length( overlap ) > 0 ){            #### locus already represented
      if( length( overlap ) > 1 ){
        print( 'Multiple hit matches!!!' )  ### this is not necessarily an error...just probably
        df_out[overlap[1],'bpi'] <- min( df_out[overlap,'bpi'],  df_in[i,'bpi'] )
        df_out[overlap[1],'bpf'] <- max( df_out[overlap,'bpf'],  df_in[i,'bpf'] )
        df_out[overlap[1],'logP']<- max( df_out[overlap,'logP'], df_in[i,'logP'], na.rm=TRUE )
        df_out  <- df_out[ -overlap[-1], ]
      } else {
        df_out[overlap,'bpi'] <- min( df_out[overlap,'bpi'],  df_in[i,'bpi'] )
        df_out[overlap,'bpf'] <- max( df_out[overlap,'bpf'],  df_in[i,'bpf'] )
        df_out[overlap,'logP']<- max( df_out[overlap,'logP'], df_in[i,'logP'], na.rm=TRUE )
      }
    } else {                                #### add new row to df_out
      window  <- which( df_dist( df_in[i,], df_in )<= width )
      df_out  <- rbind( df_out, data.frame(
        phen=unique(  df_in$phen[window] ),
        chr =unique(  df_in$chr [window] ),
        bpi =min(     df_in$bpi [window] ),
        bpf =max(     df_in$bpf [window] ),
        logP=max(     df_in$logP[window], na.rm=TRUE )
      ))
    }
  }
  df_out
}

pval2df <- function( pval, thresh, phen, index2chr, index2bp, width=3 ){

  #### snps above threshold
  hit_inds  <- which( -log10(pval) > thresh )
  if( length( hit_inds ) ==  0 )
    return( data.frame( phen=NULL, chr=NULL, bpi=NULL, bpf=NULL, logP=NULL ) )

  #### create data frame out of hit snps
  pval_df   <- data.frame( 
    phen=rep( phen, length(hit_inds) ),
    chr =index2chr    [hit_inds],
    bpi =index2bp     [hit_inds],
    bpf =index2bp     [hit_inds],
    logP=-log10( pval [hit_inds] )
  )

  #### trim for redundancies
  if( nrow( pval_df ) > 1 )
    pval_df   <- df_trim( pval_df, width=width )
  pval_df

}

hits2all_hits <- function( df1, df2, width=0, p1s, p2s, bp, chr ){
  all_hits  <- df_trim( rbind( df1, df2 ), width=width )  ### loses p-value information
  p1        <- numeric( nrow(all_hits) )
  p2        <- numeric( nrow(all_hits) )
  for( i in 1:nrow(all_hits) ){
    hit     <- all_hits[i,]
    phen    <- as.character( hit$phen )
    window  <- which( chr == as.character(hit$chr) & bp/10^6 >= hit$bpi & bp/10^6 <= hit$bpf )
    p1[i]   <- -log10( min( p1s[ phen, window ], na.rm=TRUE ) )
    p2[i]   <- -log10( min( p2s[ phen, window ], na.rm=TRUE ) )
  }

  return( final_df  <- data.frame( 
    phen  =all_hits$phen,
    chr   =all_hits$chr,
    bpi   =all_hits$bpi,
    bpf   =all_hits$bpf,
    logP1 =p1,
    logP2 =p2
  ) )
}

trim_badlocs  <- function( df, width, bad_df ){
  bads  <- list()
  if( length( df ) == 0 )
    return( df )
  for( i in 1:nrow( bad_df ) ){
    bad_reg   <- bad_df[i,]
    bads[[i]] <- which( sapply( 1:nrow(df), function(j) dist( c( bad_reg, list( phen=df$phen[j] ) ), df[j,] ) ) <= width )
  }
  if( any( sapply( bads, length ) > 0 ) )
    df  <- df[ -unique(unlist(bads)), ]
  df
}

#################################################
######### plotting functions ####################
#################################################
lm.text <- function(x,y,x0,y0,col,cex=.7,sigdig=2,trace=TRUE,weirdtol=NULL,labs){

  if( is.null( weirdtol ) | length( sub <- which( abs(x-y) > weirdtol ) ) == 0 ){
    points( x, y, pch=16, col=col, cex=.5 )
  } else {
    points( x[-sub],y[-sub],pch=16, col=col, cex=.5 )
    text(   x[sub], y[sub], pch=16, col=col, cex=.5, labs[sub] )
  }

  fit <- lm( y ~ x )
  abline( fit, col=col )
  betas <- coef( fit )
  text( x0, y0, paste0( 'y=', round(betas[1],sigdig), '+', round(betas[2],sigdig), 'x' ), col=col, cex=cex )
  if( trace ){
    print( dimnames(acc)[[2]] [col] )
    print( summary( lm( y-x ~ x ) ) )
  }
}



################ extract phenix info
U_fxn  <- function( out, logscale=F, h2=F, cutoff ){

  S     <- out$vb_out$vb_pars$mu_S
  beta  <- out$vb_out$vb_pars$mu_beta
  
  sample_comp_cov <- 1/length( out$U ) * ( t(S) %*% S ) * ( beta %*% t(beta) )
  sing_vals       <- sqrt( diag( sample_comp_cov ) )
  if( logscale )
    sing_vals <- log10( sing_vals )
  if( missing( cutoff ) ){
    if( h2 ){
      sing_vals <- ( out$h2 )[ sort.list( sing_vals ) ]
      sing_vals <- sapply( 1:length(sing_vals), function(i){
        len   <- 10
        inds  <- intersect( 1:length(sing_vals), (i-len) + 0:(2*len) )
        mean( sing_vals[ inds ] )
      })
    } else {
      sing_vals <- sort( sing_vals, dec=T )
    }
  } else {
    print( paste( 'cutoff=', cutoff ) )
    M   <- sum( sing_vals > cutoff )
    print( paste( 'P=', ncol(beta) ) )
    print( paste( 'M=', M ) )
    S   <- S[   , sort.list( sing_vals, dec=T )[1:M]   ]
    beta<- beta[  sort.list( sing_vals, dec=T )[1:M],  ]
    sing_vals <- sort( sing_vals, dec=T )[1:M]
  }
  return( list( d=sing_vals, S=S, beta=beta ) )

}
