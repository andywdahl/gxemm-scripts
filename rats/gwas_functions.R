create_X  <- function( y, phen, all_Xs, model_menu_file='data/model_menu_final.txt', rm_sel=TRUE, noadj ){

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
