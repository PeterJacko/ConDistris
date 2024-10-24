###################################################################################################
# Define Functions
###################################################################################################
#
# PETER ONE-ARM (BEGIN)
#
#    var  <- function( t1 ){ # variance of the proportion estimator, defined p_hat*(1-p_hat)/n = t*(n-t)/n^3 ; maximum at t=n/2 being 1/(4*n)
#      ifelse( max_var , 1 / ( 4 * n[1] ) , # use the maximum Bernoulli variance
#              max( ( t1 * ( n[1] - t1 ) ) / ( n[1]^3 ) , boundary ) ) # use the boundary if var is too small
#    }
#
# The next ".one" functions are obtained from the 2-arm variant setting t2/n[2]=ref and N=n[2]=infty, using only n[1]
glf.boundary.one <- function( n , s ){ # obtained by using t = s in the formula of estimator variance p_hat*(1-p_hat)/n = t*(n-t)/n^3
  # if( n == 0 ){ return( Inf ) }
  # assuming n > 0, s >= 0
  ifelse( n == 1 | n < 2 * s , 1 / ( 4 * n ) , s * ( n - s ) / ( n^3 ) )
}
#
glf.atinterim.successbyinterimrespondersandinterimsamplesize.one <- function( alpha , beta , Delta , ref , p_eval , s_step , n_interim_step , n , var_option , var_parameter){
    sub.atinterim.successprobability <- function( s , n_interim ){
    if ( n_interim < s ){ # mis-specified s
      return( NaN )
    }
    if ( n_interim > n ){ # mis-specified interim
      return( NaN )
    }
    # assuming n >= n_interim >= s for the rest
    # the next 3 if statements define the variance of the proportion estimator, defined p_hat*(1-p_hat)/n = t*(n-t)/n^3 ; maximum at t=n/2 being 1/(4*n)
    if ( var_option == "boundary" ) { # use the lower boundary if var is too small
      var  <- function( t1 ){ max( ( t1 * ( n - t1 ) ) / ( n^3 ) , boundary ) }
    }
    if ( var_option == "fixed" ) { # use a fixed variance specified by var_parameter
      var  <- function( t1 ){ var_parameter / n }
    }
    if ( var_option == "plus" ) {  # use the variance with additional var_parameter successes and var_parameter failures if less than 10 successes or less than 10 failures
      var  <- function( t1 ){ 
        ifelse( ( t1 < 10 | ( n - t1 < 10 ) ) , ( ( t1 + var_parameter ) * ( n - t1 + var_parameter ) ) / ( ( n + 2 * var_parameter )^3 ) , 
                ( t1 * ( n - t1 ) ) / ( n^3 ) ) 
      }
    }
    phi1 <- function( t1 ){ pnorm( ( ref - t1 / n ) / sqrt( var( t1 ) ) ) } 
    phi2 <- function( t1 ){ pnorm( ( ( ref + Delta ) - t1 / n ) / sqrt( var( t1 ) ) ) }
    con  <- function( t1 ){ ifelse( phi1( t1 ) < alpha , ifelse( phi2( t1 ) < beta , 1 , 0 ) , 0 ) }
    prob <- function( t ){ dbinom( t , n - n_interim , p_eval ) }
    t_vec   <- c( 0 : ( n - n_interim ) )
    sum( sapply( t_vec + s , con ) * sapply( t_vec , prob ) )
  }
  if ( var_option == "boundary" ) { boundary <- glf.boundary.one( n , var_parameter ) }
  s_vec <- seq( 0 , n , by = s_step )
  n_interim_vec <- seq( 0 , n , by = n_interim_step )
  s_vec1 <- expand.grid( s_vec , n_interim_vec )[ 1 ]
  n_interim_vec1 <- expand.grid( s_vec , n_interim_vec )[ 2 ]
  mapply( sub.atinterim.successprobability , t( s_vec1 ) , t( n_interim_vec1 ) )
}
#
glf.atinterim.successbyinterimsamplesize.one <- function( alpha , beta , Delta , ref , p_eval , s , n_interim_step , n , var_option , var_parameter){
  sub.atinterim.successprobability <- function( n_interim ){
    if ( n_interim < s ){ # mis-specified s
      return( NaN )
    }
    if ( n_interim > n ){ # mis-specified interim
      return( NaN )
    }
    # assuming n >= n_interim >= s for the rest
    # the next 3 if statements define the variance of the proportion estimator, defined p_hat*(1-p_hat)/n = t*(n-t)/n^3 ; maximum at t=n/2 being 1/(4*n)
    if ( var_option == "boundary" ) { # use the lower boundary if var is too small
      var  <- function( t1 ){ max( ( t1 * ( n - t1 ) ) / ( n^3 ) , boundary ) }
    }
    if ( var_option == "fixed" ) { # use a fixed variance specified by var_parameter
      var  <- function( t1 ){ var_parameter / n }
    }
    if ( var_option == "plus" ) {  # use the variance with additional var_parameter successes and var_parameter failures if less than 10 successes or less than 10 failures
      var  <- function( t1 ){ 
        ifelse( ( t1 < 10 | ( n - t1 < 10 ) ) , ( ( t1 + var_parameter ) * ( n - t1 + var_parameter ) ) / ( ( n + 2 * var_parameter )^3 ) , 
                ( t1 * ( n - t1 ) ) / ( n^3 ) ) 
      }
    }
    phi1 <- function( t1 ){ pnorm( ( ref - t1 / n ) / sqrt( var( t1 ) ) ) } 
    phi2 <- function( t1 ){ pnorm( ( ( ref + Delta ) - t1 / n ) / sqrt( var( t1 ) ) ) }
    con  <- function( t1 ){ ifelse( phi1( t1 ) < alpha , ifelse( phi2( t1 ) < beta , 1 , 0 ) , 0 ) }
    prob <- function( t ){ dbinom( t , n - n_interim , p_eval ) }
    t_vec   <- c( 0 : ( n - n_interim ) )
    sum( sapply( t_vec + s , con ) * sapply( t_vec , prob ) )
  }
  if ( var_option == "boundary" ) { boundary <- glf.boundary.one( n , var_parameter ) }
  n_interim_vec <- seq( 0 , n , by = n_interim_step )
  sapply( n_interim_vec , sub.atinterim.successprobability )
}
#
glf.atinterim.successbyinterimresponders.one <- function( alpha , beta , Delta , ref , p_eval , s_step , n_interim , n , var_option , var_parameter){
  sub.atinterim.successprobability <- function( s ){
    if ( n_interim < s ){ # mis-specified s
      return( NaN )
    }
    if ( n_interim > n ){ # mis-specified interim
      return( NaN )
    }
    # assuming n >= n_interim >= s for the rest
    # the next 3 if statements define the variance of the proportion estimator, defined p_hat*(1-p_hat)/n = t*(n-t)/n^3 ; maximum at t=n/2 being 1/(4*n)
    if ( var_option == "boundary" ) { # use the lower boundary if var is too small
      var  <- function( t1 ){ max( ( t1 * ( n - t1 ) ) / ( n^3 ) , boundary ) }
    }
    if ( var_option == "fixed" ) { # use a fixed variance specified by var_parameter
      var  <- function( t1 ){ var_parameter / n }
    }
    if ( var_option == "plus" ) {  # use the variance with additional var_parameter successes and var_parameter failures if less than 10 successes or less than 10 failures
      var  <- function( t1 ){ 
        ifelse( ( t1 < 10 | ( n - t1 < 10 ) ) , ( ( t1 + var_parameter ) * ( n - t1 + var_parameter ) ) / ( ( n + 2 * var_parameter )^3 ) , 
                ( t1 * ( n - t1 ) ) / ( n^3 ) ) 
      }
    }
    phi1 <- function( t1 ){ pnorm( ( ref - t1 / n ) / sqrt( var( t1 ) ) ) } 
    phi2 <- function( t1 ){ pnorm( ( ( ref + Delta ) - t1 / n ) / sqrt( var( t1 ) ) ) }
    con  <- function( t1 ){ ifelse( phi1( t1 ) < alpha , ifelse( phi2( t1 ) < beta , 1 , 0 ) , 0 ) }
    prob <- function( t ){ dbinom( t , n - n_interim , p_eval ) }
    t_vec   <- c( 0 : ( n - n_interim ) )
    sum( sapply( t_vec + s , con ) * sapply( t_vec , prob ) )
  }
  if ( var_option == "boundary" ) { boundary <- glf.boundary.one( n , var_parameter ) }
  s_vec <- seq( 0 , n_interim , by = s_step )
  sapply( s_vec , sub.atinterim.successprobability )
}
#
glf.atinterim.successbysamplesize.one <- function( alpha , beta , Delta , ref , p_eval , s , n_interim , n_max , n_step , var_option , var_parameter){
  sub.atinterim.successprobability <- function( n ){
    if ( n_interim < s ){ # mis-specified s
      return( NaN )
    }
    if ( n_interim > n ){ # mis-specified interim
      return( NaN )
    }
    # assuming n >= n_interim >= s for the rest
    # the next 3 if statements define the variance of the proportion estimator, defined p_hat*(1-p_hat)/n = t*(n-t)/n^3 ; maximum at t=n/2 being 1/(4*n)
    if ( var_option == "boundary" ) { # use the lower boundary if var is too small
      var  <- function( t1 ){ max( ( t1 * ( n - t1 ) ) / ( n^3 ) , boundary ) }
    }
    if ( var_option == "fixed" ) { # use a fixed variance specified by var_parameter
      var  <- function( t1 ){ var_parameter / n }
    }
    if ( var_option == "plus" ) {  # use the variance with additional var_parameter successes and var_parameter failures if less than 10 successes or less than 10 failures
      var  <- function( t1 ){ 
        ifelse( ( t1 < 10 | ( n - t1 < 10 ) ) , ( ( t1 + var_parameter ) * ( n - t1 + var_parameter ) ) / ( ( n + 2 * var_parameter )^3 ) , 
                ( t1 * ( n - t1 ) ) / ( n^3 ) ) 
      }
    }
    phi1 <- function( t1 ){ pnorm( ( ref - t1 / n ) / sqrt( var( t1 ) ) ) } 
    phi2 <- function( t1 ){ pnorm( ( ( ref + Delta ) - t1 / n ) / sqrt( var( t1 ) ) ) }
    con  <- function( t1 ){ ifelse( phi1( t1 ) < alpha , ifelse( phi2( t1 ) < beta , 1 , 0 ) , 0 ) }
    prob <- function( t ){ dbinom( t , n - n_interim , p_eval ) }
    t_vec   <- c( 0 : ( n - n_interim ) )
    sum( sapply( t_vec + s , con ) * sapply( t_vec , prob ) )
  }
  if ( var_option == "boundary" ) { boundary <- glf.boundary.one( n , var_parameter ) }
  n_vec <- seq( n_interim , n_max , by = n_step )
  sapply( n_vec , sub.atinterim.successprobability )
}
#
glf.successbyparameter.one <- function( alpha , beta , Delta , ref , step , n , var_option , var_parameter){
  sub.successprobability <- function( p1 ){
    if ( n == 0 ){ # to avoid NaN in t1/n (and infinity var)
      return( 0 )
    }
    # assuming n > 0 for the rest
    # the next 3 if statements define the variance of the proportion estimator, defined p_hat*(1-p_hat)/n = t*(n-t)/n^3 ; maximum at t=n/2 being 1/(4*n)
    if ( var_option == "boundary" ) { # use the lower boundary if var is too small
      var  <- function( t1 ){ max( ( t1 * ( n - t1 ) ) / ( n^3 ) , boundary ) }
    }
    if ( var_option == "fixed" ) { # use a fixed variance specified by var_parameter
      var  <- function( t1 ){ var_parameter / n }
    }
    if ( var_option == "plus" ) {  # use the variance with additional var_parameter successes and var_parameter failures if less than 10 successes or less than 10 failures
      var  <- function( t1 ){ 
        ifelse( ( t1 < 10 | ( n - t1 < 10 ) ) , ( ( t1 + var_parameter ) * ( n - t1 + var_parameter ) ) / ( ( n + 2 * var_parameter )^3 ) , 
                ( t1 * ( n - t1 ) ) / ( n^3 ) ) 
      }
    }
    phi1 <- function( t1 ){ pnorm( ( ref - t1 / n ) / sqrt( var( t1 ) ) ) } 
    phi2 <- function( t1 ){ pnorm( ( ( ref + Delta ) - t1 / n ) / sqrt( var( t1 ) ) ) }
    con  <- function( t1 ){ ifelse( phi1( t1 ) < alpha , ifelse( phi2( t1 ) < beta , 1 , 0 ) , 0 ) }
    prob <- function( t1 ){ dbinom( t1 , n , p1 ) }
    t1_vec   <- c( 0 : n )
    sum( sapply( t1_vec , con ) * sapply( t1_vec , prob ) )
  }
  if ( var_option == "boundary" ) { boundary <- glf.boundary.one( n , var_parameter ) }
  p1_vec <- seq( ref , 1 , by = step )
  sapply( p1_vec , sub.successprobability )
}
#
glf.successbysamplesize.one <- function(alpha,beta,Delta,p,maxss,step,var_option,var_parameter){
  sub.successprobability <- function( n ){
    if ( n == 0 ){ # to avoid NaN in t1/n (and infinity var)
      return( 0 )
    }
    # assuming n > 0 for the rest
    if ( var_option == "boundary" ) { boundary <- glf.boundary.one( n , var_parameter ) }
    # the next 3 if statements define the variance of the proportion estimator, defined p_hat*(1-p_hat)/n = t*(n-t)/n^3 ; maximum at t=n/2 being 1/(4*n)
    if ( var_option == "boundary" ) { # use the lower boundary if var is too small
      var  <- function( t1 ){ max( ( t1 * ( n - t1 ) ) / ( n^3 ) , boundary ) }
    }
    if ( var_option == "fixed" ) { # use a fixed variance specified by var_parameter
      var  <- function( t1 ){ var_parameter / n }
    }
    if ( var_option == "plus" ) {  # use the variance with additional var_parameter successes and var_parameter failures if less than 10 successes or less than 10 failures
      var  <- function( t1 ){ 
        ifelse( ( t1 < 10 | ( n - t1 < 10 ) ) , ( ( t1 + var_parameter ) * ( n - t1 + var_parameter ) ) / ( ( n + 2 * var_parameter )^3 ) , 
                ( t1 * ( n - t1 ) ) / ( n^3 ) ) 
      }
    }
    phi1 <- function( t1 ){ pnorm( ( ref - t1 / n ) / sqrt( var( t1 ) ) ) }
    phi2 <- function( t1 ){ pnorm( ( ( ref + Delta ) - t1 / n ) / sqrt( var( t1 ) ) ) }
    con  <- function( t1 ){ ifelse( phi1( t1 ) < alpha , ifelse( phi2( t1 ) < beta , 1 , 0 ) , 0 ) }
    prob <- function( t1 ){ dbinom( t1 , n , p1 ) }
    t1_vec   <- c( 0 : n )
    sum( sapply( t1_vec , con ) * sapply( t1_vec , prob ) )
  }
  n_vec <- seq( 1 , maxss , by = step )
  ref <- p[1]
  p1 <- p[2]
  sapply( n_vec , sub.successprobability )
}
#
glf.failurebyparameter.one <- function(gamma,ref,step,n,var_option,var_parameter){
  sub.failureprobability <- function( p1 ){
    if ( n == 0 ){ # to avoid NaN in t1/n (and infinity var)
      return( 0 )
    }
    # assuming n > 0 for the rest
    # the next 3 if statements define the variance of the proportion estimator, defined p_hat*(1-p_hat)/n = t*(n-t)/n^3 ; maximum at t=n/2 being 1/(4*n)
    if ( var_option == "boundary" ) { # use the lower boundary if var is too small
      var  <- function( t1 ){ max( ( t1 * ( n - t1 ) ) / ( n^3 ) , boundary ) }
    }
    if ( var_option == "fixed" ) { # use a fixed variance specified by var_parameter
      var  <- function( t1 ){ var_parameter / n }
    }
    if ( var_option == "plus" ) {  # use the variance with additional var_parameter successes and var_parameter failures if less than 10 successes or less than 10 failures
      var  <- function( t1 ){ 
        ifelse( ( t1 < 10 | ( n - t1 < 10 ) ) , ( ( t1 + var_parameter ) * ( n - t1 + var_parameter ) ) / ( ( n + 2 * var_parameter )^3 ) , 
                ( t1 * ( n - t1 ) ) / ( n^3 ) ) 
      }
    }
    phi1 <- function( t1 ){ pnorm( ( ref - t1 / n ) / sqrt( var( t1 ) ) ) }
    con  <- function( t1 ){ ifelse( phi1( t1 ) > gamma , 1 , 0 ) }
    prob <- function( t1 ){ dbinom( t1 , n , p1 ) }
    t1_vec   <- c( 0 : n )
    sum( sapply( t1_vec , con ) * sapply( t1_vec , prob ) )
  }
  if ( var_option == "boundary" ) { boundary <- glf.boundary.one( n , var_parameter ) }
  p1_vec <- seq( ref , 1 , by = step )
  sapply( p1_vec , sub.failureprobability )
}
#
# PETER ONE-ARM (END)
#
#
#
###################################################################################################
