

library(lmtest)
library(ADGofTest)
library(abind)

# dwtest, bgtest

simARproc <- function(  
  nn = 1000, cf = list( ar=c(0.2,0.1) ),
  dof = Inf
){
  pp <- length(cf);
  if(is.infinite(dof)){
    xx <- arima.sim(
      model=cf, n = nn, rand.gen = rnorm
    )
  } else {
    xx <- arima.sim(
      model=cf, n = nn, rand.gen = rt, df=dof
    )
  }
  return(xx)
}

doDWTest <- function( dat ){
  out = dwtest( dat~1 );
  return(out$p.value);
}

doBGTest <- function( dat, lag=1 ){
  out = bgtest( dat~1, order=lag );
  return(out$p.value);
}

doSOTest <- function( dat, lag=1 ){
  nn    <- length(dat);
  # center & scale
  dat.s <- scale(dat)/sqrt(nn-1);
  # comp diff
  dat.d <- dat.s[1:(nn-lag)]*dat.s[(lag+1):nn]
  # comp t-stat
  st    <- sum(dat.d)
  # comp pv conc
  pv1   <- exp( -(nn-2)*st^2/64 );
  # comp corrected pv
  pr1   <- 32*nn*(nn+2)/( (nn-lag)*(nn-2) )
  c0    <- sqrt( pr1 )*exp(
             lgamma( pr1 ) - lgamma(pr1+0.5)
           );
  pv2   <- c0*pbeta( pv1, pr1, 0.5 );
  return(c(pv1,pv2))
}

runAR1Tests <- function(
  nn = 1000, cf.max = 0.15, inc=0.01, 
  lag= 1,reps=100,offset=0,multLag=5,
  dof=Inf
){
  phi = seq(0,cf.max,inc);
  len = length(phi);
  tot = 1+multLag+2*multLag
  pvs = array( 0, dim=c(len,tot,reps) );
  
  pb = txtProgressBar( min=0,max=reps,style=3 )
  for( rr in 1:reps ){
    for( ll in 1:len ){
      # Generate Data
      if(ll==1){
        dat = simARproc( nn, list(), dof=dof )
      } else {
        par = c( rep(0,offset), phi[ll] )
        dat = simARproc( nn, list(ar=par), dof=dof )
      }
      # Do Tests
      pvs[ ll,1,rr ]   = doDWTest( dat )
      for(ii in 1:multLag)
        pvs[ ll,1+ii,rr ]   = doBGTest( dat, lag=lag+ii-1 )
      for(ii in 1:multLag)
        pvs[ ll,1+multLag+2*(ii-1)+1:2,rr ] = doSOTest( dat, lag=lag+ii-1 )
    }
    setTxtProgressBar( pb, rr );
  }
  close(pb);
  
  return(pvs)
}

plotAROutputs <- function( 
  out, multLag=5
){
  pv.m = apply(out,c(1,2),mean);
  ind1 = c(1,rep(2,multLag),rep(3:4,multLag));
  #ind2 = c(1,1:multLag,rep(1:5,each=multLag));
  ind2 = c(1,rep(3,3*multLag));
  matplot(
    pv.m, col=ind1,lty=ind2,
    type='l', lwd=3,log='y',las=1,ylim=c(10^-3,1)
  )
  legend(
    "bottomleft",legend=c(
      "Durbin-Watson","Breusch-Godfrey",
      "Rotation Conc","Rotation Adj"
    ),
    lty=1,lwd=3,col=1:4
  )
}

plotARPower <- function( 
  out, multLag=5, alf=0.05
){
  phi = seq(0,0.15,0.01);
  ind1 = c(1,rep(2,multLag),rep(3:4,multLag));
  #ind2 = c(1,1:multLag,rep(1:5,each=multLag));
  ind2 = c(1,rep(3,3*multLag));
  pv.m <- (out <= alf)
  pv.m = apply(pv.m,c(1,2),mean);
  matplot(
    phi, pv.m, col=ind1,lty=ind1,
    type='l', lwd=3,las=1,ylim=c(0,1),
    ylab="Power",xlab="AR Parameter"
  )
  legend(
    "topleft",legend=c(
      "Durbin-Watson","Breusch-Godfrey",
      "Rotation Conc","Rotation Adj"
    ),
    lty=1:4,lwd=3,col=1:4
  )
  abline(h=alf,col='gray',lwd=4)
}


runNullTests <- function(
  nn = 1000, reps=100, dof=Inf, multLag=5, lag=1
){
  tot = 1+multLag+2*multLag
  pvs = array( 0, dim=c(tot,reps) );
  
  pb = txtProgressBar( min=0,max=reps,style=3 )
  for( rr in 1:reps ){
      # Generate Data
      dat = simARproc( nn, list(), dof=dof )
      # Do Tests
      pvs[ 1,rr ]   = doDWTest( dat )
      for(ii in 1:multLag)
        pvs[ 1+ii,rr ]   = doBGTest( dat, lag=lag+ii-1 )
      for(ii in 1:multLag)
        pvs[ 1+multLag+2*(ii-1)+1:2,rr ] = doSOTest( dat, lag=lag+ii-1 )
    setTxtProgressBar( pb, rr );
  }
  close(pb);
  
  return(pvs)
}

plotNull <- function( pvs ){
  hist(pvs,main="",xlab="p-values",las=1)
  print(ks.test(pvs,punif))
  print(ad.test(pvs,punif))
}

betaMOM <- function( pvs ){
  mu = mean(pvs)
  vr = var(pvs)
  aa = mu^2*(1-mu)/vr - mu;
  bb = ( mu*(1-mu)/vr - 1 )*(1-mu)
  return( c(aa,bb) )
}



runARMultLagTests <- function(
  nn = 1000, cf.max = 0.15, inc=0.01, 
  lag= 1,reps=100,offset=1,multLag=5,
  dof=Inf
){
  phi = seq(0,cf.max,inc);
  len = length(phi);
  tot = 1+multLag+2*multLag
  pvs = array( 0, dim=c(len,tot,reps) );
  
  pb = txtProgressBar( min=0,max=reps,style=3 )
  for( rr in 1:reps ){
    for( ll in 1:len ){
      # Generate Data
      if(ll==1){
        dat = simARproc( nn, list(), dof=dof )
      } else {
        par = c( 0,phi[ll], rep(0,offset), -phi[ll] )
        dat = simARproc( nn, list(ar=par), dof=dof )
      }
      # Do Tests
      pvs[ ll,1,rr ]   = doDWTest( dat )
      for(ii in 1:multLag)
        pvs[ ll,1+ii,rr ]   = doBGTest( dat, lag=lag+ii-1 )
      for(ii in 1:multLag)
        pvs[ ll,1+multLag+2*(ii-1)+1:2,rr ] = doSOTest( dat, lag=lag+ii-1 )
    }
    setTxtProgressBar( pb, rr );
  }
  close(pb);
  
  return(pvs)
}

plotARPowerMult <- function( 
  out, multLag=5, alf=0.05
){
  phi = seq(0,0.15,0.01);
  ind1 = c(1,rep(2,multLag),rep(3:4,multLag),5);
  #ind2 = c(1,1:multLag,rep(1:5,each=multLag));
  ind2 = c(1,rep(3,3*multLag));
  fshpv <- pchisq(
    q = array(
      data = -2*(log( out[,10,] ) + log( out[,14,] )),
      dim  = c(16,5000)
    ), df= 4,lower.tail=F
  )
  out  <- abind( out, fshpv, along=2 )
  pv.m <- (out <= alf)
  pv.m = apply(pv.m,c(1,2),mean);
  matplot(
    phi, pv.m, col=ind1,lty=ind1,
    type='l', lwd=3,las=1,ylim=c(0,1),
    ylab="Power",xlab="AR Parameter"
  )
  legend(
    "topleft",legend=c(
      "Durbin-Watson","Breusch-Godfrey",
      "Rotation Conc","Rotation Adj",
      "Rotation Fisher"
    ),
    lty=1:5,lwd=3,col=1:5
  )
  abline(h=alf,col='gray',lwd=4)
}

plotARPowerMult2 <- function( 
  out, multLag=5, alf=0.05
){
  phi = seq(0,0.15,0.01);
  ind1 = c(1,rep(2,multLag),rep(3:4,multLag),rep(5,multLag));
  #ind2 = c(1,1:multLag,rep(1:5,each=multLag));
  ind2 = c(1,rep(3,3*multLag));
  
  out  <- abind( out, -2*log(out[,8,]), along=2 );
  indx <- c(8,10,12,14,16);
  for( ii in 2:multLag ){
    out  <- abind( 
      out, out[,17+ii-2,] - 2*log(out[,indx[ii],]), 
      along=2 
    )
  }
  for( ii in 1:multLag ){
    fshpv <- 
      pchisq(
      q = array(
        data = out[,17+ii-1,],
        dim  = c(16,5000)
      ), df= 2*ii,lower.tail=F
    )
    out[,17+ii-1,]  <- fshpv;
  }
  pv.m <- (out <= alf)
  pv.m = apply(pv.m,c(1,2),mean);
  matplot(
    phi, pv.m, col=ind1,lty=ind1,
    type='l', lwd=3,las=1,ylim=c(0,1),
    ylab="Power",xlab="AR Parameter"
  )
  legend(
    "topleft",legend=c(
      "Durbin-Watson","Breusch-Godfrey",
      "Rotation Conc","Rotation Adj",
      "Rotation Fisher"
    ),
    lty=1:5,lwd=3,col=1:5
  )
  abline(h=alf,col='gray',lwd=4)
}

###########
# JSON // sunspots 
###########

library(jsonlite)

dat <- read_json( 
  "./observed-solar-cycle-indices.json",
  simplifyVector=T
)

xx = dat[[2]];
yy = resid( 
  ar(xx,order.max=1)
);

plotSolar <- function(
  dat, rng, bwd=18
){
  tmp= log2(dat+1)[rng]
  ks = ksmooth(rng,tmp,bandwidth=bwd);
  plot(
    rng,tmp, type='l', xlab="year", 
    ylab="Log Solar Intensity", las=1,
    xaxt='n', lwd=3
  )
  lines(
    ks,col='salmon',lwd=2,lty=3
  )
  axis(
    side=1,at=seq( 1813, 3313, length.out=6 ),
    labels=seq( 1900, 2025, length.out=6 )
  )
}

plotSolarRes <- function(
  dat, rng, bwd=18
){
  tmp= log2(dat+1)[rng]
  ks = ksmooth(rng,tmp,bandwidth=bwd);
  plot(
    rng,tmp-ks$y, xlab="year", 
    ylab="Residual Solar Intensity", las=1,
    xaxt='n', lwd=3
  )
  axis(
    side=1,at=seq( 1813, 3313, length.out=6 ),
    labels=seq( 1900, 2025, length.out=6 )
  )
}

getSolarResid <- function(
  dat, bwd=12, lag.mx=24, rng = 1:3314,
  ARRES=F
){
  tmp= log2(dat+1)
  ks = ksmooth(rng,tmp[rng],bandwidth=bwd);
  pvs= matrix( 0, nrow=lag.mx, ncol=2 );
  df = (tmp[rng]-ks$y)
  if(ARRES){
    df <- resid(ar(df,order.max=1))
  }
  df  <- df[-1]
  rng <- rng[-1]
  for( ii in 1:lag.mx ){
    pvs[ii,1] = doBGTest( df, ii );
    pvs[ii,2] = doSOTest( df, ii )[2];
  }
  matplot( 
    pvs, col=c('red','blue'),lty=2:3,log='y',
    type='b',las=1,lwd=3,pch=1:2, 
    ylab='p-values',xlab='lag'
  )
  abline(h=c(0.05,0.01,0.001),lwd=3,col=c('lightgray','gray','darkgray'))
  text( 
    x=rep(18,3),y=c(0.05,0.01,0.001),
    pos=3,labels=c("5%","1%","0.1%"),offset=0.2
  )
  legend(
    15,0.0001,legend=c( "Breusch-Godfrey","Rotation Adj" ), 
    lwd=3, lty=2:3, pch=1:2, col=c("red","blue")
  )
  return(pvs)
}

