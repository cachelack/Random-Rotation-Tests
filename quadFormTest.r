
# > xx.comp1 = runAR1Tests(COMPSIM=T,offset=0,BB=100,multLag=3,n=100,reps=1000,cf.max=0.4,inc=0.1); xx.comp3 = runAR1Tests(COMPSIM=T,offset=2,BB=100,multLag=3,n=100,reps=1000,cf.max=0.5,inc=0.5/4)


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

simMAproc <- function(  
    nn = 1000, cf = list( ma=c(0.2,0.1) ),
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

doLBTest <- function( dat, lag=1 ){
  out = Box.test(
    x=dat, lag=lag, type=c("Ljung-Box")
  )
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

doSnTest <- function( dat, lag=1 ){
  nn    <- length(dat);
  # center & scale
  dat.s <- scale(dat)/sqrt(nn-1);
  # comp diff
  dat.d <- dat.s[1:(nn-lag)]*dat.s[(lag+1):nn]
  # comp t-stat
  st    <- sum(dat.d)
  # comp pv conc
  pv1   <- exp( -st^2/48 );
  # comp corrected pv
  pr1   <- 24*nn*(nn+2)/( (nn-lag) )
  c0    <- sqrt( pr1 )*exp(
    lgamma( pr1 ) - lgamma(pr1+0.5)
  );
  pv2   <- c0*pbeta( pv1, pr1, 0.5 );
  return(c(pv2))
}

doBnTest <- function( dat, lag=1 ){
  nn    <- length(dat);
  # center & scale
  dat.s <- scale(dat)/sqrt(nn-1);
  # comp diff
  dat.d <- dat.s[1:(nn-lag)]*dat.s[(lag+1):nn]
  # comp t-stat
  st    <- sum(dat.d)
  # comp pv conc
  pv1   <- exp( -st^2/64 );
  # comp corrected pv
  pr1   <- 32*nn*(nn+2)/( (nn-lag) )
  c0    <- sqrt( pr1 )*exp(
    lgamma( pr1 ) - lgamma(pr1+0.5)
  );
  pv2   <- c0*pbeta( pv1, pr1, 0.5 );
  return(c(pv2))
}


doRandTest <- function( dat, lag=1, BB=5000, type="PERM" ){
  nn     <- length(dat);
  # center & scale
  dat.s  <- scale(dat)/sqrt(nn-1);
  # comp diff
  dat.d  <- dat.s[1:(nn-lag)]*dat.s[(lag+1):nn]
  # comp t-stat
  st0    <- sum(dat.d)
  st.rnd <- rep(0,BB);
  for( bb in 1:BB ){
    if(type=="ROT"){
      mat = matrix( rnorm(nn^2),nn,nn );
      tmp = qr.Q(qr(mat));
      dat.rnd <- tmp%*%dat.s;
    } else {
      dat.rnd <- switch( type,
        PERM = sample(dat.s),
        SIGN = dat.s*rsign(nn),
        BOTH = sample(dat.s)*rsign(nn)
      )
    }
    dat.d   <- dat.rnd[1:(nn-lag)]*dat.rnd[(lag+1):nn]
    st.rnd[bb] <- sum(dat.d)
  }
  pv.rnd = (sum( st0 < st.rnd ) + 1)/( BB + 1 )
  #return(list(st0,st.rnd))
  return(pv.rnd)
}

rsign <- function( nn ){
  return(
    2*rbinom(nn,1,0.5)-1
  )
}

runAR1Tests <- function(
  nn = 1000, cf.max = 0.15, inc=0.01, 
  lag= 1,reps=100,offset=0,multLag=5,
  dof=Inf, BB=100, HET=F, COMPSIM=F
){
  phi = seq(0,cf.max,inc);
  len = length(phi);
  tot = 1+multLag+2*multLag+3*multLag+multLag
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
      if(HET)
        dat <- dat*seq(1,5,length.out=nn)
      # Do Tests
      # Durbin-Watson
      pvs[ ll,1,rr ]   = doDWTest( dat )
      # Breusch-Godfrey
      for(ii in 1:multLag)
        pvs[ ll,1+ii,rr ]   = doBGTest( dat, lag=lag+ii-1 )
      # SO(n)
      for(ii in 1:multLag)
        pvs[ ll,1+multLag+2*(ii-1)+1:2,rr ] = doSOTest( dat, lag=lag+ii-1 )
      if(COMPSIM){
        for(ii in 1:multLag)
          pvs[ ll,1+3*multLag+ii,rr ] = doRandTest( dat, lag=lag+ii-1, BB,"PERM" )
        for(ii in 1:multLag)
          pvs[ ll,1+4*multLag+ii,rr ] = doRandTest( dat, lag=lag+ii-1, BB,"SIGN" )
        for(ii in 1:multLag)
          pvs[ ll,1+5*multLag+ii,rr ] = doRandTest( dat, lag=lag+ii-1, BB,"ROT" )
      } else {
        for(ii in 1:multLag)
          pvs[ ll,1+3*multLag+ii,rr ] = doSnTest( dat, lag=lag+ii-1 )
        for(ii in 1:multLag)
          pvs[ ll,1+4*multLag+ii,rr ] = doBnTest( dat, lag=lag+ii-1 )
      }
      for(ii in 1:multLag)
        pvs[ ll,1+6*multLag+ii,rr ] = doLBTest( dat, lag=lag+ii-1 )
      
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
  ind1 = c(
    1,rep(2,multLag),rep(3:4,multLag),
    rep(5,multLag),rep(6,multLag),rep(7,multLag)
  );
  #ind2 = c(1,1:multLag,rep(1:5,each=multLag));
  ind2 = c(1,rep(3,6*multLag));
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
  out, multLag=5, alf=0.05, phi = seq(0,0.15,0.01),
  Plot_Comp = F, MA_PROC=F, leg = 0.6
){
  ind1 = c(
    1,rep(2,multLag),rep(3:4,multLag),
    rep(5,multLag),rep(6,multLag),rep(7,multLag)
  );
  indPt1 = c(
    1,1:multLag,rep(1:multLag,each=2),rep(1:multLag,3)
  )
  indPt2 = c(
    1,rep(1:multLag,each=2),rep(1:multLag,3)
  )
  ind3 = c(
    1,rep(3:4,multLag),
    rep(5,multLag),rep(6,multLag),rep(7,multLag)
  );
  #ind2 = c(1,1:multLag,rep(1:5,each=multLag));
  ind2 = c(1,rep(3,6*multLag));
  pv.m <- (out <= alf)
  pv.m = apply(pv.m,c(1,2),mean);
  len = length(phi)
  if(MA_PROC){
    lab.x = "MA Parameter"
  } else {
    lab.x = "AR Parameter"
  }
  if(Plot_Comp){
    matplot(
      phi, pv.m[,c(1,(2+multLag):(1+6*multLag))], col=ind3,lty=ind3,
      type='b', lwd=3,las=1,ylim=c(0,1),
      ylab="Power",xlab=lab.x,pch=indPt2
    )
    legend(
      "topleft",legend=c(
        "Durbin-Watson",
        "Rotation Conc","Rotation Adj",
        "Permutation", "Reflection", "Rotational Comp"
      ),
      lty=c(1,3:7),lwd=3,col=c(1,3:7)
    )
    legend(
      0,leg,legend=paste("lag",1:multLag),pch=1:multLag,col=1
    )
  } else {
    if(MA_PROC){
      lbIndx = (2+3*multLag):(1+4*multLag) 
    } else {
      lbIndx = (2+6*multLag):(1+7*multLag) 
    }
    matplot(
      phi, pv.m[1:len,c( 1:(1+3*multLag),lbIndx)], 
      col=ind1,lty=ind1,pch=indPt1,
      type='b', lwd=3,las=1,ylim=c(0,1),
      ylab="Power",xlab=lab.x
    )
    matplot(
      phi, pv.m[1:len,1+multLag+c( 2,4,6 )], 
      col=ind1[1+multLag+c( 2,4,6 )],
      lty=ind1[1+multLag+c( 2,4,6 )],
      pch=indPt1[1+multLag+c( 2,4,6 )],
      type='b', lwd=3, add=T
    )
    if(MA_PROC){
      sumLogPv <- apply( 
        -2*log(out[,1+multLag+c( 2,4,6 ),]), c(1,3), sum
      )
      fshpv <- pchisq(
        q = sumLogPv,
        df= 2*multLag,lower.tail=F
      )
      pv.f <- (fshpv <= alf)
      pv.f = apply(pv.f,1,mean);
      lines(  
        phi,pv.f[1:len],lwd=3,lty=6,col=6
      )
      legend(
        "topleft",legend=c(
          "Durbin-Watson","Breusch-Godfrey",
          "Rotation Conc","Rotation Adj",
          "Ljung-Box", "Rotation Fisher"
        ),
        lty=c(1:6),lwd=3,col=c(1:6)
      )
    } else {
      legend(
        "topleft",legend=c(
          "Durbin-Watson","Breusch-Godfrey",
          "Rotation Conc","Rotation Adj",
          "Ljung-Box"
        ),
        lty=c(1:4,5),lwd=3,col=c(1:4,5)
      )
    }
    legend(
      0,leg,legend=paste("lag",1:multLag),pch=1:multLag,col=1
    )
  }
  abline(h=alf,col='gray',lwd=4)
}

plotCompPower <- function( 
    out, multLag=4, alf=0.05, phi = seq(0,0.45,0.09)
){
  ind1 = c(
    1,rep(2,multLag),rep(3:4,multLag),
    rep(5,multLag),rep(6,multLag),rep(7,multLag)
  );
  #ind2 = c(1,1:multLag,rep(1:5,each=multLag));
  ind2 = c(1,rep(3,6*multLag));
  pv.m <- (out <= alf)
  pv.m = apply(pv.m,c(1,2),mean);
  pv.m <- pv.m[,which( (ind1!=2)&(ind1!=3))]
  matplot(
    phi, pv.m, col=ind1[which( (ind1!=2)&(ind1!=3))],
    lty=ind1[which( (ind1!=2)&(ind1!=3))],
    type='l', lwd=3,las=1,ylim=c(0,1),
    ylab="Power",xlab="AR Parameter"
  )
  legend(
    "topleft",legend=c(
      "Durbin-Watson", #"Breusch-Godfrey",
      #"Rotation Conc",
      "Rotation Adj",
      "Permutation", "Reflection", "Rotation Comp"
    ),
    lty=c(1,4:7),lwd=3,col=c(1,4:7)
  )
  abline(h=alf,col='gray',lwd=4)
  return(pv.m)
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
  tot = 1+multLag+2*multLag+multLag
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
      for(ii in 1:multLag)
        pvs[ ll,1+3*multLag+ii,rr ] = doLBTest( dat, lag=lag+ii-1 )
    }
    setTxtProgressBar( pb, rr );
  }
  close(pb);
  
  return(pvs)
}

plotARPowerMult <- function( 
  out, multLag=5, alf=0.05, lags=c(2,4), cf.max = 0.15, inc=0.01, leg = 0.6
){
  reps = dim(out)[3]
  phi = seq(0,cf.max,inc);
  len = length(phi);
  ind1 = c(1,rep(2,multLag),rep(3:4,multLag),rep(5,multLag),6);
  indPt1 = c(
    1,1:multLag,rep(1:multLag,each=2),rep(1:multLag,1),-1
  )
  #ind2 = c(1,1:multLag,rep(1:5,each=multLag));
  ind2 = c(1,rep(3,3*multLag));
  fshpv <- pchisq(
    q = array(
      data = -2*(
        log( out[,3+multLag+lags[1],] ) + 
        log( out[,3+multLag+lags[1]+lags[2],] )
      ),
      dim  = c(len,reps)
    ), df= 4,lower.tail=F
  )
  out  <- abind( out, fshpv, along=2 )
  pv.m <- (out <= alf)
  pv.m = apply(pv.m,c(1,2),mean);
  matplot(
    phi, pv.m[,-(2+4*multLag)], col=ind1,lty=ind1,pch=indPt1,
    type='b', lwd=3,las=1,ylim=c(0,1),
    ylab="Power",xlab="AR Parameter"
  )
  lines(
    phi, pv.m[,2+4*multLag],col=6,lty=6,lwd=3
  )
  legend(
    "topleft",legend=c(
      "Durbin-Watson","Breusch-Godfrey",
      "Rotation Conc","Rotation Adj",
      "Ljung-Box",
      "Rotation Fisher"
    ),
    lty=1:6,lwd=3,col=1:6
  )
  legend(
    0,leg,legend=paste("lag",1:multLag),pch=1:multLag,col=1
  )
  abline(h=alf,col='gray',lwd=4)
}

plotARPowerMult2 <- function( 
  out, multLag=5, alf=0.05, cf.max = 0.15, inc=0.01
){
  reps = dim(out)[3]
  phi = seq(0,cf.max,inc);
  len = length(phi);
  ind1 = c(1,rep(2,multLag),rep(3:4,multLag),rep(5,multLag),rep(6,multLag));
  #ind2 = c(1,1:multLag,rep(1:5,each=multLag));
  ind2 = c(1,rep(3,3*multLag));
  
  out  <- abind( out, -2*log(out[,3+multLag,]), along=2 );
  indx <- 1+multLag + 2*(1:multLag); #c(8,10,12,14,16);
  for( ii in 2:multLag ){
    out  <- abind( 
      out, out[,2+4*multLag+ii-2,] - 2*log(out[,indx[ii],]), 
      along=2 
    )
  }
  for( ii in 1:multLag ){
    fshpv <- 
      pchisq(
      q = array(
        data = out[,2+4*multLag+ii-1,],
        dim  = c(len,reps)
      ), df= 2*ii,lower.tail=F
    )
    out[,2+4*multLag+ii-1,]  <- fshpv;
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
      "Ljung-Box",
      "Rotation Fisher"
    ),
    lty=1:6,lwd=3,col=1:6
  )
  abline(h=alf,col='gray',lwd=4)
}

runARPosCorTests <- function(
    nn = 1000, cf.max = 0.15, inc=0.01, 
    lag= 1,reps=100,offset=1,multLag=5,
    dof=Inf,ord=3,BB=100
){
  phi = seq(0,cf.max,inc);
  len = length(phi);
  tot = 1+multLag+2*multLag+1*multLag;
  pvs = array( 0, dim=c(len,tot,reps) );
  
  pb = txtProgressBar( min=0,max=reps,style=3 )
  for( rr in 1:reps ){
    for( ll in 1:len ){
      # Generate Data
      if(ll==1){
        dat = simMAproc( nn, list(), dof=dof )
      } else {
        par = c( rep(0,offset), rep(phi[ll],ord) )
        dat = simMAproc( nn, list(ma=par), dof=dof )
      }
      # Do Tests
      pvs[ ll,1,rr ]   = doDWTest( dat )
      for(ii in 1:multLag)
        pvs[ ll,1+ii,rr ]   = doBGTest( dat, lag=lag+ii-1 )
      for(ii in 1:multLag)
        pvs[ ll,1+multLag+2*(ii-1)+1:2,rr ] = doSOTest( dat, lag=lag+ii-1 )
      for(ii in 1:multLag)
        pvs[ ll,1+3*multLag+ii,rr ] = doLBTest( dat, lag=lag+ii-1 )
      if(F){
       for(ii in 1:multLag)
         pvs[ ll,1+3*multLag+ii,rr ] = doRandTest( dat, lag=lag+ii-1, BB,"PERM" )
       for(ii in 1:multLag)
         pvs[ ll,1+4*multLag+ii,rr ] = doRandTest( dat, lag=lag+ii-1, BB,"SIGN" )
       for(ii in 1:multLag)
         pvs[ ll,1+5*multLag+ii,rr ] = doRandTest( dat, lag=lag+ii-1, BB,"BOTH" )
       }
    }
    setTxtProgressBar( pb, rr );
  }
  close(pb);
  
  return(pvs)
}

runHeteroTests <- function(
    nn = 1000, cf.max = 4, inc=0.5, 
    lag= 1,reps=100,offset=1,multLag=5,
    dof=Inf,ord=3,BB=100
){
  phi = seq(0,cf.max,inc);
  len = length(phi);
  tot = 1+multLag+2*multLag+3*multLag;
  pvs = array( 0, dim=c(len,tot,reps) );
  
  pb = txtProgressBar( min=0,max=reps,style=3 )
  for( rr in 1:reps ){
    for( ll in 1:len ){
      # Generate Data
      #if(ll==1){
        dat = simMAproc( nn, list(), dof=dof )
      #} else {
      #  par = c( rep(0,offset), rep(phi[ll],ord) )
      #  dat = simMAproc( nn, list(), dof=dof )
      #}
      dat <- dat*seq(1,1+phi[ll],length.out=nn)
      # Do Tests
      pvs[ ll,1,rr ]   = doDWTest( dat )
      for(ii in 1:multLag)
        pvs[ ll,1+ii,rr ]   = doBGTest( dat, lag=lag+ii-1 )
      for(ii in 1:multLag)
        pvs[ ll,1+multLag+2*(ii-1)+1:2,rr ] = doSOTest( dat, lag=lag+ii-1 )
      for(ii in 1:multLag)
        pvs[ ll,1+3*multLag+ii,rr ] = doRandTest( dat, lag=lag+ii-1, BB,"PERM" )
      for(ii in 1:multLag)
        pvs[ ll,1+4*multLag+ii,rr ] = doRandTest( dat, lag=lag+ii-1, BB,"SIGN" )
      for(ii in 1:multLag)
        pvs[ ll,1+5*multLag+ii,rr ] = doRandTest( dat, lag=lag+ii-1, BB,"BOTH" )
    }
    setTxtProgressBar( pb, rr );
  }
  close(pb);
  
  return(pvs)
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
  ARRES=F, ylm = NULL, arOrd=1
){
  tmp= log2(dat+1)
  ks = ksmooth(rng,tmp[rng],bandwidth=bwd);
  pvs= matrix( 0, nrow=lag.mx, ncol=2 );
  df = (tmp[rng]-ks$y)
  if(ARRES){
    md <- arima(df,order=c(arOrd,0,0));
    df <- resid( md )
  }
  df  <- df[-1]
  rng <- rng[-1]
  for( ii in 1:lag.mx ){
    pvs[ii,1] = doBGTest( df, ii );
    pvs[ii,2] = doSOTest( df, ii )[2];
  }
  if(is.null(ylm)){
    matplot( 
      pvs, col=c('red','blue'),lty=2:3,log='y',
      type='b',las=1,lwd=3,pch=1:2, 
      ylab='p-values',xlab='lag'
    )
  } else {
    matplot( 
      pvs, col=c('red','blue'),lty=2:3,log='y',
      type='b',las=1,lwd=3,pch=1:2, 
      ylab='p-values',xlab='lag',ylim=ylm
    )
  }
  abline(h=c(0.05,0.01,0.001),lwd=3,col=c('lightgray','gray','darkgray'))
  text( 
    x=rep(18,3),y=c(0.05,0.01,0.001),
    pos=3,labels=c("5%","1%","0.1%"),offset=0.2
  )
  legend(
    15,0.0001,legend=c( "Breusch-Godfrey","Rotation Adj" ), 
    lwd=3, lty=2:3, pch=1:2, col=c("red","blue")
  )
  pvRot  <- pvs[,2];
  pvSort <- sort(pvRot);
  pvOrd  <- order(pvRot);
  pvAdj  <- p.adjust(pvSort, method="BH")
  print( cbind(pvOrd[1:10],pvSort[1:10],pvAdj[1:10]) )
  return(pvs)
}


############
# Timing Test
############

runTimingTest <- function( 
    nn = c(1e2,1e3,1e4,1e5), pp = 1000  
){
  len = length(nn)
  rtme= matrix(0,len,3)
  for( ii in 1:len ){
    print( paste("test: n =",nn[ii])  )
    rtme[ii,1] = system.time(
      for(jj in 1:nn[ii])
        tmp = 2*rbinom(pp,0,1)-1
    )[3]
    rtme[ii,2] = system.time(
      for(jj in 1:nn[ii])
        tmp = sample(pp)
    )[3]
    rtme[ii,3] = system.time(
      for(jj in 1:nn[ii]){
        mat = matrix( rnorm(pp^2),pp,pp )
        tmp = qr(mat)
      }
    )[3]
  }
  return(rtme)
}

if(F){
  matplot( 
    c(1e2,1e3,1e4,1e5),t.out4,type='b',log='xy',
    las=1,ylab='',xlab="Replications",
    main="Runtimes by Dimension", pch=1:3,lwd=2,lty=1,
    col='blue',ylim = c(1e-2,1e4)
  )
  matplot( 
    c(1e2,1e3,1e4,1e5),t.out2,type='b',
    pch=1:3,lwd=2,lty=2,
    col='red',add = T
  )
  matplot( 
    c(1e2,1e3,1e4,1e5),t.out,type='b',
    pch=1:3,lwd=2,lty=3,
    col='green',add = T
  )
  legend(
    x=1e2,y=1e4,legend = c( "Reflect","Permute","Rotate" ),
    pch = 1:3,bty = 'n'
  )
  legend(
    x=5e2,y=1e4,legend = c( "n = 100","n = 200","n = 400" ),
    lty = 3:1,lwd=2, col=c("green","red","blue"),bty = 'n'
  )
}


##
constCheck <- function(nn=100,hh=1){
  val1 = 32*nn*(nn+2)/( (nn-hh)*(nn-2) )
  val2 = 32*nn*(nn+2)/( nn-hh )
  return(
    c(
      0.5*log(val1) + lgamma(val1) - lgamma(0.5+val1),
      0.5*log(val2) + lgamma(val2) - lgamma(0.5+val2)
    )
  )
}
##
