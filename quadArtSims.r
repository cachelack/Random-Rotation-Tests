
source("quadFormTest.")

# Small Sample Size Tests

xx.comp1 = runAR1Tests(
  COMPSIM=F,offset=0,BB=100,multLag=3,n=100,reps=1000,cf.max=0.4,inc=0.1
); 
plotARPower( xx.comp1, multLag=3, alf=0.05, phi=seq(0,0.4,0.1) )

xx.comp3 = runAR1Tests(
  COMPSIM=F,offset=2,BB=100,multLag=3,n=100,reps=1000,cf.max=0.5,inc=0.5/4
)
plotARPower( xx.comp3, multLag=3, alf=0.05, phi=seq(0,0.4,0.1) )

# AR(1) Tests

xx.ar1.df1 <- runAR1Tests(
  COMPSIM=F,offset=0,BB=100,multLag=3,
  n=1000,reps=2000,cf.max=0.15,inc=0.01, dof=1
);

xx.ar1.df2 <- runAR1Tests(
  COMPSIM=F,offset=0,BB=100,multLag=3,
  n=1000,reps=2000,cf.max=0.15,inc=0.01, dof=2
);

xx.ar1.df4 <- runAR1Tests(
  COMPSIM=F,offset=0,BB=100,multLag=3,
  n=1000,reps=2000,cf.max=0.15,inc=0.01, dof=4
);

xx.ar1.dfI <- runAR1Tests(
  COMPSIM=F,offset=0,BB=100,multLag=3,
  n=1000,reps=2000,cf.max=0.15,inc=0.01, dof=Inf
);

pdf("PICS/powAR1dof1_new.pdf",width=6,height=5)
 plotARPower( xx.ar1.df1, multLag=3, alf=0.05, phi=seq(0,0.15,0.01) )
 title("Cauchy Distribution")
dev.off()
pdf("PICS/powAR1dof2_new.pdf",width=6,height=5)
 plotARPower( xx.ar1.df2, multLag=3, alf=0.05, phi=seq(0,0.15,0.01) )
 title("t(2) Distribution")
dev.off()
pdf("PICS/powAR1dof4_new.pdf",width=6,height=5)
 plotARPower( xx.ar1.df4, multLag=3, alf=0.05, phi=seq(0,0.15,0.01) )
 title("t(4) Distribution")
dev.off()
pdf("PICS/powAR1dofI_new.pdf",width=6,height=5)
 plotARPower( xx.ar1.dfI, multLag=3, alf=0.05, phi=seq(0,0.15,0.01) )
 title("Gaussian Distribution")
dev.off()

# AR(3) Tests

xx.ar3.df1 <- runAR1Tests(
  COMPSIM=F,offset=2,BB=100,multLag=3,
  n=1000,reps=2000,cf.max=0.15,inc=0.01, dof=1
);

xx.ar3.df2 <- runAR1Tests(
  COMPSIM=F,offset=2,BB=100,multLag=3,
  n=1000,reps=2000,cf.max=0.15,inc=0.01, dof=2
);

xx.ar3.df4 <- runAR1Tests(
  COMPSIM=F,offset=2,BB=100,multLag=3,
  n=1000,reps=2000,cf.max=0.15,inc=0.01, dof=4
);

xx.ar3.dfI <- runAR1Tests(
  COMPSIM=F,offset=2,BB=100,multLag=3,
  n=1000,reps=2000,cf.max=0.15,inc=0.01, dof=Inf
);

pdf("PICS/powAR3dof1_new.pdf",width=6,height=5)
 plotARPower( xx.ar3.df1, multLag=3, alf=0.05, phi=seq(0,0.15,0.01) )
 title("Cauchy Distribution")
dev.off()
pdf("PICS/powAR3dof2_new.pdf",width=6,height=5)
 plotARPower( xx.ar3.df2, multLag=3, alf=0.05, phi=seq(0,0.15,0.01) )
 title("t(2) Distribution")
dev.off()
pdf("PICS/powAR3dof4_new.pdf",width=6,height=5)
 plotARPower( xx.ar3.df4, multLag=3, alf=0.05, phi=seq(0,0.15,0.01) )
 title("t(4) Distribution")
dev.off()
pdf("PICS/powAR3dofI_new.pdf",width=6,height=5)
 plotARPower( xx.ar3.dfI, multLag=3, alf=0.05, phi=seq(0,0.15,0.01) )
 title("Gaussian Distribution")
dev.off()

# MultLag Test

xx.mult.24.dofI <- runARMultLagTests(
  nn=1000,reps=2000,cf.max=0.15,inc=0.01, dof=Inf,
  multLag=4,offset=1
);

xx.mult.24.dof2 <- runARMultLagTests(
  nn=1000,reps=2000,cf.max=0.15,inc=0.01, dof=2,
  multLag=4,offset=1
);

pdf("PICS/powARMultdofI_new.pdf",width=6,height=5)
 plotARPowerMult( xx.mult.24.dofI, multLag=4, alf=0.05, leg=0.5 )
 title("Gaussian Distribution")
dev.off()
pdf("PICS/powARMultdof2_new.pdf",width=6,height=5)
 plotARPowerMult( xx.mult.24.dof2, multLag=4, alf=0.05, leg=0.5 )
 title("t(2) Distribution")
dev.off()

# MA Process

xx.ma1.dfI <- runARPosCorTests(
  offset=0,multLag=3,
  n=1000,reps=2000,cf.max=0.15,inc=0.01, dof=Inf
)

xx.ma2.dfI <- runARPosCorTests(
  offset=1,multLag=3,
  n=1000,reps=2000,cf.max=0.15,inc=0.01, dof=Inf
)

pdf("PICS/powMA1dofI_new.pdf",width=6,height=5)
 plotARPower( 
   xx.ma1.dfI, multLag=3, alf=0.05, phi=seq(0,0.12,0.01),
   MA_PROC=T,leg=0.5
 )
 title("MA(1) Process")
dev.off()
pdf("PICS/powMA2dofI_new.pdf",width=6,height=5)
 plotARPower( 
   xx.ma2.dfI, multLag=3, alf=0.05, phi=seq(0,0.14,0.01),
   MA_PROC=T,leg=0.5
 )
 title("MA(2) Process")
dev.off()


# Computational Randomization Tests

xx.comp1 = runAR1Tests(
  COMPSIM=T,offset=0,BB=100,multLag=3,n=100,reps=1000,cf.max=0.4,inc=0.1
);
xx.comp3 = runAR1Tests(
  COMPSIM=T,offset=2,BB=100,multLag=3,n=100,reps=1000,cf.max=0.4,inc=0.1
);

pdf("PICS/powComp1dofI_new.pdf",width=6,height=5)
 plotARPower( xx.comp1, multLag=3, alf=0.05, phi=seq(0,0.4,0.1),Plot_Comp=T )
 title("AR(1) Process")
dev.off()
pdf("PICS/powComp3dofI_new.pdf",width=6,height=5)
 plotARPower( xx.comp3, multLag=3, alf=0.05, phi=seq(0,0.4,0.1),Plot_Comp=T )
 title("AR(3) Process")
dev.off()

# Solar Data

sol.bwd3  = getSolarResid( xx, bwd=3,  lag.mx=24, rng = 1:3314, ARRES=F )
sol.bwd6  = getSolarResid( xx, bwd=6,  lag.mx=24, rng = 1:3314, ARRES=F )
sol.bwd12 = getSolarResid( xx, bwd=12, lag.mx=24, rng = 1:3314, ARRES=F )
sol.bwd24 = getSolarResid( xx, bwd=24, lag.mx=24, rng = 1:3314, ARRES=F )
sol.bwd60 = getSolarResid( xx, bwd=60, lag.mx=24, rng = 1:3314, ARRES=F )


sol.bwd3.x  = getSolarResid( 
  xx, bwd=3,  lag.mx=160, rng = 1813:3313, ARRES=F, ylm=c(10^-20,1)
)
sol.bwd6.x  = getSolarResid( 
  xx, bwd=6,  lag.mx=160, rng = 1813:3313, ARRES=F, ylm=c(10^-20,1)
)
sol.bwd12.x = getSolarResid( 
  xx, bwd=12, lag.mx=160, rng = 1813:3313, ARRES=F, ylm=c(10^-20,1)
)
sol.bwd18.x = getSolarResid( 
  xx, bwd=18, lag.mx=160, rng = 1813:3313, ARRES=F, ylm=c(10^-20,1)
)
sol.bwd24.x = getSolarResid( 
  xx, bwd=24, lag.mx=160, rng = 1813:3313, ARRES=F, ylm=c(10^-20,1)
)
sol.bwd60.x = getSolarResid( 
  xx, bwd=60, lag.mx=160, rng = 1813:3313, ARRES=F, ylm=c(10^-20,1)
)
