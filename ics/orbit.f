C Subroutines for use with 'maketrail'
C    * Integrates orbit from current time (t=tnow)  
C     backwards to t=0 and forwards to t=2*tnow.
C    * Fills arrays at points equally spaced in angle 
C     psi along the trail.
C    * Notes times and distances of pericenters.
C    * Subroutine 'accel' can be replaced with any Galactic model
C
C Input from maketrail:
C     initial phase-space coords in x(6)
C     tnow       
C
C Outputs in 'orb.dat' and to maketrail:
C At points equally spaced in angle along orbit
C     tout       time in Gyears
C     y          phase-space coords (kpc,km/s)
C     psi        angle along orbit in radians
C     pdot       angular speed in radians/Gyr
C     e          energy (to check conservation)
C
C Requires in common blocks:
C 'reals':
C     x(6)=(x,v) phase-space coords
C     t          time in simulations units
C     dt         time step for integration in simulation units
C     tu         unit time in Gyears
C 'coords':
C     dxdt(6)=(v,a) 
C     pot        total potential
C
C***********************************************************************
C 
C
        SUBROUTINE accel
C
C 
C***********************************************************************
C Finds acceleration in Spergel's model of the Galaxy
C ---------------------------------------------------
C Unit velocity is  km/s and length is kpc.
C
      INCLUDE 'maketrail.h'
        REAL*8 x,t,dt,dxdt,pot,Z,
     &         a,b,GM,GMs,rs,
     &         tsrad,sqz2b2,tdr,tdz,phim,phis,
     &         r,r2,phih,thrad,
     &         vcirc2,tu,
     &         qx,qy,qz,x2law,y2law,z2law,thtemp,
     &         thxlaw,thylaw,thzlaw,c,vh2,dfx,dfy,dfz,rothalo,
     &         rc1,rc2,rc3,
     &         tsradSUN,tdrSUN,dxdtSUN,tempSUN
        LOGICAL firstc
        COMMON/reals/x(6),t,dt,tu,c,vh2,dfx,dfy,dfz
        COMMON/coords/dxdt(6),pot

        DATA firstc/.TRUE./
        SAVE firstc,rs,GMs,b,a,GM

        IF(firstc)THEN
           firstc=.FALSE.
           a=6.5
           b=.26
           rs=.7
C           c=13. is best fit for q=0.9, 11 for q=1.25
           c=dhalo
C           vcirc=220 matches masses given in Law 2005
           vcirc2=220.**2

C vc_halo=225, approx=200
C           vh2 = 0.4*vcirc2
C vc_halo=200, approx=200
C           vh2 = 0.325*vcirc2
C           vh2 = 0.325*vcirc2
C vc_halo=175, approx=200km/s
C           vh2=0.2*vcirc2
C           GM = 8.887*vcirc2
C           GMs = 3.0*vcirc2
C Modifications for testing flattening 6/20/04
           GM = 1.0*8.887*vcirc2
           GMs = 3.0*vcirc2
C          129 for oblate, 115 for prolate
           vh2=vhalo*vhalo
C           vh2 = 163.**2
C     vh2 now set by code itself so 220 km/s LSR speed
        ENDIF

C Find dynamical friction vector
C        CALL dynfric
C                                       -------

C SPHEROID
        r2=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
        r=SQRT(r2)
        tsrad = GMs/(r+rs)**2/r
        phis = -GMs/(r+rs)
C DISK
        R=SQRT(x(1)*x(1)+x(2)*x(2))
        Z=ABS(x(3))
        sqz2b2 = SQRT(Z*Z + b*b)
        tdr = GM/(R*R + (a + sqz2b2)**2)**1.5
        tdz = tdr*(a/sqz2b2 + 1.)
        phim = -GM/sqrt(R*R+(a+sqz2b2)**2)

C Accel from spheroid and disk
        dxdt(4) = -(tsrad+tdr)*x(1)
        dxdt(5) = -(tsrad + tdr)*x(2)
        dxdt(6) = -(tsrad+tdz)*x(3)

C HALO
C David editing code: allow for flattened or triaxial halo
C q_law is flattening, p_law is triaxiality
C Keep qy=1.00, qx is my variable to go with phi
        qx=q1
        qy=q2
        qz=q3
C Rotation angle of halo- set in param2 for ease of scripting
        rothalo=drlphi*3.14159/180.
C Rotation constants
        rc1=cos(rothalo)*cos(rothalo)/qx/qx 
     & + sin(rothalo)*sin(rothalo)/qy/qy
        rc2=cos(rothalo)*cos(rothalo)/qy/qy
     & + sin(rothalo)*sin(rothalo)/qx/qx
        rc3=2*sin(rothalo)*cos(rothalo)*(1/qx/qx -1/qy/qy)

        x2law=x(1)*x(1)
        y2law=x(2)*x(2)
        z2law=x(3)*x(3)

C Original fragment
C        thtemp=2.*vh2/(x2law/qx/qx+y2law/qy/qy+z2law/qz/qz+c*c)
C        thxlaw=thtemp/qx/qx
C        thylaw=thtemp/qy/qy
C        thzlaw=thtemp/qz/qz
C        phih=vh2*LOG(x2law/qx/qx+y2law/qy/qy+z2law/qz/qz+c*c)
C Modified fragment rotating in XY plane
        thtemp=vh2/(x2law*rc1+y2law*rc2+x(1)*x(2)*rc3+z2law/qz/qz+c*c)
        thxlaw=thtemp*(2*rc1*x(1)+rc3*x(2))
        thylaw=thtemp*(2*rc2*x(2)+rc3*x(1))
        thzlaw=thtemp*(2*x(3)/qz/qz)
        phih=vh2*LOG(x2law*rc1+y2law*rc2+rc3*x(1)*x(2)+z2law/qz/qz+c*c)


        thrad = 2.*vh2/(r2+c*c)

C Add in accel from halo
        dxdt(4) = dxdt(4)-thxlaw
        dxdt(5) = dxdt(5)-thylaw
        dxdt(6) = dxdt(6)-thzlaw

        pot=phih+phim+phis
   
        RETURN
        END


C***********************************************************************
C
C
      SUBROUTINE integrate(nsteps)
C
C
C***********************************************************************
C Integrate using simple leapfrog alogorithm
        INTEGER i,nsteps
        REAL*8 dt,x,t,tu,c,vh2,dfx,dfy,dfz
        COMMON/reals/x(6),t,dt,tu,c,vh2,dfx,dfy,dfz

        CALL accel
C            -----
        CALL stepvel(dt/2.)
C            -------
        DO 10 i=1,nsteps
         
           CALL steppos
C               -------
           CALL accel
C               -----
           CALL stepvel(dt/2.)
C               ------
           CALL out
C               ---
           IF(i.LT.nsteps)CALL stepvel(dt/2.)
C               ------
 10     CONTINUE

        RETURN
        END

C***********************************************************************
C
C
      SUBROUTINE orbit
C
C
C***********************************************************************

      INCLUDE 'maketrail.h'
      INTEGER nsteps,i,j
      REAL*8 r,v
      REAL*8 x,t,dt,dxdt,pot,kmpkpc,secpgyr,tu,c,vh2,dfx,dfy,dfz
      PARAMETER(kmpkpc=3.085678e16,secpgyr=60.*60.*24.*365.*1.e9)
      COMMON/reals/x(6),t,dt,tu,c,vh2,dfx,dfy,dfz
      COMMON/coords/dxdt(6),pot

C tu=unit time in Gyrs
      tu=kmpkpc/secpgyr

C Set timestep to 1/100th of initial crossing time
      r=SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      v=SQRT(x(4)*x(4)+x(5)*x(5)+x(6)*x(6))
      dt=-r/v/100.
      nsteps=-INT(tnow/tu/dt)+1
      nout=0
      nenc=0
C
      t=tnow/tu
C Integrate backwards
      CALL integrate(nsteps)
C          ---------
C Reorder arrays to start from t=0
      CALL reorder
C          -------
C Integrate forwards
      t=tnow/tu
      dt=-dt
      DO 10 i=1,6
         x(i)=y(i,nout)
 10   CONTINUE
      CALL integrate(nsteps)
C          ---------

      OPEN(UNIT=10,FILE='orb.dat',STATUS='UNKNOWN')
      DO 15 i=1,nout
         IF(tout(i).GT.tenc(1))
     &    WRITE(10,99)tout(i)-tnow,(y(j,i),j=1,6),psi(i),pdot(i),e(i)
         rr(i)=SQRT(y(1,i)**2+y(2,i)**2+y(3,i)**2)
         vrr(i)=(y(1,i)*y(4,i)+y(2,i)*y(5,i)+y(3,i)*y(6,i))/rr(i)
 15   CONTINUE
 99   FORMAT(10(1pe12.4))

      RETURN
      END

C***********************************************************************
C
C
        SUBROUTINE out
C
C
C***********************************************************************
C saves at angles spaced by pi/save
C and finds peri time and distance
C -------------------------------------------
C
        INCLUDE 'maketrail.h'
        INTEGER i
        REAL*8 r,dpsi,dp,rp(5),r1,r2,rdotr,vr,save
        PARAMETER(save=200.)
        REAL*8 x,t,dt,dxdt,pot,tu,c,vh2,dfx,dfy,dfz
        COMMON/reals/x(6),t,dt,tu,c,vh2,dfx,dfy,dfz
        COMMON/coords/dxdt(6),pot

        SAVE rp

        dp=4.*ATAN(1.0)/save

C Find change in angle since last output
        IF(nout.GT.0)THEN
           rdotr=x(1)*y(1,nout)+x(2)*y(2,nout)+x(3)*y(3,nout)
           r1=SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
           r2=SQRT(y(1,nout)*y(1,nout)+y(2,nout)*y(2,nout)
     &          +y(3,nout)*y(3,nout))
           dpsi=ACOS(rdotr/r1/r2)
        END IF

        r=SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
        IF(dpsi.GT.dp.OR.nout.EQ.0)THEN
           nout=nout+1
           psi(nout)=0.
           IF(nout.GT.1)psi(nout)=psi(nout-1)+dpsi
           y(1,nout)=x(1)
           y(2,nout)=x(2)
           y(3,nout)=x(3)
           y(4,nout)=x(4)
           y(5,nout)=x(5)
           y(6,nout)=x(6)
           tout(nout)=t*tu
C
           e(nout)=(x(4)*x(4)+x(5)*x(5)+x(6)*x(6))/2.+pot
C
           r=SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
           vr=(x(1)*x(4)+x(2)*x(5)+x(3)*x(6))/r
           pdot(nout)=SQRT(x(4)*x(4)+x(5)*x(5)+x(6)*x(6)-vr*vr)/r/tu
C
        END IF

C Peri
        IF(dt.LT.0.)THEN
           DO 10 i=1,4
              rp(i)=rp(i+1)
 10        CONTINUE
           rp(5)=r
           IF(rp(4).GT.rp(3).AND.rp(2).GT.rp(3))THEN
              IF(rp(5).GT.rp(4).AND.rp(1).GT.rp(2))THEN
                 nenc=nenc+1
                 peri(nenc)=rp(3)
                 tenc(nenc)=(t-2.*dt)*tu
              ENDIF
           END IF
        END IF

        RETURN
        END

C***********************************************************************
C
C
      SUBROUTINE reorder
C
C
C***********************************************************************
C Reorders arrays to start from t=0
      INCLUDE 'maketrail.h'
      INTEGER i,j,k
      REAL*8 dum1(nmax),dum2(nmax),dum3(nmax),dum4(nmax),dum(6,nmax)
      REAL*8 x,t,dt,tu,c,vh2,dfx,dfy,dfz
      COMMON/reals/x(6),t,dt,tu,c,vh2,dfx,dfy,dfz

      DO 10 i=1,nout
         dum1(i)=tout(i)
         dum2(i)=psi(i)
         dum3(i)=pdot(i)
         dum4(i)=e(i)
         DO 15 j=1,6
            dum(j,i)=y(j,i)
 15      CONTINUE
 10   CONTINUE

      DO 20 i=1,nout
         k=nout-i+1
         tout(i)=dum1(k)
         psi(i)=-dum2(k)
         pdot(i)=dum3(k)
         e(i)=dum4(k)
         DO 25 j=1,6
            y(j,i)=dum(j,k)
 25      CONTINUE
 20   CONTINUE

      DO 30 i=1,nenc
         dum1(i)=tenc(i)
         dum2(i)=peri(i)
 30   CONTINUE
      DO 35 i=1,nenc
         k=nenc-i+1
         tenc(i)=dum1(k)
         peri(i)=dum2(k)
 35   CONTINUE

      RETURN
      END
C***********************************************************************
C
C
        SUBROUTINE steppos
C
C
C***********************************************************************
C
        INTEGER i
        REAL*8 x,t,dt,dxdt,pot,tu,c,vh2,dfx,dfy,dfz
        COMMON/reals/x(6),t,dt,tu,c,vh2,dfx,dfy,dfz
        COMMON/coords/dxdt(6),pot

        DO 5 i=1,3
           x(i)=x(i)+dt*dxdt(i)
 5      CONTINUE
        t=t+dt

        RETURN
        END

C***********************************************************************
C
C
        SUBROUTINE stepvel(dt2)
C
C
C***********************************************************************
C Changes velocities for integration and output
C ----------------------------
        INTEGER i
        REAL*8 x,t,dt,dt2,dxdt,pot,tu,c,vh2,dfx,dfy,dfz
        COMMON/reals/x(6),t,dt,tu,c,vh2,dfx,dfy,dfz
        COMMON/coords/dxdt(6),pot

        DO 10 i=4,6
                x(i)=x(i)+dxdt(i)*dt2
                dxdt(i-3)=x(i)

 10     CONTINUE

        RETURN
        END

C***********************************************************************
C
C
      SUBROUTINE dynfric
C
C       
C***********************************************************************
C Subroutine to calculate dynamical friction due to halo component alone
      INCLUDE 'maketrail.h'
      REAL*8 x,t,dt,tu,r,mr,v,
     &       coulog,rho,sigma2,xx,df,erff,mfricnow,tfrac,tstart,thalf,
     &       ot,xlim,rtide,jfit,G,pi,c,vh2,dfx,dfy,dfz,massloss,ttemp,
     &       ttemp2
      COMMON/reals/x(6),t,dt,tu,c,vh2,dfx,dfy,dfz
      LOGICAL firstc
      DATA firstc/.TRUE./
      SAVE ot,sigma2,firstc,G,pi,tstart,massloss

      IF(firstc)THEN
         firstc=.FALSE.
         ot=1.0d0/3.0d0
         sigma2=vh2
C G in units of kpc*km*km/s/s/Msun
         G=4.3e-6
         pi=4.0*ATAN(1.0)
C tstart is time of apo at which full sim will be started
         tstart=2.75
C thalf is time at which half of original mass left (for Hernquist model)
         thalf=1.0
         ttemp2=-tstart+thalf
C massloss is fraction of original mass lost over interaction
         massloss=0.9
      ENDIF

C Allow a linear mass-loss rate for the dynamical friction calculation
C Only start mass decay once reach the starting point for the full scf run

C 1-part linear decay in frictional mass
      IF((t*tu-tnow).GE.(-tstart))THEN
       tfrac=(t*tu-tnow+tstart)/tstart
       mfricnow=mfric*(1-tfrac*massloss)

C 2-part linear decay in frictional mass
C      ttemp=t*tu-tnow
C      IF(((ttemp).GE.(-tstart)).AND.((ttemp).LE.(ttemp2)))THEN
C        tfrac=(t*tu-tnow+tstart)/thalf
C        mfricnow=mfric*(1-tfrac*0.5)
C      ELSE IF(((ttemp).GE.(-tstart)).AND.((ttemp).GE.(ttemp2)))THEN
C        tfrac=(t*tu-tnow+tstart-thalf)/(tstart-thalf)
C        mfricnow=0.5*mfric*(1-tfrac*0.9)
C      ELSE
C        mfricnow=mfric

      ENDIF

      r=SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))

C Calculate mass enclosed within r (=(1/r^2 G)(d phi /dr))
      mr=2.0d0*vh2*r*r*r/(r*r+c*c)/G

C Calculate speed - note: strictly should be updated so speed coincides with position
      v=SQRT(x(4)*x(4)+x(5)*x(5)+x(6)*x(6))

C Coulomb logarithm via James' fitting function
      rtide=r*(ot*mfric/mr)**ot
      xlim=rtide
      coulog=DLOG(r/rtide)+jfit(xlim)+0.135d0

      IF(coulog.GT.0)THEN
C Local density (=(1/4 pi r^2)(d mr/dr))
         rho=(vh2/2.0d0/pi/G)*(r*r+3.0d0*c*c)/(r*r+c*c)/(r*r+c*c)

         xx=v/SQRT(2.d0*sigma2)
 
C B+T, eqn 7-17
         df=4.0d0*pi*coulog*(G**2)*mfricnow*rho*
     &          (ERFF(xx)-2.0d0*xx*EXP(-(xx*2))/SQRT(pi))/(v**3)
      ELSE
         df=0
      ENDIF
         dfx=-df*x(4)
         dfy=-df*x(5)
         dfz=-df*x(6)
      RETURN
      END

C***********************************************************************
C***********************************************************************
C From Numerical Recipes
C***********************************************************************

      FUNCTION erff(x)
      REAL*8 erff,x
CU    USES gammp
      REAL*8 gammp
      if(x.lt.0.)then
        erff=-gammp(0.5d0,x**2)
      else
        erff=gammp(0.5d0,x**2)
      endif
      return
      END
C***********************************************************************
      FUNCTION gammp(a,x)
      REAL*8 a,gammp,x
CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      END
C***********************************************************************
      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END

C***********************************************************************
      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammln
      INTEGER n
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END

C***********************************************************************
C
C
        FUNCTION GAMMLN(XX)
C
C
C***********************************************************************
C
C
C     A routine to compute the natural logarithm of the gamma
C     function.  (Taken from numerical recipes.)
C
C
C=======================================================================

        INTEGER j

        REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER,gammln,xx

        DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     &      -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
        DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/

        X=XX-ONE
        TMP=X+FPF
        TMP=(X+HALF)*LOG(TMP)-TMP
        SER=ONE

        DO 11 J=1,6
          X=X+ONE
          SER=SER+COF(J)/X
11      CONTINUE

        GAMMLN=TMP+LOG(STP*SER)

        RETURN
        END

C***********************************************************************
C
C
      FUNCTION jfit(x)
C
C       
C***********************************************************************

      REAL*8 jfit, x, num, den, fx
      REAL*8 a,b,c,d,e,f,g,h,i,j

CCCC Fitting parameters.
      
      PARAMETER (a=0.2021d0,b=0.04474d0,c=3.1621d0,
     & d=0.5488d0,e=2.9965d0,f=0.7375d0,g=0.9217d0,
     & h=0.2749d0,i=1.2583d0,j=2.3069d0)

CCCC

      num = DLOG(1.d0+a*x)*(b*(x**c)+d*(x**e))
      den = (1.d0+f*(x**g)+h*(x**i))**j
      fx=log(1.d0+x)-x/(1.d0+x)

      jfit = num/den
      jfit=jfit/fx/fx

      RETURN
      END
