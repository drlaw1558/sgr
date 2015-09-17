C***************************************************************
C
C
c	PROGRAM scfltr
C
C This variant codes leading/trailing, and permits the halo
C to be rotated in the disk plane
C
C
C***************************************************************
C
C Heads the scf program to run parameter space explorations 
C
C ---------------------------------------------------------------
	INCLUDE 'scfltr.h'

	INTEGER i,nruns
	REAL*8 u,v,w,rr

C Read in parameters.
	OPEN(UNIT=9,FILE='RUNPAR',STATUS='OLD')
	READ(9,*)nruns
	DO 10 i=1,nruns
           nproc=i
	   READ(9,*)
	   READ(9,*)tend
	   READ(9,*)mu,ru
	   READ(9,*)drlphi
	   READ(9,*)q1,q2,q3
	   READ(9,*)dhalo
	   READ(9,*)vhalo
C	   READ(9,*)mfric
	   READ(9,*)xframe,yframe,zframe,vxframe,vyframe,vzframe
C Stopping point for large msat
	   READ(9,*)delx,delv
	   READ(9,*)xend,yend,zend,u,v,w
	   rr=SQRT(xend*xend+yend*yend+zend*zend)
	   vend=(xend*u+yend*v+zend*w)/rr
C	   mfric=mfric/mu
C Run scf.
	   CALL scfm
C               ----
 10     CONTINUE
	CLOSE(9)

	STOP
	END


C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

C**************************************************************
C
C
                             SUBROUTINE scfm
C
C
C***********************************************************************
C
C
C    A code to evolve self-gravitating systems using a self-consistent
C    field approach.  This version has been optimized for supercomputers
C    and is fully vectorized.  The code is written in standard FORTRAN,
C    although some CRAY-specific vector intrinsic routines have been
C    used.
C
C    The computational system of units is determined by the input data.
C    No explicit assumtions about the value of the gravitational 
C    constant have been made; it is read in as a parameter.
C    Particles are not required to have identical masses.
C
C
C                       Version 1: January 1, 1991
C
C
C                    Lars Hernquist, U.C. Santa Cruz
C
C
C=======================================================================
C
C
C     This is the top-level evolution program scfm.  Its tasks are:
C
C          1) to initialize file structures and global variables;
C          2) to input parameters and the initial system state;
C          3) to advance the state of the system for a given number
C             of timesteps;
C          4) to perform a diagnostic analysis of the system at
C             each time step (energy, angular momentum, etc.);
C          5) to periodically record the state of the system;
C          6) and to terminate the simulation and close data files.
C
C
C=======================================================================
C
C
C     Basic global variables/parameters:
C
C          ax,ay,az    : accelerations of bodies.
C          clm, dlm,   : radial functions used to evaluate expansions.
C          elm, flm
C          cputime     : cpu time (secs) used during the simulation.
C          cputime0    : cumulative cpu time at start of run.
C          cputime1    : cumulative cpu time at end of run.
C          dtime       : the timestep.
C          fixacc      : option to force conservation of linear
C                        momentum by setting acceleration of c.o.m.=0.
C          G           : the gravitational constant.
C          headline    : identification string for the run.
C          inptcoef    : option to read in expansion coefficients.
C          lmax        : number of angular eigenfunctions.
C          mass        : masses of bodies.
C          nbodies     : total number of bodies.
C          nbodsmax    : maximum number of bodies.
C          nmax        : number of radial eigenfunctions.
C          noutbod     : frequency of system record outputs.
C          noutlog     : frequency of outputs to log file.
C          nsteps      : number of time-steps to integrate the system.
C          one         : the constant 1.
C          onesixth    : the constant 1/6.
C          outpcoef    : option to write out expansion coefficients.
C          pi          : the constant pi.
C          pot         : potentials of bodies (self-gravity).
C          potext      : potentials of bodies (external field).
C          selfgrav    : option to turn off (.FALSE.) system self-
C                        gravity.
C          tnow        : current system time.
C          tpos        : current position time.
C          tvel        : current velocity time.
C          twoopi      : the constant 2./pi.
C          vx,vy,vz    : velocities of bodies.
C          x,y,z       : positions of bodies.
C          zeroeven    : option to zero out all even terms in the
C                        basis function expansions.
C          zeroodd     : option to zero out all odd terms in the
C                        basis function expansions.
C
C
C-----------------------------------------------------------------------
C
C   Definitions specific to input/output.
C
C          uterm, upars, ulog, ubodsin,   : logical i/o unit numbers.
C            ubodsout,utermfil,uoutcoef,
C            uincoef,ubodsel
C          parsfile, logfile, ibodfile,   : character names of files.
C            obodfile,termfile,outcfile,
C            incfile,elfile
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER icheck,n

C=======================================================================
C   Initialize state of the system.
C   -------------------------------
        CALL initsys
C            -------
	WRITE(uterm,*)'initsys'
	WRITE(utermfil,*)'initsys'
C   If tidal interaction turn on field.
C   -----------------------------------
        IF(david.AND.ntide.GT.0) CALL tidalstart
C                                     ----------
	WRITE(uterm,*)'tidalstart'
	WRITE(utermfil,*)'tidalstart'
C   Advance system state for a given number of steps.
C   -------------------------------------------------
	icheck=0
        DO 100 n=1,nsteps

           IF((tnow.GT.tend).OR.(icheck.GT.0))THEN
              CALL stepsys(nsteps)
C                  -------
              GO TO 110
	   ELSE 
              CALL stepsys(n)
C                  -------
           ENDIF
C
	   IF(david)CALL checkpos(n,icheck)
C                        --------
 100    CONTINUE

 110    CONTINUE
C   Terminate the simulation.
C   -------------------------
        CALL endrun
C            ------

	RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE accp_LH
C
C
C***********************************************************************
C
C
C     Subroutine to compute accelerations, potential, and density.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

        INTEGER k,l,m,n,lmin,lskip
        LOGICAL firstc
        REAL*8 anltilde,knl,sinth,sinmphi,cosmphi,phinltil,deltam0,
     &         gammln,arggam,coeflm,factrl,
     &         dblfact,ttemp5,ar,ath,aphi,temp3,temp4,
     &         temp5,temp6,plm,dplm,ultrasp,ultrasp1,ultraspt,clm,
     &         dlm,elm,flm,xi,costh,phi,r,twoalpha,c1,c2,c3,un,unm1,
     &         plm1m,plm2m,cosp,sinp

        DIMENSION ultrasp(0:nmax,0:lmax),
     &            ultraspt(0:nmax,0:lmax),ultrasp1(0:nmax,0:lmax),
     &            anltilde(0:nmax,0:lmax),dblfact(lmax+1),
     &            coeflm(0:lmax,0:lmax),
     &            twoalpha(0:lmax),c1(1:nmax,0:lmax),c2(1:nmax,0:lmax),
     &            c3(1:nmax),cosmphi(0:lmax),sinmphi(0:lmax),
     &            plm(0:lmax,0:lmax),dplm(0:lmax,0:lmax)

        DATA firstc/.TRUE./

        SAVE firstc,dblfact,anltilde,coeflm,lmin,
     &       lskip,twoalpha,c1,c2,c3

C=======================================================================

        IF(firstc) THEN

           firstc=.FALSE.

           dblfact(1)=1.

           DO 5 l=2,lmax
              dblfact(l)=dblfact(l-1)*(2.*l-1.)
 5         CONTINUE

           DO 20 n=0,nmax
              DO 10 l=0,lmax
                 knl=0.5*n*(n+4.*l+3.)+(l+1.)*(2.*l+1.)
                 anltilde(n,l)=-2.**(8.*l+6.)*FACTRL(n)*(n+2.*l+1.5)
                 arggam=2.*l+1.5
                 anltilde(n,l)=anltilde(n,l)*(EXP(GAMMLN(arggam)))**2
                 anltilde(n,l)=anltilde(n,l)/(4.*pi*knl*FACTRL(n+4*l+2))
 10           CONTINUE
 20        CONTINUE

           DO 25 l=0,lmax

              twoalpha(l)=2.0*(2.*l+1.5)

              DO 23 m=0,l
                 deltam0=2.
                 IF(m.EQ.0) deltam0=1.
                 coeflm(l,m)=(2.*l+1.)*deltam0*FACTRL(l-m)/FACTRL(l+m)
 23           CONTINUE
 25        CONTINUE

           DO 30 n=1,nmax
              c3(n)=1.0/(n+1.0)

              DO 27 l=0,lmax
                 c1(n,l)=2.0*n+twoalpha(l)
                 c2(n,l)=n-1.0+twoalpha(l)
 27           CONTINUE

 30        CONTINUE

           lskip=1
           IF(zeroodd.OR.zeroeven) lskip=2

           lmin=0
           IF(zeroeven) lmin=1

        ENDIF

        DO 60 l=0,lmax
           DO 50 m=0,l
              DO 40 n=0,nmax
                 sinsum(n,l,m)=0.0
                 cossum(n,l,m)=0.0
 40           CONTINUE
 50        CONTINUE
 60     CONTINUE

        DO 120 k=1,nbodies

           IF(ibound(k).GT.0)THEN
              r=SQRT(x(k)**2+y(k)**2+z(k)**2)
              costh=z(k)/r
              phi=ATAN2(y(k),x(k))
              xi=(r-1.)/(r+1.)
              
              DO 105 m=0,lmax
                 cosmphi(m)=COS(m*phi)           
                 sinmphi(m)=SIN(m*phi)           
 105          CONTINUE
              
              DO 113 l=0,lmax

                 ultrasp(0,l)=1.0
                 ultrasp(1,l)=twoalpha(l)*xi

                 un=ultrasp(1,l)
                 unm1=1.0
                 
                 DO 111 n=1,nmax-1
                    ultrasp(n+1,l)=(c1(n,l)*xi*un-c2(n,l)*unm1)*c3(n)
                    unm1=un
                    un=ultrasp(n+1,l)
 111             CONTINUE
                 
                 DO 112 n=0,nmax
                    ultraspt(n,l)=ultrasp(n,l)*anltilde(n,l)
 112             CONTINUE

 113          CONTINUE

              DO 1132 m=0,lmax
                 
                 plm(m,m)=1.0
                 IF(m.GT.0) plm(m,m)=(-1.)**m*dblfact(m)*SQRT(1.-
     &                costh*costh)**m
                 plm1m=plm(m,m)
                 plm2m=0.0
                 
                 DO 1131 l=m+1,lmax
                    plm(l,m)=(costh*(2.*l-1.)*plm1m-
     &                        (l+m-1.)*plm2m)/(l-m)
                    plm2m=plm1m
                    plm1m=plm(l,m)
 1131            CONTINUE

 1132         CONTINUE

              DO 118 l=lmin,lmax,lskip
                 
                 temp5=r**l/((1.+r)**(2*l+1))*mass(k)

                 DO 116 m=0,l

                    ttemp5=temp5*plm(l,m)*coeflm(l,m)
                    temp3=ttemp5*sinmphi(m)
                    temp4=ttemp5*cosmphi(m)

                    DO 114 n=0,nmax
                       sinsum(n,l,m)=sinsum(n,l,m)+temp3*ultraspt(n,l)
                       cossum(n,l,m)=cossum(n,l,m)+temp4*ultraspt(n,l)
 114                CONTINUE

 116             CONTINUE
 118          CONTINUE
           END IF
 120    CONTINUE

        
        DO 200 k=1,nbodies

           r=SQRT(x(k)**2+y(k)**2+z(k)**2)
           costh=z(k)/r
           phi=ATAN2(y(k),x(k))
           xi=(r-1.)/(r+1.)

           DO 130 m=0,lmax
              cosmphi(m)=COS(m*phi)           
              sinmphi(m)=SIN(m*phi)           
 130       CONTINUE

           pot(k)=0.0
           ar=0.0
           ath=0.0
           aphi=0.0

           DO 148 l=0,lmax

              ultrasp(0,l)=1.0
              ultrasp(1,l)=twoalpha(l)*xi
              ultrasp1(0,l)=0.0
              ultrasp1(1,l)=1.0

              un=ultrasp(1,l)
              unm1=1.0

              DO 144 n=1,nmax-1
                 ultrasp(n+1,l)=(c1(n,l)*xi*un-c2(n,l)*unm1)*c3(n)
                 unm1=un
                 un=ultrasp(n+1,l)
                 ultrasp1(n+1,l)=((twoalpha(l)+(n+1)-1.)*unm1-(n+1)*xi*
     &                    ultrasp(n+1,l))/(twoalpha(l)*(1.-xi*xi))
 144          CONTINUE

 148       CONTINUE

           DO 1482 m=0,lmax

              plm(m,m)=1.0
              IF(m.GT.0) plm(m,m)=(-1.)**m*dblfact(m)*SQRT(1.-
     &                            costh*costh)**m
              plm1m=plm(m,m)
              plm2m=0.0

              DO 1481 l=m+1,lmax
                 plm(l,m)=(costh*(2.*l-1.)*plm1m-(l+m-1.)*plm2m)/(l-m)
                 plm2m=plm1m
                 plm1m=plm(l,m)
 1481         CONTINUE

 1482      CONTINUE

           dplm(0,0)=0.0

           DO 1486 l=1,lmax

              DO 1484 m=0,l

                 IF(l.EQ.m) THEN
                    dplm(l,m)=l*costh*plm(l,m)/(costh*costh-1.0)
                 ELSE
                    dplm(l,m)=(l*costh*plm(l,m)-(l+m)*plm(l-1,m))/
     &                        (costh*costh-1.0)
                 ENDIF

 1484         CONTINUE
 1486      CONTINUE

           DO 190 l=lmin,lmax,lskip

              temp3=0.0
              temp4=0.0
              temp5=0.0
              temp6=0.0

              DO 180 m=0,l

                 clm=0.0
                 dlm=0.0
                 elm=0.0
                 flm=0.0

                 DO 150 n=0,nmax
                    clm=clm+ultrasp(n,l)*cossum(n,l,m)
                    dlm=dlm+ultrasp(n,l)*sinsum(n,l,m)
                    elm=elm+ultrasp1(n,l)*cossum(n,l,m)
                    flm=flm+ultrasp1(n,l)*sinsum(n,l,m)
 150             CONTINUE

                 temp3=temp3+plm(l,m)*(clm*cosmphi(m)+dlm*sinmphi(m))
                 temp4=temp4-plm(l,m)*(elm*cosmphi(m)+flm*sinmphi(m))
                 temp5=temp5-dplm(l,m)*(clm*cosmphi(m)+dlm*sinmphi(m))
                 temp6=temp6-m*plm(l,m)*(dlm*cosmphi(m)-clm*sinmphi(m))
 180          CONTINUE

              phinltil=r**l/((1.+r)**(2*l+1))
              pot(k)=pot(k)+temp3*phinltil
              ar=ar+phinltil*(-temp3*(l/r-(2.*l+1.)/
     &              (1.+r))+temp4*4.*(2.*l+1.5)/(1.+r)**2)
              ath=ath+temp5*phinltil
              aphi=aphi+temp6*phinltil

 190       CONTINUE

	   cosp=COS(phi)
	   sinp=SIN(phi)

           sinth=SQRT(1.-costh**2)
           ath= -sinth*ath/r
           aphi=aphi/(r*sinth)
           ax(k)=G*(sinth*cosp*ar+costh*cosp*ath-
     &           sinp*aphi)
           ay(k)=G*(sinth*sinp*ar+costh*sinp*ath+
     &           cosp*aphi)
           az(k)=G*(costh*ar-sinth*ath)
           pot(k)=pot(k)*G

 200    CONTINUE

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE accp_LHa
C
C
C***********************************************************************
C
C
C     Subroutine to compute accelerations, potential, and density.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

        INTEGER k
        REAL*8 r

C=======================================================================

        DO 10 k=1,nbodies
           r=SQRT(x(k)**2+y(k)**2+z(k)**2)
           ax(k)=ax(k)-G*x(k)/(r*(1.+r)**2)
           ay(k)=ay(k)-G*y(k)/(r*(1.+r)**2)
           az(k)=az(k)-G*z(k)/(r*(1.+r)**2)
           potext(k)=potext(k)-G/(1.+r)
 10     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE accpot(n)
C
C
C***********************************************************************
C
C
C     Subroutine to compute accelerations, potential, and density.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

        INTEGER i,n

C=======================================================================

        IF(selfgrav) THEN
           CALL accp_LH
C               -------
           IF(fixacc) CALL corracc
C                          -------
           IF(david) CALL accp_DAVID(n)
C                         ----------
        ENDIF

        IF(.NOT.selfgrav) THEN
           DO 10 i=1,nbodies
              ax(i)=0.0
              ay(i)=0.0
              az(i)=0.0
              pot(i)=0.0
 10        CONTINUE

           IF(david) CALL accp_DAVID(n)
C                         ----------
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE checkinp
C
C
C***********************************************************************
C
C
C     Subroutine to check consistency of input parameters and data,
C     output warnings to the terminal and/or log file, and terminate
C     the simulation if necessary.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

C=======================================================================

        IF(nsteps.LT.0.OR.nsteps.GT.10000000)
     &     CALL terror(' input error for parameter nsteps ')
C               ------

        IF(noutbod.LT.0)
     &     CALL terror(' input error for parameter noutbod ')
C               ------

        IF(noutlog.LT.0)
     &     CALL terror(' input error for parameter noutlog ')
C               ------

        IF(dtime.LE.-1.e20.OR.dtime.GT.1.e20)
     &     CALL terror(' input error for parameter dtime ')
C               ------

        IF(G.LE.0.0.OR.G.GT.1.e20)
     &     CALL terror(' input error for parameter G ')
C               ------

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE corracc
C
C
C***********************************************************************
C
C
C     Subroutine to correct accelerations so that the center of
C     mass remains fixed at the origin.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i
        REAL*8 axcm,aycm,azcm,mtot

C=======================================================================

        axcm=0.0
        aycm=0.0
        azcm=0.0
        mtot=0.0

        DO 10 i=1,nbodies
           mtot=mtot+mass(i)
           axcm=axcm+mass(i)*ax(i)
           aycm=aycm+mass(i)*ay(i)
           azcm=azcm+mass(i)*az(i)
 10     CONTINUE

        axcm=axcm/mtot
        aycm=aycm/mtot
        azcm=azcm/mtot

        DO 20 i=1,nbodies
           ax(i)=ax(i)-axcm
           ay(i)=ay(i)-aycm
           az(i)=az(i)-azcm
 20     CONTINUE

        RETURN
        END

C***********************************************************************
C
C
                         SUBROUTINE corrvel(rc)
C
C
C***********************************************************************
C
C
C     Subroutine to synchronize particle coordinates when outputing
C     particle data to body data file or when computing diagnostics.
C  
C
C=======================================================================

        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 rc
        INTEGER p
        REAL*8 rcsign

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------
  
        IF(rc.EQ.'correct') THEN
           rcsign=-1.
        ELSE
           rcsign=1.
        ENDIF

        DO 10 p=1,nbodies
           vx(p)=vx(p)+rcsign*ax(p)*0.5*dtime
           vy(p)=vy(p)+rcsign*ay(p)*0.5*dtime
           vz(p)=vz(p)+rcsign*az(p)*0.5*dtime
 10     CONTINUE

C   Update velocity time, system time.
C   ----------------------------------
        tvel=tvel+rcsign*0.5*dtime
        tnow=tvel

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE endrun
C
C
C***********************************************************************
C
C
C     Subroutine to end the simulation.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------

C        REAL*8 rtc
        REAL tarray(2),etime
C=======================================================================

c	CALL SECOND(sec)
c        cputime1=RTC()
	cputime1=ETIME(tarray)

        CALL stopout
C            -------

        RETURN
        END
C***********************************************************************
C
C
        FUNCTION FACTRL(N)
C
C
C***********************************************************************
C
C
C     A function to compute factorials.  (From numerical recipes.)
C
C
C=======================================================================

        INTEGER n,ntop,j
        REAL*8 factrl,a,gammln,arggam

        DIMENSION A(33)

        DATA NTOP,A(1)/0,1./

        IF (N.LT.0) THEN
          PAUSE 'negative factorial'
        ELSE IF (N.LE.NTOP) THEN
          FACTRL=A(N+1)
        ELSE IF (N.LE.32) THEN
          DO 11 J=NTOP+1,N
            A(J+1)=J*A(J)
11        CONTINUE
          NTOP=N
          FACTRL=A(N+1)
        ELSE
          arggam=n+1.
          FACTRL=EXP(GAMMLN(arggam))
        ENDIF

        RETURN
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
        SUBROUTINE inbods
C
C
C***********************************************************************
C
C
C     Subroutine to read phase coordinates.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

        INTEGER i
        CHARACTER*8 filepar
	CHARACTER*6 filename

        filepar=ibodfile
        IF(ntide.GT.0)THEN
	   filename=filepar(1:2)
	ELSE
	   filename=outfile(1:4)//filepar(1:2)
	END IF

c        OPEN(ubodsin,FILE=filename,STATUS='OLD')
       OPEN(ubodsin,FILE='BI',
C       OPEN(ubodsin,FILE='/export/radegund/kvj/video/data/test3/scfbi',
     &  STATUS='OLD')

        READ(ubodsin,*) nbodies,tnow
C	nbodies=nbodies/30.
        DO 10 i=1,nbodies
           READ(ubodsin,*) mass(i),x(i),y(i),z(i),vx(i),vy(i),vz(i)
C	   mass(i)=mass(i)*30.0d0
 10     CONTINUE

        CLOSE(ubodsin)

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initpars
C
C
C***********************************************************************
C
C
C     Subroutine to initialize system parameters that depend on
C     either the input data or defined PARAMETERS.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

        INTEGER i,istring
C
        CHARACTER*4 sstring
        CHARACTER*10 nstring

        DATA nstring/'0123456789'/
	SAVE nstring

C=======================================================================

C   Initialize misc. useful numbers.
C   --------------------------------
        one=1.0
        two=2.0
        pi=4.0*ATAN(one)
        twoopi=2./pi
        onesixth=1./6.
        tiny=1.e-30
        zero=0.0

	dtime=dt0
C Define relevant subdirectory.
	sstring(1:1)='R'
        sstring(2:2)=nstring(1+nproc/100:1+nproc/100)
        istring=1+MOD(nproc,100)/10
        sstring(3:3)=nstring(istring:istring)
        istring=1+MOD(nproc,10)
        sstring(4:4)=nstring(istring:istring)
        outfile=sstring(1:4)

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initsys
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the state of the system.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,ival,ierr,iclktck
c        REAL*8 rtc
	REAL etime,tarray(2)
C=======================================================================

C   Begin timing. First find 'tick' value.
C   -------------
c	ival=0
c	iclktck=0
c	CALL PXFCONST('CLK_TCK',iclktck,ierr)
C            --------
c	IF(ierr.NE.0)CALL terror('PXFCONST error')
C                         ------
c	CALL PXFSYSCONF(iclktck,ival,ierr)
C            ----------
c	IF(ierr.NE.0)CALL terror('PXFSYSCONF error')
C                         ------
c	tick=DBLE(REAL(ival))	
c        cputime0=RTC()
	tick=1.
c	CALL SECOND(sec)
	cputime0=ETIME(tarray)
C
        CALL inparams
C            --------
        CALL initpars
C            --------
        CALL inbods
C            ------
        tpos=tnow
        tvel=tnow

C All particles initially bound.
        DO 10 i=1,nbodies
           ibound(i)=1
	   tub(i)=zero
	   lmflag(i)=zero
 10     CONTINUE

        CALL startout
C            --------
        CALL checkinp
C            --------
        CALL accpot(0)
C            ------
        IF(david) CALL frame(0)
C                      -----
        CALL outstate(0)
C            --------
        CALL initvel
C            -------

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initvel
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the velocities of the bodies for the
C     initial timestep.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------
  
        DO 10 p=1,nbodies
           vx(p)=vx(p)+0.5*dtime*ax(p)
           vy(p)=vy(p)+0.5*dtime*ay(p)
           vz(p)=vz(p)+0.5*dtime*az(p)
 10     CONTINUE

C   Update velocity time, system time.
C   ----------------------------------
        tvel=tvel+0.5*dtime
        tnow=tvel

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE inparams
C
C
C***********************************************************************
C
C
C     Subroutine to read in parameters.
C
C     Input parameters:
C
C        headline  : identification string for the run.
C        nsteps    : number of timesteps.
C        noutbod   : output system state once every nsteps/noutbod 
C                    steps.
C        noutlog   : output logfile data once every nsteps/noutlog
C                    steps.
C        dtime     : the timestep.
C        G         : value of gravitational constant, in appropriate
C                    units.
C        selfgrav  : option to turn off (.FALSE.) system self-gravity.
C        inptcoef  : option to read-in expansion coefficients.
C        outpcoef  : option to write-out expansion coefficients.
C        zeroodd   : option to zero all odd terms in the expansion.
C        zeroeven  : option to zero all even terms in the expansion.
C        fixacc    : option to force conservation of linear
C                    momentum by subtracting acceleration of c.o.m.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

        CHARACTER *1 pcomment
        CHARACTER*8 filepar,filename
C=======================================================================


        filepar=parsfile
        filename=filepar(1:6)

        OPEN(UNIT=upars,FILE=filename,STATUS='OLD')

C   Read parameters, close the file.
C   --------------------------------

        READ(upars,'(a)') pcomment
        READ(upars,'(a)') headline
	READ(upars,*) nsteps
	READ(upars,*) noutlog
	READ(upars,*) noutbod
        READ(upars,*) dt0
        READ(upars,*) G
        READ(upars,*) selfgrav
        READ(upars,*) inptcoef
        READ(upars,*) outpcoef
        READ(upars,*) zeroodd
        READ(upars,*) zeroeven
        READ(upars,*) fixacc
        READ(upars,*) david
	READ(upars,*) dfric
        READ(upars,*) disk
	READ(upars,*) flat
	READ(upars,*) growing
	READ(upars,*) flattening
	READ(upars,*) nsort
	READ(upars,*) ntide
        CLOSE(UNIT=upars)
 
        RETURN
        END
C***********************************************************************
C
C
                     SUBROUTINE iocoef
C
C
C***********************************************************************
C
C
C     Subroutine to input and output expansion coefficients.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

        INTEGER n,l,m
        LOGICAL firstc
        REAL*8 tt
	CHARACTER*8 filename
	CHARACTER*8 filepar

        DATA firstc/.TRUE./

        SAVE firstc

C=======================================================================

        IF(firstc) THEN

           firstc=.FALSE.

	   filepar=outcfile
	   filename=outfile(1:4)//filepar(1:4)
           IF(outpcoef) OPEN(uoutcoef,FILE=filename,STATUS='NEW')
C	   filepar=incfile
C	   filename=outdir(1:7)//filepar(1:8)
C           IF(inptcoef) OPEN(uincoef,FILE=filename,STATUS='OLD')

        ENDIF

        IF(outpcoef) THEN

           WRITE(uoutcoef,100) tnow

           DO 30 n=0,nmax
              DO 20 l=0,lmax
                 DO 10 m=0,l
                    WRITE(uoutcoef,100) sinsum(n,l,m),cossum(n,l,m)
 10              CONTINUE
 20           CONTINUE
 30        CONTINUE

 100       FORMAT(1x,10(1pe22.13))

        ENDIF

        IF(inptcoef) THEN

           READ(uincoef,*) tt

           IF(tt.NE.tnow) CALL terror(' input error in iocoef ')
C                              ------

           DO 130 n=0,nmax
              DO 120 l=0,lmax
                 DO 110 m=0,l
                    READ(uoutcoef,*) sinsum(n,l,m),cossum(n,l,m)
 110             CONTINUE
 120          CONTINUE
 130       CONTINUE

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE outbods(n)
C
C
C***********************************************************************
C
C
C     Subroutine to output phase coordinates.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

        CHARACTER*8 sstring
        CHARACTER*40 filename
        CHARACTER*10 nstring

        INTEGER i,istring,nsnap,n

        SAVE nsnap,nstring

        DATA nsnap/0/,nstring/'0123456789'/

C=======================================================================

	IF(n.EQ.0)nsnap=0
        nsnap=nsnap+1

	sstring(1:1)='S'
        sstring(2:2)=nstring(1+nsnap/100:1+nsnap/100)
        istring=1+MOD(nsnap,100)/10
        sstring(3:3)=nstring(istring:istring)
        istring=1+MOD(nsnap,10)
        sstring(4:4)=nstring(istring:istring)
        filename=outfile(1:4)//sstring(1:4)

	OPEN(UNIT=ubodsout,FILE=filename,STATUS='NEW')
	WRITE(ubodsout,20) nbodies,tnow
 20     FORMAT(1x,1i6,1pe14.6)

        DO 30 i=1,nbodies
           IF(david)THEN
	      IF(selfgrav)
     &              WRITE(ubodsout,111)mass(i),x(i),y(i),z(i),
     &              vx(i)-vxframe,vy(i)-vyframe,vz(i)-vzframe,
     &                            tub(i),ep(i),potext(i),lmflag(i)
	      IF(.NOT.selfgrav)WRITE(ubodsout,111) mass(i),
     &            x(i),y(i),z(i),vx(i),vy(i),vz(i),potext(i)
           ELSE
              WRITE(ubodsout,111) mass(i),x(i),y(i),z(i),
     &                            vx(i),vy(i),vz(i),pot(i)
           END IF
 30     CONTINUE
C	ENDIF
 111    FORMAT(1x,11(1pe14.6))

        CLOSE(ubodsout)
        
        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE outlog(n)
C
C
C***********************************************************************
C
C
C     Subroutine to output phase coordinates.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

        INTEGER i,n
        LOGICAL firstc
        REAL*8 rtc
	REAL etime,tarray(2)
        REAL*8 lxtot,lytot,lztot,mtot,vxcm,vycm,vzcm,etot,ektot,epext,
     &         m2tw,t1,clausius,m2claus,cpux,xcm,ycm,zcm,epselfg
	CHARACTER*7 filename
	CHARACTER*8 filepar

        DATA firstc/.TRUE./

        SAVE firstc

C=======================================================================

C        IF(firstc) THEN
C           firstc=.FALSE.
	IF(n.EQ.0)THEN
           cputime=cputime0

	   filepar=logfile
	   filename=outfile(1:4)//filepar(1:3)
           OPEN(UNIT=ulog,FILE=filename,STATUS='UNKNOWN')
        ENDIF

C        DO 10 i=1,2*nsteps
C           READ(ulog,120,end=20) t1
C 10     CONTINUE
 
C 20     CONTINUE

        mtot=0.0
        etot=0.0
        ektot=0.0
        epext=0.0
        epselfg=0.0

	IF(selfgrav)THEN
        IF(.NOT.david) THEN
           xcm=0.0
           ycm=0.0
           zcm=0.0
           vxcm=0.0
           vycm=0.0
           vzcm=0.0
           lxtot=0.0
           lytot=0.0
           lztot=0.0
           clausius=0.0
        END IF

        DO 30 i=1,nbodies

           mtot=mtot+mass(i)
           epext=epext+mass(i)*potext(i)
           epselfg=epselfg+0.5*mass(i)*pot(i)
           ektot=ektot+0.5*mass(i)*((vx(i))**2+
     &         (vy(i))**2+(vz(i))**2)
           IF(.NOT.david) THEN
              xcm=xcm+mass(i)*x(i)
              ycm=ycm+mass(i)*y(i)
              zcm=zcm+mass(i)*z(i)
              vxcm=vxcm+mass(i)*vx(i)
              vycm=vycm+mass(i)*vy(i)
              vzcm=vzcm+mass(i)*vz(i)
              lxtot=lxtot+mass(i)*(y(i)*vz(i)-z(i)*vy(i))
              lytot=lytot+mass(i)*(z(i)*vx(i)-x(i)*vz(i))
              lztot=lztot+mass(i)*(x(i)*vy(i)-y(i)*vx(i))
              clausius=clausius+mass(i)*(x(i)*ax(i)+y(i)*ay(i)
     &                                  +z(i)*az(i))
           END IF
 30     CONTINUE
 99     FORMAT(5(1pe14.4))
        etot=ektot+epext+epselfg

        IF(david)THEN
           xcm=xframe
           ycm=yframe
           zcm=zframe
           vxcm=vxframe
           vycm=vyframe
           vzcm=vzframe
           IF(ABS(eprem).GT.0.)THEN
	      m2tw= -2.*ekrem/eprem
	   ELSE
	      m2tw=0.0
	   END IF
        ELSE
           xcm=xcm/mtot
           ycm=ycm/mtot
           zcm=zcm/mtot
           vxcm=vxcm/mtot
           vycm=vycm/mtot
           vzcm=vzcm/mtot
           m2claus= -2.*ektot/clausius
	   m2tw=-2.*ektot/epselfg
        END IF

        WRITE(ulog,120) tnow,dtime
        WRITE(ulog,120) mtot,mrem
        WRITE(ulog,120) xcm,ycm,zcm
        WRITE(ulog,120) vxcm,vycm,vzcm
        WRITE(ulog,120) ektot,epext,epselfg,etot
C
        IF(david)THEN
           WRITE(ulog,120) ekrem,eprem,m2tw
        ELSE
           WRITE(ulog,120) m2tw,clausius,m2claus
        END IF
        ELSE
           DO 40 i=1,nbodies
              mtot=mtot+mass(i)
              epext=epext+mass(i)*potext(i)
              ektot=ektot+0.5*mass(i)*((vx(i))**2+
     &         (vy(i))**2+(vz(i))**2)
 40        CONTINUE
      	   etot=ektot+epext
           WRITE(ulog,120)tnow,dtime
           WRITE(ulog,120)mtot,mrem
	   WRITE(ulog,120)ektot,epext,etot
        END IF
C
c	CALL SECOND(sec)
	cpux=ETIME(tarray)
c        cpux=RTC()

        WRITE(ulog,120) (cpux-cputime)/tick
        WRITE(ulog,130)

 120    FORMAT(5(1pe18.10))
 130    FORMAT(/)


        cputime=cpux

        RETURN
        END
C***********************************************************************
C
C
                        SUBROUTINE outstate(n)
C
C
C***********************************************************************
C
C
C     Subroutine to output information about the system state to
C     the log and body data files.
C
C
C=======================================================================
 
        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n,i

C=======================================================================

        CALL outterm(' step completed: ',n)
C            -------

        IF(n.EQ.0) THEN

           IF(david)THEN
C	      CALL stepsize
C                  --------
              IF(ntide.EQ.0)CALL findrem(n)
C                                -------
                    IF(.NOT.selfgrav)THEN
                       DO 10 i=1,nbodies
                          x(i)=x(i)+xframe
                          y(i)=y(i)+yframe
                          z(i)=z(i)+zframe
 10                    CONTINUE
                       xframe=zero
                       yframe=zero
                       zframe=zero
                    END IF
           ELSE
              CALL outbods(n)
C                  -------
              CALL iocoef
C                  ------
           END IF
C
           CALL outlog(n)
C               ------
        ELSE

           IF(MOD(n,noutlog).EQ.0) THEN

              CALL corrvel('correct')
C                  -------
              IF(david)THEN
	         IF(selfgrav)THEN
C                    CALL stepsize
C                        --------
		    CALL findrem(n)
C                        -------
		    IF(MOD(n,noutbod).EQ.0)CALL outbods(n)
C                                                -------
                    IF(.NOT.selfgrav)THEN
                       DO 20 i=1,nbodies
                          x(i)=x(i)+xframe
                          y(i)=y(i)+yframe
                          z(i)=z(i)+zframe
 20                    CONTINUE
                       xframe=zero
                       yframe=zero
                       zframe=zero
                    END IF
	         ELSE
                    IF(n.EQ.nsteps)CALL outbods(n)
C                                       -------
                 END IF
              ELSE IF(n.EQ.nsteps) THEN
                 CALL outbods(n)
C                     -------
C                 CALL iocoef
C                     ------
              END IF
C
              CALL outlog(n)
C                  ------
              CALL corrvel('reset  ')
C                  -------
           ENDIF

        ENDIF
 
        RETURN
        END
C***********************************************************************
C
C
                      SUBROUTINE outterm(message,n)
C
C
C***********************************************************************
C
C
C     Subroutine to output a message to the terminal and to the
C     terminal emulation file.
C
C
C=======================================================================
 
        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message
        INTEGER n


C   Write the message.
C   ------------------

        IF(n.GE.0) THEN
           WRITE(uterm,*) nproc,message,n
           WRITE(utermfil,*) message,n
        ELSE
           WRITE(utermfil,40)
           WRITE(utermfil,50) message 
           WRITE(utermfil,40)
        ENDIF

 40     FORMAT(/,1x,72('*'))
 50     FORMAT(/,a)


        RETURN
        END

C***********************************************************************
C
C
                          SUBROUTINE startout
C
C
C***********************************************************************
C
C
C     Subroutine to open disk files for subsequent input/output.
C
C
C=======================================================================
 
        INCLUDE 'scfltr.h'

	CHARACTER*7 filename
	CHARACTER*8 filepar

C=======================================================================
 
C   Create terminal emulation file.
C   -------------------------------
        filepar=termfile
        filename=outfile(1:4)//filepar(1:3)
        OPEN(UNIT=utermfil,FILE=filename,STATUS='NEW')
        WRITE(utermfil,*) ' Start of output '

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
C
C     Subroutine to advance the positions of the bodies for a timestep.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------

        DO 10 p=1,nbodies
           x(p)=x(p)+vx(p)*dtime
           y(p)=y(p)+vy(p)*dtime
           z(p)=z(p)+vz(p)*dtime
 10     CONTINUE

C   Update position time, system time.
C   ----------------------------------
        tpos=tpos+dtime
        tnow=tpos

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE stepsys(n)
C
C
C***********************************************************************
C
C
C     Subroutine to advance the state of the system by one timestep.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n

C=======================================================================

        CALL steppos
C            -------
        IF (david.AND.selfgrav) CALL frame(n)
C                                    -----
        CALL accpot(n)
C            ------
        CALL stepvel
C            -------

        CALL outstate(n)
C            --------

        RETURN
        END

C***********************************************************************
C
C
                         SUBROUTINE stepvel
C
C
C***********************************************************************
C
C
C     Subroutine to advance the velocities of the bodies for timestep.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p

C=======================================================================

C   Loop over all velocity components for all bodies.
C   -------------------------------------------------

        DO 10 p=1,nbodies
           vx(p)=vx(p)+ax(p)*dtime
           vy(p)=vy(p)+ay(p)*dtime
           vz(p)=vz(p)+az(p)*dtime
 10     CONTINUE

C   Update velocity time, system time.
C   ----------------------------------
        tvel=tvel+dtime
        tnow=tvel

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE stopout
C
C
C***********************************************************************
C
C
C     Subroutine to end output.
C
C
C=======================================================================
 
        INCLUDE 'scfltr.h'

        INTEGER i
        REAL*8 t1

C=======================================================================
C        DO 10 i=1,2*nsteps
C           READ(ulog,120,end=20) t1
C 10     CONTINUE
 
C 20     CONTINUE

        WRITE(ulog,30) (cputime1-cputime0)/tick
 30     FORMAT(//,15x,' Total cpu time used (secs) =',1pe15.7)

 120    FORMAT(15(1pe13.5))

        CLOSE(ulog)

        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE terror(message)
C
C
C***********************************************************************
C
C
C     Subroutine to terminate the program as the result of a fatal
C     error, close the output files, and dump timing information.
C
C
C=======================================================================

        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message
        INTEGER ierror
        REAL*8 rtc
        REAL etime,tarray(2)

C=======================================================================

C   Write error message to the log file and to the terminal.
C   --------------------------------------------------------
        ierror=-1

        CALL outterm(message,ierror)
C            -------

C-----------------------------------------------------------------------
C   Stop timing, output timing data, close files, terminate the
C   simulation.
C-----------------------------------------------------------------------
c        cputime1=RTC()
c	CALL SECOND(sec)
	cputime1=ETIME(tarray)

        STOP
        END

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

C***********************************************************************
C
C
                          SUBROUTINE accp_DAVID(n)
C
C       
C***********************************************************************

        INCLUDE 'scfltr.h'

C   Declaration of local variables.
C   -------------------------------
        INTEGER i,n
        REAL*8 r2,z2,thrad,rad,tsrad,sqz2b2,tdr,
     &         tdz,phim,phis,phih,dx,dy,dz,
     &         tu,gee,msun,cmpkpc,secperyr,
     &         a,b,rs,vcirc2,GMs,GM,growth,thtemp,
     &         x2,y2,thx,thy,thz,rpart,qx,qy,qz
     &         rothalo,rc1,rc2,rc3,x1,y1,z1
        PARAMETER(gee=6.67e-8,msun=1.989e33,cmpkpc=3.085678e21,
     &            secperyr=3.1536d7)

        LOGICAL firstc
        DATA firstc/.TRUE./
        SAVE firstc,a,b,GM,GMs,rs

C -----------------------------------------------------------------------
C       Three Component Galaxy model:
C
C               Miyamoto-Nagai Disk + Hernquist Spheroid + Logarithmic Halo
C
C
C       Features:
C
C               Nearly flat rotation curve from 1 - 30 kpc
C               Disk z scale height of 200 pc
C               kappa_z has correct radial dependance
C
C
C       Mass inside solar circle: (roughly)  
C
C
C               Disk 42%
C               Spheroid 28%
C               Halo 30%
C --------------------------------------------------------------
C Convert parameters to appropriate units....
        
C        IF(firstc)THEN
C           firstc=.FALSE.
	IF(n.EQ.0)THEN
C
           tu=SQRT((cmpkpc*ru)**3/(msun*mu*gee))
           vu=(ru*cmpkpc*1.e-5)/tu
C
	   tu=tu/secperyr/1.e9
	   WRITE(utermfil,*)ru,vu,tu
C
           tend=tend/tu
	   dtime=tend/REAL(nsteps)
	   WRITE(utermfil,*)dtime
	   WRITE(uterm,*)dtime
           strength=1.
           IF(ntide.GT.0)THEN
              strength=0.
C
              xframe=xframe/ru
              yframe=yframe/ru
              zframe=zframe/ru
              vxframe=vxframe/vu
              vyframe=vyframe/vu
              vzframe=vzframe/vu
           END IF
C
           DO 5 i=1,nbodies
              vx(i)=vx(i)+vxframe
              vy(i)=vy(i)+vyframe
              vz(i)=vz(i)+vzframe
 5         CONTINUE
C
           a=6.5/ru
           b=.26/ru
           c=dhalo/ru
           rs=.7/ru
           vcirc2=220.**2/vu/vu
C           vh2 = 0.325*vcirc2
           GM = 8.887*vcirc2/ru
           GMs = 3.0*vcirc2/ru
	   vh2 = vhalo*vhalo/vu/vu
        END IF

      IF(strength.EQ.one)THEN
C Find dynamical friction vector
         IF((dfric).AND.(mfric.GT.zero))CALL dynfric(n)
C                                       -------
      ENDIF

C Grow the halo
	growth=1.
	IF(growing)THEN
	   growth=0.5
	   IF(strength.EQ.1.)growth=(1.+tnow/tend)/2.
	ENDIF
C And flatten it..
	qx=q1
	qy=q2
	qz=q3
C Rotation angle of halo- set in RUNPAR for ease of scripting
        rothalo=drlphi*3.14159/180.
C Rotation constants
        rc1=cos(rothalo)*cos(rothalo)/qx/qx 
     & + sin(rothalo)*sin(rothalo)/qy/qy
        rc2=cos(rothalo)*cos(rothalo)/qy/qy
     & + sin(rothalo)*sin(rothalo)/qx/qx
        rc3=2*sin(rothalo)*cos(rothalo)*(1/qx/qx -1/qy/qy)

C       q=0.75
C	IF(flattening)THEN
C	   IF(strength.EQ.1.)q=(1.-0.25*tnow/tend)
C	ENDIF
C Triaxiality
c	p=0.95

        DO 10 i=1,nbodies
C
	   rpart=SQRT(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
C     compute acceleration due to spheroid
C
           x2 = (x(i)+xframe)*(x(i)+xframe)
           y2 = (y(i)+yframe)*(y(i)+yframe)
           z2 = (z(i)+zframe)*(z(i)+zframe)
           r2=x2+y2
           rad = sqrt(r2+z2)
           tsrad = GMs/(rad+rs)**2/rad
           phis = -GMs/(rad+rs)

           ax(i) = ax(i) -strength*(tsrad)*(x(i)+xframe)
           ay(i) = ay(i) -strength*(tsrad)*(y(i)+yframe)
           az(i) = az(i) -strength*(tsrad)*(z(i)+zframe)

C     
C
C      compute acceleration due to halo
           x1 = (x(i)+xframe)
           y1 = (y(i)+yframe)
           z1 = (z(i)+zframe)
	   thtemp=growth*vh2/(x2*rc1+y2*rc2+x1*y1*rc3+z2/qz/qz+c*c)
	   thx=thtemp*(2*rc1*x1+rc3*y1)
           thy=thtemp*(2*rc2*y1+rc3*x1)
           thz=thtemp*(2*z1/qz/qz)
           phih=growth*vh2*LOG(x2*rc1+y2*rc2+rc3*x1*y1+z2/qz/qz+c*c)

           ax(i) = ax(i) -strength*(thx)
           ay(i) = ay(i) -strength*(thy)
           az(i) = az(i) -strength*(thz)
C
C
           potext(i)= strength*(phis + phih)
C
C      compute acceleration due to disk
C
           IF(disk) THEN
C
              sqz2b2 = sqrt(z2 + b*b)
              tdr = GM/(r2 + (a + sqz2b2)**2)**1.5
              tdz = tdr*(a/sqz2b2 + 1.)
              phim = -GM/sqrt(r2+(a+sqz2b2)**2)
C
              ax(i)=ax(i)-strength*tdr*(x(i)+xframe)                         
              ay(i)=ay(i)-strength*tdr*(y(i)+yframe)                          
              az(i)=az(i)-strength*tdz*(z(i)+zframe) 
              potext(i)=potext(i)+strength*phim
           END IF
C Dynamical friction on particles within 2 tidal radii
        IF((dfric).AND.(strength.EQ.one))THEN
C        IF((dfric).AND.(strength.EQ.one).AND.
C     &                     (rpart.LE.4.0*rtide))THEN
C        IF((dfric).AND.(strength.EQ.one).AND.
C     &                     (tub(i).EQ.zero))THEN
            ax(i)=ax(i)+dfx
            ay(i)=ay(i)+dfy
            az(i)=az(i)+dfz
        ENDIF

 10     CONTINUE
                       
        RETURN
        END

C***********************************************************************
C
C
      SUBROUTINE accp_FLAT(rad,r2,i,dx,dy,dz,phih)
C
C       
C***********************************************************************

        INCLUDE 'scfltr.h'

        INTEGER l,n,nhalo,lhalo,lskip,i
        PARAMETER(nhalo=4,lhalo=4)
        LOGICAL firstc
        REAL*8 gmor,rh,phih,dx,dy,dz,rad,r2,r
        REAL*8 knl,sinth,phinltil,
     &         gammln,arggam,factrl,
     &         ar,ath,temp3,temp4,
     &         temp5,clm,
     &         elm,xi,costh,cosp,sinp,un,unm1,
     &         plm1,plm2

        REAL*8  ultrasp(0:nhalo,0:lhalo),
     &          ultrasp1(0:nhalo,0:lhalo),
     &          anltilde(0:nhalo,0:lhalo),
     &          twoalpha(0:lhalo),c1(1:nhalo,0:lhalo),
     &          c2(1:nhalo,0:lhalo),
     &          c3(1:nhalo),
     &          plm(0:lhalo),dplm(0:lhalo),
     &          cosum(0:nhalo,0:lhalo)

        DATA firstc/.TRUE./

        SAVE firstc,anltilde,
     &       lskip,twoalpha,c1,c2,c3,cosum,rh


        IF(firstc) THEN
           firstc=.FALSE.

           DO 20 n=0,nhalo
              DO 10 l=0,lhalo
                 knl=0.5*n*(n+4.*l+3.)+(l+1.)*(2.*l+1.)
                 anltilde(n,l)=-2.**(8.*l+6.)*FACTRL(n)*(n+2.*l+1.5)
                 arggam=2.*l+1.5
                 anltilde(n,l)=anltilde(n,l)*(EXP(GAMMLN(arggam)))**2
                 anltilde(n,l)=anltilde(n,l)/(4.*pi*knl*FACTRL(n+4*l+2))
 10           CONTINUE
 20        CONTINUE
           DO 25 l=0,lhalo
              twoalpha(l)=2.0*(2.*l+1.5)
 25        CONTINUE
           DO 30 n=1,nhalo
              c3(n)=1.0/(n+1.0)
              DO 27 l=0,lhalo
                 c1(n,l)=2.0*n+twoalpha(l)
                 c2(n,l)=n-1.0+twoalpha(l)
 27           CONTINUE
 30        CONTINUE

           lskip=2
C Hernquist oblate model, q=.7, mass=2e12, scale length 80kpc
           gmor=4.296d-6*2.d12/80./vu/vu
           rh=80./ru

           cosum(0,0)=-1.052243E+00*gmor
           cosum(0,2)=8.055071E-01*gmor
           cosum(0,4)=-1.194018E+00*gmor
           cosum(1,0)=1.882629E-02*gmor
           cosum(1,2)=2.759267E-02*gmor
           cosum(1,4)=-5.932488E-02*gmor
           cosum(2,0)=-7.688231E-04*gmor
           cosum(2,2)=2.033584E-02*gmor
           cosum(2,4)=-3.041664E-02*gmor
           cosum(3,0)=4.084624E-05*gmor
           cosum(3,2)=1.986118E-03*gmor
           cosum(3,4)=-3.922398E-03*gmor
           cosum(4,0)=-2.470607E-06*gmor
           cosum(4,2)=1.829174E-03*gmor
           cosum(4,4)=-2.362483E-03*gmor

        ENDIF

C HALO
        r=rad/rh
        costh=(z(i)+zframe)/rad
        xi=(r-1.)/(r+1.)
        cosp=(x(i)+xframe)/SQRT(r2)
        sinp=(y(i)+yframe)/SQRT(r2)

        phih=0.0
        ar=0.0
        ath=0.0

        DO 148 l=0,lhalo
           ultrasp(0,l)=1.0
           ultrasp(1,l)=twoalpha(l)*xi
           ultrasp1(0,l)=0.0
           ultrasp1(1,l)=1.0
           un=ultrasp(1,l)
           unm1=1.0
           DO 144 n=1,nhalo-1
              ultrasp(n+1,l)=(c1(n,l)*xi*un-c2(n,l)*unm1)*c3(n)
              unm1=un
              un=ultrasp(n+1,l)
              ultrasp1(n+1,l)=((twoalpha(l)+(n+1)-1.)*unm1-(n+1)*xi*
     &             ultrasp(n+1,l))/(twoalpha(l)*(1.-xi*xi))
 144       CONTINUE
 148    CONTINUE

        plm(0)=1.
        plm1=1.0
        plm2=0.0
        DO 1481 l=1,lhalo
           plm(l)=(costh*(2.*l-1.)*plm1-(l-1.)*plm2)/l
           plm2=plm1
           plm1=plm(l)
 1481   CONTINUE
        
        dplm(0)=0.0
        DO 1486 l=1,lhalo
           dplm(l)=(l*costh*plm(l)-(l)*plm(l-1))/
     &          (costh*costh-1.0)
 1486   CONTINUE

        DO 190 l=0,lhalo,lskip
           clm=0.0
           elm=0.0
           DO 150 n=0,nhalo
              clm=clm+ultrasp(n,l)*cosum(n,l)
              elm=elm+ultrasp1(n,l)*cosum(n,l)
 150       CONTINUE
           temp3=plm(l)*clm
           temp4=-plm(l)*elm
           temp5=-dplm(l)*clm
           phinltil=r**l/((1.+r)**(2*l+1))
           phih=phih+temp3*phinltil
           ar=ar+phinltil*(-temp3*(l/r-(2.*l+1.)/
     &          (1.+r))+temp4*4.*(2.*l+1.5)/(1.+r)**2)
           ath=ath+temp5*phinltil
 190    CONTINUE
        
        sinth=SQRT(1.-costh**2)
        ath= -sinth*ath/r

        dx=(sinth*cosp*ar+costh*cosp*ath)/rh
        dy=(sinth*sinp*ar+costh*sinp*ath)/rh
        dz=(costh*ar-sinth*ath)/rh

        RETURN
        END
C *********************************************************************
C
C
	SUBROUTINE checkpos(n,icheck)
C
C
C *********************************************************************
C Checks if we're close to the ending phase-space position...
C note - velocities offset from position....

	INCLUDE 'scfltr.h'
	INTEGER icheck,n
	REAL*8 dx,dy,dz,dd,vlos,dv,rframe,lastcheck,dcheck
	SAVE lastcheck

	IF(n.EQ.1)THEN
C
	   xend=xend/ru
	   yend=yend/ru
	   zend=zend/ru
	   vend=vend/vu
C
	   delx=delx/ru
	   delv=delv/vu
C
	   lastcheck=delx
	ENDIF

	rframe=SQRT(xframe*xframe+yframe*yframe+zframe*zframe)
C
	dx=(xframe-xend)
	dy=(yframe-yend)
	dz=(zframe-zend)
	dd=SQRT(dx*dx+dy*dy+dz*dz)
C
	vlos=(xframe*vxframe+yframe*vyframe+zframe*vzframe)/rframe
	dv=ABS(vlos-vend)

	icheck=0
	IF((dd.LT.delx).AND.(dv.LT.delv))THEN
	   dcheck=dd-lastcheck
	   IF(dcheck.GT.zero)icheck=1	
	   lastcheck=dd
	ENDIF

	WRITE(utermfil,99)icheck,n,dd,delx,dv,delv
 99	FORMAT(2I6,4(1pe10.2))

        RETURN
        END
C *********************************************************************
C
C
      SUBROUTINE findrem(n)
C
C
C *********************************************************************
      INCLUDE 'scfltr.h'

      INTEGER i,j,first,n
      REAL*8 msafe,
     &       ek(nbodsmax),
     &       ax0(nbodsmax),ay0(nbodsmax),az0(nbodsmax),
     &       vxrel,vyrel,vzrel,Lx,Ly,Lz,Lcore

      DATA first/0/
      SAVE first

C =====================================================================
      mrem=1.d10
      Lx=yframe*vzframe-zframe*vyframe
      Ly=-xframe*vzframe+zframe*vxframe
      Lz=xframe*vyframe-yframe*vxframe
      Lcore=SQRT(Lx*Lx+Ly*Ly+Lz*Lz)

C Save variables.....
      DO 5 i=1,nbodies
         ep(i)=pot(i)
         ax0(i)=ax(i)
         ay0(i)=ay(i)
         az0(i)=az(i)
         vxrel=vx(i)-vxframe
         vyrel=vy(i)-vyframe
         vzrel=vz(i)-vzframe
         ek(i)=.5d0*(vxrel*vxrel+vyrel*vyrel+vzrel*vzrel)
 5    CONTINUE

      DO 10 i=1,30
         msafe=mrem
         mrem=zero

C Find mass of remnant.....
         DO 20 j=1,nbodies
            IF(ek(j).GT.ABS(pot(j)))THEN
               ibound(j)=0
	       Lx=(y(j)+yframe)*vz(j)-(z(j)+zframe)*vy(j)
	       Ly=-(x(j)+xframe)*vz(j)+(z(j)+zframe)*vx(j)
	       Lz=(x(j)+xframe)*vy(j)-(y(j)+yframe)*vx(j)
	       IF(tub(j).EQ.zero)THEN
		  tub(j)=tnow
		  IF(SQRT(Lx*Lx+Ly*Ly+Lz*Lz).LE.Lcore)lmflag(j)=1.0e0
		  IF(SQRT(Lx*Lx+Ly*Ly+Lz*Lz).GT.Lcore)lmflag(j)=-1.0e0
	       ENDIF
            ELSE
               mrem=mrem+mass(j)
            END IF
 20      CONTINUE
C
         IF (msafe.LE.mrem) GO TO 25
         CALL accp_LH
C             -------
 10   CONTINUE

      mrem=zero
      
 25   	CALL outterm(' found rem in: ',i)
	IF((mrem.EQ.zero.AND.first.LT.2).
     &                               OR.tnow.LT.1.e-6) THEN
         first=first+1
         CALL outbods(n)
C             -------
         CALL iocoef
C             ------
      END IF

      eprem=zero
      ekrem=zero
      DO 30 i=1,nbodies
         IF(ibound(i).GT.0)THEN
            ekrem=ekrem+mass(i)*ek(i)
            eprem=eprem+0.5d0*mass(i)*pot(i)
         END IF
         ax(i)=ax0(i)
         ay(i)=ay0(i)
         az(i)=az0(i)
         pot(i)=ep(i)
         ibound(i)=1
 30   CONTINUE

      IF(mrem.EQ.zero)selfgrav=.FALSE.
      RETURN

      END
C *********************************************************************
C
C
      SUBROUTINE frame(n)
C
C
C *********************************************************
C
C Shifts phase space coordinates to be centred on minimum potential.
C
C-------------------------------------------------------------
      INCLUDE 'scfltr.h'
      INTEGER n,i,j,nend,rank(nbodsmax)
      REAL*8 xmin,ymin,zmin,mred

      SAVE rank
C ================================================================
C
C Find array that labels the position of each particle in a potential
C ranking.
        IF (MOD(n,nsort).EQ.0.OR.n.EQ.0) CALL indexx(nbodies,pot,rank)
C                                             ------
        nend=nbodies/100
C Find minimum of potential and shift origin.
        xmin=0.0d0
        ymin=0.0d0
        zmin=0.0d0
        mred=0.0d0

        DO 5 i=1,nend
           j=rank(i)
           xmin=xmin+mass(j)*x(j)
           ymin=ymin+mass(j)*y(j)
           zmin=zmin+mass(j)*z(j)
           mred=mred+mass(j)
 5      CONTINUE

        xmin=xmin/mred
        ymin=ymin/mred
        zmin=zmin/mred

C Update frame and shift to centre on minimum potential. 
        xframe=xframe+xmin
        yframe=yframe+ymin
        zframe=zframe+zmin

        DO 10 i=1,nbodies
           x(i)=x(i)-xmin
           y(i)=y(i)-ymin
           z(i)=z(i)-zmin
 10     CONTINUE

C Find v-frame
C        IF(MOD(n,noutlog).EQ.0)THEN
           vxframe=0.0d0
           vyframe=0.0d0
           vzframe=0.0d0
           DO 15 i=1,nend
              j=rank(i)
              vxframe=vxframe+mass(j)*vx(j)
              vyframe=vyframe+mass(j)*vy(j)
              vzframe=vzframe+mass(j)*vz(j)
 15        CONTINUE

           vxframe=vxframe/mred
           vyframe=vyframe/mred
           vzframe=vzframe/mred
C        ENDIF

        RETURN
        END

C***********************************************************************
C
C
                          SUBROUTINE stepsize
C
C
C***********************************************************************
C
C Subroutine to change stepsize
C
       INCLUDE 'scfltr.h'

       REAL*8 rfsq,dt1

       rfsq=xframe*xframe+yframe*yframe+zframe*zframe

C      dtime=dt0*rfsq/(rfsq+rch*rch)
	
       dt1=SQRT((rfsq+c*c)/tfac2/vh2)/1.d2

       dtime=MIN(dt0,dt1)

       RETURN
       END

C ********************************************************************
C
C
      SUBROUTINE tidalstart
C
C
C ********************************************************************
      INCLUDE 'scfltr.h'

      INTEGER i,j
      REAL*8 vxcm,vycm,vzcm,axcm,aycm,azcm,mtot,tee

      DO 10 i=1,ntide
C Advance position by one step.....
         CALL steppos
C             -------
C Find com velocity and acceleration.....
         vxcm=zero    
         vycm=zero
         vzcm=zero
         axcm=zero    
         aycm=zero
         azcm=zero
         mtot=zero

         DO 15 j=1,nbodies
            vxcm=vxcm+mass(j)*vx(j)
            vycm=vycm+mass(j)*vy(j)
            vzcm=vzcm+mass(j)*vz(j)
            axcm=axcm+mass(j)*ax(j)
            aycm=aycm+mass(j)*ay(j)
            azcm=azcm+mass(j)*az(j)
            mtot=mtot+mass(j)
 15      CONTINUE
C
         vxcm=vxcm/mtot      
         vycm=vycm/mtot      
         vzcm=vzcm/mtot      
         axcm=axcm/mtot      
         aycm=aycm/mtot      
         azcm=azcm/mtot 
C
         DO 20 j=1,nbodies
C Retard position by one step......
            x(j)=x(j)-vxcm*dtime 
            y(j)=y(j)-vycm*dtime
            z(j)=z(j)-vzcm*dtime
C Retard velocity by one step......
            vx(j)=vx(j)-axcm*dtime
            vy(j)=vy(j)-aycm*dtime
            vz(j)=vz(j)-azcm*dtime
 20      CONTINUE
C Increase tidal field......
         tee=REAL(i)/REAL(ntide)
         strength=(-2.*tee+3.)*tee*tee
C Find new accelerations......
	 IF(i.EQ.ntide)tnow=zero
         CALL accpot(i)
C             ------
         CALL stepvel
C             -------
	CALL outterm(' tidal start: ',i)
 10   CONTINUE

      tvel=0.5*dtime
      tpos=zero
      tnow=zero

      i=0
      CALL corrvel('correct')
C          -------
      CALL frame(i)
C          -----
      CALL findrem(i)
C          -------
      CALL outlog(i)
C          ------
      CALL corrvel('reset  ')
C          -------

      RETURN
      END
C *******************************************************************
C
C
      SUBROUTINE INDEXX(N,ARRIN,INDX)
C
C
C *******************************************************************
C
C Numerical Recipes subroutine to index arrays.
C
      INTEGER N
      INTEGER INDX(N),INDXT,I,IR,J,L
      REAL*8 Q,ARRIN(N)
C
C ==================================================================
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END

C***********************************************************************
C
C
      SUBROUTINE dynfric(n)
C
C       
C***********************************************************************
C Subroutine to calculate dynamical friction due to halo component alone
C
      INCLUDE 'scfltr.h'
      INTEGER n
      REAL*8 r,p,mr,aror,v,
     &       coulog,rho,sigma2,xx,df,erff,
     &       ot,xlim,jfit,mfricnow
      LOGICAL firstc
      DATA firstc/.TRUE./
      SAVE ot,sigma2,firstc,Gfric

      IF(firstc)THEN
         firstc=.FALSE.
         ot=1.0d0/3.0d0
         sigma2=vh2
      ENDIF

C Allow a linear mass-loss rate for the dynamical friction calculation
      mfricnow=mfric*(1-(tnow/tend)*0.8)
C      mfricnow=mfric

      r=SQRT(xframe*xframe+yframe*yframe+zframe*zframe)

C Calculate mass enclosed within r (=(1/r^2 G)(d phi /dr)
      mr=2.0d0*vh2*r*r*r/(r*r+c*c)/G

C Calculate speed - note: strictly should be updated so speed coincides with position
      v=SQRT(vxframe*vxframe+vyframe*vyframe+vzframe*vzframe)

C Coulomb logarithm via James' fitting function
      rtide=r*(ot*mfricnow/mr)**ot            
      xlim=rtide
      coulog=DLOG(r/rtide)+jfit(xlim)+0.135d0
      IF(coulog.GT.zero)THEN
C Local density (=(1/4 pi r^2)(d mr/dr))
         rho=(vh2/2.0d0/pi/G)*(r*r+3.0d0*c*c)/(r*r+c*c)/(r*r+c*c)

         xx=v/SQRT(two*sigma2)

C B+T, eqn 7-17
         df=4.0d0*pi*coulog*(G**2)*mfricnow*rho*
     &          (ERFF(xx)-two*xx*EXP(-(xx*2))/SQRT(pi))/(v**3)
      ELSE
         df=zero
      ENDIF

      dfx=-df*vxframe
      dfy=-df*vyframe
      dfz=-df*vzframe
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
