C***********************************************************************
C
C
                             PROGRAM maketrail
C
C
C***********************************************************************
C See K.V. Johnston, 1998 ApJ, 495, 297 for details
C (equation numbers refer to this paper)
C
C Program to find density, distance and velocity distribution
C from heliocentric viewpoint along tidal streamers from
C Galactic satellite.
C Assumes that an equal fractional mass instantaneously on 
C each pericentric passage.
C Requires subroutines 'orbit.f' to integrate orbit
C
C Input files 'param' contains
C     msat      - current mass of the satellite in solar masses
C     f         - fractional mass loss
C     mfric     - Mass of the satellite which acts for dynamical friction
C     tnow      - time in Gyrs over which this mass is lost
C                 (current time is labelled as tnow)
C     vc        - circular velocity of logarithmic halo that best
C     x0        - current Galactocentric position of satellite in kpc
C     v0        - current Galactocentric velocity of satellite in km/s
C                 approximates chosen Galactic model
C
C Outputs:
C 'orb.dat' contains (columns)
C     t         - time in Gyrs
C     x         - Galactocentric position in kpc
C     v         - Galactocentric velocity in km/s
C     psi       - angle along orbit in radians
C     pdot      - angular velocity in radians/Gyr
C     e         - energy (checks conservation) 
C
C 'encounters.dat'
C line 1:
C nencounters, msat, msat0, tnow
C lines 1-nencounters: 
C pericenter, instantaneous msat, mass lost, rtide/peri
C
C 'npsi.dat'
C line 1: 
C nencounters
C For encounters 1-nencounter:
C nlead, ntrail - number of points inleading trailing streamer
C                 from the encounter
C next nlead+ntrail lines (in columns):
C     psi       - angular position in radians from satellite along orbit
C     dmdp      - fraction of total debris/radian at this point
C     de        - offset in orbital energy from satellite of
C                 this material
C     (l,b)     - Heliocentric longitude and latitude in degrees
C     r         - Heliocentric distance in kpc
C     v         - line of sight velocity from Sun in Galactocentric 
C                 rest-frame
C =================================================================
      INCLUDE 'maketrail.h'
      INTEGER n

C Input parameters and define what's needed from these
      CALL get_param

C          ---------
C Integrate backwards and forwards from current time to
C find position and velocity at each point along the orbit,
C pericentric distances and number of encounters
      CALL orbit
C          -----
C Define quantities at encounters
      CALL encounter
C          ---------
C For each encounter, find  mass surface density
C distribution of debris today 
      
      DO 10 n=1,nenc
C
         CALL debris(n)
C             ------
 10   CONTINUE

C Output for trail produced at each encounter
C from Heliocentric viewpoint:
C     mass surface density
C     (l,b) pair along trail
C     average distance
C     average radial velocity
C
      CALL output
C          ------

      STOP
      END



C***********************************************************************
C
C
      SUBROUTINE debris(n)
C
C
C***********************************************************************
      INCLUDE 'maketrail.h'
      INTEGER ipsi,n
      REAL*8 eps,w,u,tau,q,dmdp,dndq

      nlead(n)=0
      ntrail(n)=0
C equation (8)
      eps=vc4o3*gmoro3(n)
C coefficient in equation (15)
      w=vc2o3/gmoro3(n)
      DO 10 ipsi=1,nout
         u=tout(ipsi)-tenc(n)
         IF(u.GT.0.)THEN
C equation (11)
            tau=(tnow-tenc(n))/u
            de(n,ipsi)=vc2*DLOG(tau)
            q=ABS(de(n,ipsi)/eps)
            IF(dndq(q).GT.0.)THEN
               IF(de(n,ipsi).LT.0)THEN
                  ntrail(n)=ntrail(n)+1
               ELSE
                  nlead(n)=nlead(n)+1
               END IF
            END IF
C equation (15)
            dmdp(n,ipsi)=w*dndq(q)/u/pdot(ipsi)
C            
            dmdp(n,ipsi)=dmdp(n,ipsi)*dm(n)/(msat0-msat)
C
         ELSE
            dmdp(n,ipsi)=0.
            de(n,ipsi)=0.
         ENDIF
 10   CONTINUE

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION dndq(q)
C
C
C***********************************************************************
      REAL*8 q,dndq,a,b,c,h
      PARAMETER(a=0.6,b=1.1,c=3.0,h=1./(c-a))

C equation (13)
      IF(q.GT.c)THEN
         dndq=0.
      ELSEIF(q.GT.b)THEN
         dndq=(1.-(q-b)/(c-b))
      ELSEIF(q.GT.a)THEN
         dndq=-(q-a)/(a-b)
      ELSE
         dndq=0.
      ENDIF

      dndq=dndq*h

      RETURN
      END
 

C***********************************************************************
C
C
      SUBROUTINE encounter
C
C
C***********************************************************************
      INCLUDE 'maketrail.h'
      INTEGER n
c      REAL*8 f,G,width
      REAL*8 G,width
      PARAMETER(G=4.296e-6)

c      f=1.-(msat/msat0)**(1./REAL(nenc))
      msat0=msat/(1.-f)**REAL(nenc)
      OPEN(UNIT=10,FILE='encounters.dat',STATUS='UNKNOWN')
C      WRITE(10,99)nenc,msat,msat0,tnow
      WRITE(10,99)nenc,msat,msat0,f,tnow

      mass(1)=msat0
      DO 20 n=1,nenc
         dm(n)=f*mass(n)
         mass(n+1)=mass(n)-dm(n)
         gmoro3(n)=(G*mass(n)/peri(n))**(1./3.)
         width=gmoro3(n)/vc2o3
         WRITE(10,100)tenc(n),peri(n),mass(n),dm(n),width
 20   CONTINUE
      CLOSE(10)

 99   FORMAT(1I5,4(1pe10.2))
 100  FORMAT(5(1pe12.2))

      RETURN
      END

C***********************************************************************
C
C
      SUBROUTINE get_param
C
C
C***********************************************************************
      INCLUDE 'maketrail.h'
      INTEGER i
      REAL*8 x,t,dt,vc,tu,c,vh2,dfx,dfy,dfz
      COMMON/reals/x(6),t,dt,tu,c,vh2,dfx,dfy,dfz

      OPEN(UNIT=10,FILE='param_auto',STATUS='OLD')
      READ(10,*)msat
      READ(10,*)drlphi
      READ(10,*)q1,q2,q3
      READ(10,*)dhalo
      READ(10,*)vhalo
      READ(10,*)f
      READ(10,*)mfric
      READ(10,*)tnow
      READ(10,*)vc
      READ(10,*)(x(i),i=1,3)
      READ(10,*)(x(i),i=4,6)
      CLOSE(10)

      vc2=vc*vc
      vc4o3=vc2**(2./3.)
      vc2o3=SQRT(vc4o3)
C

      RETURN
      END

C***********************************************************************
C
C
      SUBROUTINE output
C
C
C***********************************************************************
      INCLUDE 'maketrail.h'
      INTEGER ipsi,n
      REAL*8 dr,l,b,ex,why,zed,r,vr,pi,deg,rsun
      PARAMETER(rsun=8.0)

      pi=4.*ATAN(1.0)
      deg=180./pi

      OPEN(UNIT=10,FILE='npsi.dat',STATUS='UNKNOWN')
      WRITE(10,98)nenc
      DO 5 n=1,nenc
         WRITE(10,98)nlead(n),ntrail(n)
         DO 10 ipsi=1,nout
            IF(dmdp(n,ipsi).GT.0.)THEN
C Trail offset in distance from orbit (as shown in LMC paper with
C Jamie, Steve and Bill)
               dr=rr(ipsi)*de(n,ipsi)/vc2
C Heliocentric viewpoint
               ex=(rr(ipsi)+dr)*y(1,ipsi)/rr(ipsi)
C+rsun
               why=(rr(ipsi)+dr)*y(2,ipsi)/rr(ipsi)
               zed=(rr(ipsi)+dr)*y(3,ipsi)/rr(ipsi)
               WRITE(10,99)psi(ipsi),dmdp(n,ipsi),de(n,ipsi),
     &               ex,why,zed,y(4,ipsi),y(5,ipsi),y(6,ipsi)

C               r=SQRT(ex*ex+why*why+zed*zed)
C               vr=(ex*y(4,ipsi)+why*y(5,ipsi)+zed*y(6,ipsi))/r
C               l=ATAN2(why,ex)*deg
C               IF(l.LT.0.)l=l+360.
C               b=ASIN(zed/r)*deg
C               WRITE(10,99)psi(ipsi),dmdp(n,ipsi),de(n,ipsi),
C     &                    l,b,r,vr
            ENDIF
 10      CONTINUE
 5    CONTINUE

 98   FORMAT(2I10)
 99   FORMAT(10(1pe12.4))
      CLOSE(10)

      RETURN
      END
