      INTEGER uterm,nenc,nmax,nperi,npsi,nout,ntrail,nlead
      PARAMETER(uterm=6,nmax=10000,nperi=20,npsi=400)
      REAL*8 tout,psi,pdot,tenc,peri,
     &       tnow,msat0,msat,vc4o3,vc2o3,gmoro3,dm,mass,vc2,
     &       y,e,x0,y0,de,rr,vrr,j2,f,mfric
      COMMON/ints/nenc,nout,ntrail(nperi),nlead(nperi)
      COMMON/param/tnow,msat0,msat,vc4o3,vc2o3,vc2,j2,f,mfric,drlphi,
     &       vhalo,dhalo,q1,q2,q3
      COMMON/orb/tout(nmax),psi(nmax),pdot(nmax),y(6,nmax),e(nmax),
     &           rr(nmax),vrr(nmax)
      COMMON/enc/tenc(nperi),gmoro3(nperi),mass(nperi),dm(nperi),
     &           peri(nperi)
      COMMON/den/dmdp(nperi,nmax),de(nperi,nmax)
      COMMON/proj/x0(3),y0(3)
