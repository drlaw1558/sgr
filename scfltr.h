C=======================================================================
C
C
C                        INCLUDE FILE scfltr.h
C
C
C=======================================================================
C
C
C     Parameter declarations, allocation of array storage, common
C     block definitions.
C
C
C=======================================================================

        INTEGER nbodsmax,nmax,lmax
 
        PARAMETER(nbodsmax=100000,nmax=6,lmax=4)

        CHARACTER*50 headline
        CHARACTER*4 outfile
        INTEGER nsteps,noutbod,nbodies,noutlog,ibound,ntide,nsort,
     &          nproc
        LOGICAL selfgrav,inptcoef,outpcoef,zeroodd,zeroeven,fixacc,
     &          david,dfric,disk,flat,growing,flattening
        REAL*8 tnow,x,y,z,vx,vy,vz,mass,pot,dtime,G,ax,ay,az,one,pi,
     &         twoopi,onesixth,tpos,tvel,
     &         potext,two,zero,tiny,
     &         mrem,ekrem,eprem,xframe,yframe,zframe,vxframe,
     &         vyframe,vzframe,sinsum,cossum,ep,strength,ru,mu,vu,dt0,
     &         tfac2,c,vh2,tend,tub,lmflag,dfx,dfy,dfz,mfric,rtide,
     &         xend,yend,zend,delx,vend,delv
	REAL cputime0,cputime1,cputime,tick
        COMMON/charcom/outfile
        COMMON/bodscom/x(nbodsmax),y(nbodsmax),z(nbodsmax),vx(nbodsmax),
     &                 vy(nbodsmax),vz(nbodsmax),mass(nbodsmax),
     &                 pot(nbodsmax),ax(nbodsmax),ay(nbodsmax),
     &                 az(nbodsmax),potext(nbodsmax),ep(nbodsmax),
     &                 tub(nbodsmax), lmflag(nbodsmax)
        COMMON/parcomi/nbodies,nsteps,noutbod,noutlog,ibound(nbodsmax),
     &                 ntide,nsort,nproc
        COMMON/parcomr/dtime,G,one,pi,twoopi,onesixth,two,tiny,zero
        COMMON/parcomc/headline
        COMMON/parcoml/selfgrav,inptcoef,outpcoef,zeroodd,zeroeven,
     &                 fixacc,david,dfric,disk,flat,growing,flattening
        COMMON/timecom/tpos,tnow,tvel,dt0,tfac2,tend
        COMMON/cpucom/cputime0,cputime1,cputime,tick
        COMMON/coefcom/sinsum(0:nmax,0:lmax,0:lmax),
     &                 cossum(0:nmax,0:lmax,0:lmax)
        COMMON/remcom/mrem,ekrem,eprem
        COMMON/orbcom/xframe,yframe,zframe,vxframe,vyframe,vzframe,
     &                strength,ru,mu,vu,c,vh2,dfx,dfy,dfz,mfric,rtide,
     &                drlphi,q1,q2,q3,dhalo,vhalo
        COMMON/checkcom/xend,yend,zend,delx,vend,delv
C=======================================================================
C   Definitions specific to input/output.
C=======================================================================
        INTEGER uterm,upars,ulog,ubodsin,ubodsout,utermfil,uoutcoef,
     &          uincoef,ubodsel,ubodsout1
        CHARACTER*8 parsfile,logfile,termfile,
     &              outcfile,incfile,elfile
        CHARACTER*8 ibodfile,obodfile
        PARAMETER(uterm=6,upars=10,ulog=11,ubodsin=12,ubodsout=13,
     &            utermfil=15,uoutcoef=16,uincoef=17,ubodsel=18,
     &		  ubodsout1=19)
        PARAMETER(parsfile='SCFPAR',logfile='LOG',
     &            ibodfile='BI',
     &            obodfile='SNAPxxxx',
     &            termfile='OUT',outcfile='COEF',
     &            incfile='SCFICOEF',elfile='SCFELxxx')






