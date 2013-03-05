      subroutine xsbexriv(ear,ne,param,ifl,photar) 

      implicit   none
      integer    ne, ifl
      real       param(12),ear(0:ne),photar(ne) 

c 
c     Driver for angle-dependent reflection from an exponentially-cutoff
c     BROKEN power law and ionized medium.
c     It needs to be compiled together with xspexrav.f (neutral reflection).
c
c     See Magdziarz & Zdziarski, 1995, MNRAS.
c     See Zdziarski et al., 1995, ApJL January 10 for description of
c     calculation of ionization (based on Done et al. 1992).
c     The abundances are of Morrison & McCammon.
c
c     number of model parameters:11
c     1: Gamma1, first power law photon index
c     2: E_break, break energy
c     3: Gamma2, second power law photon index
c     4: E_c, the cutoff energy in keV (if E_c=0 there is no cutoff;
c        one needs to change the lower limit for that)
c     5: scale, scaling factor for reflection; if <0, no direct component
c        (scale=1 for isotropic source above disk)
c     6: redshift, z
c     7: abundance of elements heavier than He relative to
c        the solar abundances
c     8: iron abundance relative to the solar abundances
c     9: cosine of inclination angle
c     10: disk temperature in K
c     11: disk ionization parameter = L/nR^2
c     algorithm:
c          a[E*(1+z)]=E^{-Gamma}*exp(-E/E_c)+scale*reflection
c     Normalization is the photon flux at 1 keV (photons keV^-1 cm^-2 s^-1)
c     of the cutoff power law only (without reflection)
c     and in the earth frame.
c
      INCLUDE 'xspec.inc'

      INTEGER ixtot, isptot, ispref

      integer    n, ne_ , nesave, ierr
      real       xn,normfac,sppbk

      Real Gamma1, Gamma2, Eb, Ec, O2p, z, Ab, AbFe, Xm, Tdisk, xi
      Logical FirstCall
      character contxt*72

      save       FirstCall
      save       ixtot, isptot, ispref, nesave

      DATA firstcall/.TRUE./
      data ixtot, isptot, ispref, nesave/4*-1/

      ierr = 0

      IF ( firstcall ) THEN
       CALL xwrite('Compton reflection from ionized medium.',5)
       CALL xwrite('See help for details.',5)
       CALL xwrite('If you use results of this model in a paper,',5)
       CALL xwrite
     & ('please refer to Magdziarz & Zdziarski 1995 MNRAS, 273, 837'
     & ,5)
       CALL xwrite('Questions: Andrzej Zdziarski, aaz@camk.edu.pl',5)
       firstcall = .FALSE.
      ENDIF

      Gamma1 = param(1)
      Eb     = param(2)/511.
      Gamma2 = param(3)
      Ec     = param(4)/511.
      O2p    = param(5)
      z      = param(6)
      Ab     = alog10(param(7))
      AbFe   = alog10(param(8))
      Xm     = param(9)
      Tdisk  = param(10)
      xi     = param(11)

c Grab the memory for the arrays

      IF ( Ne .NE. nesave ) THEN
         CALL udmget(Ne+1, 6, ixtot, ierr)
         contxt = 'Failed to get memory for xtot'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL udmget(Ne+1, 6, isptot, ierr)
         contxt = 'Failed to get memory for sptot'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL udmget(Ne+1, 6, ispref, ierr)
         contxt = 'Failed to get memory for spref'
         IF ( ierr .NE. 0 ) GOTO 999
         nesave = ne
      ENDIF

      xn = (1+z)/511.0
      normfac = 1 / sppbk(1/xn, Gamma1, Gamma2, Eb, Ec)
      ne_ = ne + 1
      do n = 0, ne_-1 
         MEMR(ixtot+n) = (1 + z) * ear(n) / 511
         MEMR(isptot+n) = sppbk(1/MEMR(ixtot+n), Gamma1, Gamma2, Eb, Ec) 
      enddo
      call GrnTeeBi(MEMR(ixtot), ne_, MEMR(isptot), MEMR(ispref),
     &              Gamma1, Gamma2, Eb, Ec,
     &              O2p, Xm,
     &              Ab, AbFe,
     &              Tdisk, xi)

c  convert to photons and integrate over the bin size,
c  sptot is normalized at 1keV at the Earth frame
        do n=0,ne_-1
           MEMR(isptot+n)=normfac*MEMR(isptot+n)/ear(n)**2
        enddo
        do n=1,ne
           photar(n)=(MEMR(isptot+n)+MEMR(isptot+n-1))
     &               *(ear(n)-ear(n-1))/2
        enddo

 999  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         CALL xwrite(contxt, 10)
      ENDIF

        return
        end

c ------------------------------------------------------------------- c
 
      subroutine GrnTeeBi(xtot,jtot,sptot,spref,
     &                    Gamma1, Gamma2, Eb, Ec,
     &                    scale, xm,
     &                    ab_, ab0,
     &                    temp, xi)
c     calculates the reflected and total spectrum from pair cascades
c     using Green's functions of Magdziarz & Zdziarski 1995
      implicit   none
      integer    jtot
      Real       Gamma1, Gamma2, Eb, Ec
      Real       xtot(jtot),sptot(jtot),scale,
     &           xm,ab_,ab0,temp,xi
      integer    j
      real       spref(jtot)  

      if(scale.gt.0.)then
         call GrnCepBi(Gamma1, Gamma2, Eb, Ec, xtot,jtot,spref,
     &                     xm,ab_,ab0,temp,xi)
         do j=1,jtot
            sptot(j)=sptot(j)+spref(j)*scale
         enddo 
      elseif(scale.lt.0.)then
         call GrnCepBi(Gamma1, Gamma2, Eb, Ec, xtot,jtot,spref,
     &                     xm,ab_,ab0,temp,xi)
         do j=1,jtot
            sptot(j)=spref(j)*abs(scale)
         enddo
      endif
 
      return
      end
 
c ------------------------------------------------------------------ c
 
      subroutine GrnCepBi(Gamma1, Gamma2, Eb, Ec, xtot,jtot,spref,
     &                              xm,ab_,ab0,temp,xi)
c     calculates the reflected specutrum according to MZ95
c     ab0 - iron log10 abundance relative to hydrogen
      implicit  none
      integer   jtot
      real      gnrj 
      real      xtot(jtot),spref(jtot)
      real      Gamma1,Gamma2,Eb,Ec,xm,ab_,ab0,temp,xi
      integer   i,j 
      real      xmin,xtransl,xtransh,xrefmax,ymin
      real      pm1(19),pmx(26), apf, ap(3),xmo
      real      xkabs,xnor,xjl,xjd,fjc
      real      y,dym,dy,ddy,y0,sr,srp
      integer   nc,ncp,jmin,jtransl,jtransh,jrefmax
      real      SigmaBK,sppbk,grxyl,grxyh
 
c  precision factor for Green's function integration
      ncp=100
c  ranges for different methods
      xmin=2.e-4
      xtransl=0.01957
      xtransh=0.02348
      ymin=.03333
c  angle dependent parameters
      xmo=max(.05,min(.95,xm))
      call pm1y(xmo,pm1)
      call pmxy(xmo,pmx)
      ap(2)=apf(xmo)
      ap(3)=0.381
      xrefmax=1/(ymin+pmx(26))
 
      jmin=0
      do 50 j=1,jtot
         if(xtot(j).gt.xmin) goto 51
         jmin=j
 50   continue
 51   jtransl=0
      do 52 j=max(1,jmin),jtot
         if(xtot(j).gt.xtransl) goto 53
         jtransl=j
 52   continue
 53   jtransh=0
      do 54 j=max(1,jtransl),jtot
         if(xtot(j).gt.xtransh) goto 55
         jtransh=j
 54   continue
 55   jrefmax=0
      do 56 j=max(1,jtransh),jtot
         if(xtot(j).gt.xrefmax) goto 57
         jrefmax=j
 56   continue
 57   continue
 
      do 10 j=1,jmin
         spref(j)=0
 10   continue
c-----------------------
      do 20 j=jmin+1,jtransl
         xkabs=SigmaBK(j,xtot,jtransh,
     &           Gamma1,temp,xi,ab_,ab0,xnor)
         spref(j)=sppbk(1/xtot(j), Gamma1, Gamma2, Eb, Ec)
     &            *gnrj(xkabs,xmo)
 20   continue
c-----------------------
      if(jmin.ge.jtransl) then
         xkabs=SigmaBK(1,xtot,max(1,jtransh),
     &                  Gamma1,temp,xi,ab_,ab0,xnor)
      endif
      ap(1)=xnor/1.21
      xjl=alog(xtransl)
      xjd=alog(xtransh)-xjl
      do 60 j=jtransl+1,jtransh
         fjc=.5*sin(3.14159*((alog(xtot(j))-xjl)/xjd-.5))
         xkabs=SigmaBK(j,xtot,jtransh,
     &           Gamma1,temp,xi,ab_,ab0,xnor)
         y=1/xtot(j)
         spref(j)=sppbk(y, Gamma1, Gamma2, Eb, Ec)*gnrj(xkabs,xmo)
         dym=y-ymin
         dy=min(2.,dym)
         ddy=(dy-pmx(26))/(ncp+1)
         y0=y-dy
         sr=sppbk(y0, Gamma1, Gamma2, Eb, Ec)*grxyl(pm1,pmx,dy,y0,ap)
         dy=pmx(26)
         y0=y-dy
         sr=.5*(sr+sppbk(y0, Gamma1, Gamma2, Eb, Ec)
     &          * grxyl(pm1,pmx,dy,y0,ap))
         do 62 i=1,ncp
            dy=dy+ddy
            y0=y-dy
            sr=sr+sppbk(y0, Gamma1, Gamma2, Eb, Ec)
     &            *grxyl(pm1,pmx,dy,y0,ap)
 62      continue
         sr=sr*ddy
         if(dym.gt.2.)then
            ddy=dym-2
            nc=int(ncp*ddy/(dym-pmx(26)))
            ddy=ddy/(nc+1)
            dy=dym
            y0=y-dy
            srp=sppbk(y0,Gamma1,Gamma2,Eb,Ec)*grxyh(pm1,pmx,dy,y0,ap)
            dy=2
            y0=y-dy
            srp=.5*(srp+sppbk(y0, Gamma1, Gamma2, Eb, Ec)
     &            *grxyh(pm1,pmx,dy,y0,ap))
            do 64 i=1,nc
               dy=dy+ddy
               y0=y-dy
               srp=srp+
     &             sppbk(y0, Gamma1, Gamma2, Eb, Ec)
     &              *grxyh(pm1,pmx,dy,y0,ap)
 64         continue
            sr=sr+srp*ddy
         endif
         spref(j)=(.5-fjc)*spref(j)+(.5+fjc)*sr*xtot(j)
 60   continue
c-----------------------
      do 30 j=jtransh+1,jrefmax
         y=1/xtot(j)
         dym=y-ymin
         dy=min(2.,dym)
         ddy=(dy-pmx(26))/(ncp+1)
         y0=y-dy
         sr=sppbk(y0, Gamma1, Gamma2, Eb, Ec)*grxyl(pm1,pmx,dy,y0,ap)
         dy=pmx(26)
         y0=y-dy
         sr=.5*(sr+sppbk(y0, Gamma1, Gamma2, Eb, Ec)
     &        *grxyl(pm1,pmx,dy,y0,ap))
         do 32 i=1,ncp
            dy=dy+ddy
            y0=y-dy
            sr=sr+sppbk(y0, Gamma1, Gamma2, Eb, Ec)
     &           *grxyl(pm1,pmx,dy,y0,ap)
 32      continue
         sr=sr*ddy
         if(dym.gt.2.)then
            ddy=dym-2
            nc=int(ncp*ddy/(dym-pmx(26)))
            ddy=ddy/(nc+1)
            dy=dym
            y0=y-dy
            srp=sppbk(y0,Gamma1,Gamma2,Eb,Ec)*grxyh(pm1,pmx,dy,y0,ap)
            dy=2
            y0=y-dy
            srp=.5*(srp+sppbk(y0, Gamma1, Gamma2, Eb, Ec)
     &            *grxyh(pm1,pmx,dy,y0,ap))
            do 34 i=1,nc
               dy=dy+ddy
               y0=y-dy
               srp=srp+
     &             sppbk(y0, Gamma1, Gamma2, Eb, Ec)
     &             *grxyh(pm1,pmx,dy,y0,ap)
 34         continue
            sr=sr+srp*ddy
         endif
         spref(j)=sr*xtot(j)
 30   continue
c----------------------
      do 40 j=jrefmax+1,jtot
        spref(j)=0
 40   continue
 
      return
      end
 
c ------------------------------------------------------------------ c
      FUNCTION SIGMABK(In,X,N,Gamma0,Temp0,Xi0,Ab_0,Ab0,Asc0)
c     bound-free cross section for ionized matter with cosmic abundances
c     and variable iron abundance (in units of solar iron abundance)
c     per hydrogen atom in units of the Thomson cross section
c     e is in units of 511 keV

      INCLUDE 'xspec.inc'

      REAL ab , Ab0 , ab_ , Ab_0 , asc , Asc0 , gamma , Gamma0 , 
     &     SIGMABK , temp , Temp0 , WARMABK , xas , xi , Xi0
      INTEGER i , In , N , nold, ierr
      INTEGER iopac, ixold
      REAL X(*)
      LOGICAL fl1 , firstflag
      CHARACTER contxt*80

      CHARACTER*4 fgsolr,fgsolr0
      EXTERNAL fgsolr

      SAVE iopac , gamma , temp , xi , ab_ , ab , asc , ixold , nold , 
     &   fl1 , firstflag,fgsolr0

      DATA nold, iopac, ixold /0,2*-1/
      DATA gamma , temp , xi , ab_ , ab/0. , 0. , 0. , 99999. , 99999./
      DATA firstflag , fl1/2*.TRUE./
      DATA fgsolr0/'    '/

      ierr = 0

c In this version, changes of ear during one xspec session are allowed.

      if(fgsolr().ne.fgsolr0)then
         fl1=.true.
         fgsolr0=fgsolr()
      endif

      IF ( (gamma.NE.Gamma0) .OR. (temp.NE.Temp0) .OR. (xi.NE.Xi0) .OR. 
     &     (ab_.NE.Ab_0) .OR. (ab.NE.Ab0) ) THEN
         gamma = Gamma0
         temp = Temp0
         xi = Xi0
         ab = Ab0
         ab_ = Ab_0
         firstflag = .TRUE.
         fl1 = .TRUE.
         xas = .01947
         asc = 1.503E2*WARMABK(xas,gamma,temp,xi,ab_,ab,firstflag,
     &                         .false.)*xas**3
      ENDIF
      IF ( (X(In).NE.MEMR(ixold+In-1)) .OR. (N.NE.nold) ) THEN
         CALL udmget(N, 6, iopac, ierr)
         contxt = 'Failed to get memory for opac'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL udmget(N, 6, ixold, ierr)
         contxt = 'Failed to get memory for xold'
         IF ( ierr .NE. 0 ) GOTO 999
         nold = N
         DO 50 i = 1 , N
            MEMR(ixold+i-1) = X(i)
 50      CONTINUE
         fl1 = .TRUE.
      ENDIF
 
      IF ( fl1 ) THEN
         DO 100 i = 1 , N
            MEMR(iopac+i-1) = 1.503E2*WARMABK(X(i),gamma,temp,xi,ab_,ab,
     &                firstflag,.false.)
 100     CONTINUE
         fl1 = .FALSE.
      ENDIF
 
      Asc0 = asc
      SIGMABK = MEMR(iopac+In-1)

 999  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         CALL xwrite(contxt, 10)
         SIGMABK = 0.
      ENDIF
 
      RETURN
      END


c -------------------------------------------------------------------

      FUNCTION WARMABK(X0,Gamma,Temp,Xi,Ab_,Ab0,Firstflag,qHHe)

      REAL WARMABK , X0 , Gamma , Temp , Xi , Ab_ , Ab0
      LOGICAL Firstflag, qHHe

      INCLUDE 'xspec.inc'

      REAL clcwrm3
      REAL aw(20), abo(20), e(800)
      INTEGER z(20), iab0, imax
      INTEGER isigma, iion, iint, iarec, inum, iratio, ierr
      CHARACTER*72 contxt

c      COMMON /aprmdt/ isigma, iion, z, aw, abo, e, iab0, imax

      SAVE isigma, iion, iint, iarec, inum, iratio
      SAVE z, aw, abo, e, iab0, imax

      DATA isigma, iion, iint, iarec, inum, iratio /6*-1/

      ierr = 0

C Get the memory for the arrays required by CLCWRM3

      IF ( isigma .EQ. -1 ) THEN
         CALL udmget(20*30*800, 6, isigma, ierr)
         contxt = 'Failed to get memory for sigma array'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL xwrite('Memory acquired for sigma array', 25)
      ENDIF
      IF ( iion .EQ. -1 ) THEN
         CALL udmget(20*30*10, 6, iion, ierr)
         contxt = 'Failed to get memory for ion array'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL xwrite('Memory acquired for ion array', 25)
      ENDIF
      IF ( iint .EQ. -1 ) THEN
         CALL udmget(20*30, 6, iint, ierr)
         contxt = 'Failed to get memory for int array'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL xwrite('Memory acquired for int array', 25)
      ENDIF
      IF ( iarec .EQ. -1 ) THEN
         CALL udmget(20*30, 6, iarec, ierr)
         contxt = 'Failed to get memory for arec array'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL xwrite('Memory acquired for arec array', 25)
      ENDIF
      IF ( inum .EQ. -1 ) THEN
         CALL udmget(20*30, 7, inum, ierr)
         contxt = 'Failed to get memory for num array'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL xwrite('Memory acquired for num array', 25)
      ENDIF
      IF ( iratio .EQ. -1 ) THEN
         CALL udmget(20*30, 7, iratio, ierr)
         contxt = 'Failed to get memory for ratio array'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL xwrite('Memory acquired for ratio array', 25)
      ENDIF

C Call the routine to do the calculation

      WARMABK = CLCWRM3(X0,Gamma,Temp,Xi,Ab_,Ab0,Firstflag,qHHe,
     &                 MEMR(isigma),MEMR(iion),MEMR(iint),MEMR(iarec),
     &                 MEMD(inum),MEMD(iratio), z, aw, abo, e, iab0,
     &                 imax)

 999  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         CALL xwrite(contxt, 10)
         WARMABK = 0.
      ENDIF

      RETURN
      END


c -------------------------------------------------------------------
 
      FUNCTION CLCWRM3(X0,Gamma,Temp,Xi,Ab_,Ab0,Firstflag,qHHe,sigma,
     &                ion,int,arec,num,ratio,z,aw,abo,e,iab0,imax)

c 	calculate cross-sections using Daltabuit and Cox, use cosmic
c       abundances, and radiative recombination rates from
c       Aldrovandi and Pequignot 1973 and explicit formula
c       for recombination to hydrogenlike atoms
c       result in units of 10^-22
      IMPLICIT NONE
      REAL CLCWRM3 , X0 , Gamma , Temp , Xi , Ab_ , Ab0
      REAL t4 , z2 , e1 , e2 , re , rp , ab_d
      REAL aw(20) , abo(20) , ab(20) , ion(20,30,*) , y
      REAL int(20,*) , arec(20,*)
      REAL denom , xp , sigmaeff , e(800) , sigma(20,30,*)
      REAL fgabnd
      REAL*8 num(20,*) , sum , ratio(20,*) , mult , mul(30) , xil
      INTEGER z(20) , i , j , imax , n , ios , k , iab0
      CHARACTER*4 fgsolr0
      LOGICAL Firstflag, qHHe

      CHARACTER fgsolr*4, el_name*2
      EXTERNAL fgsolr, el_name
c
      DATA fgsolr0/'    '/
      SAVE ab,fgsolr0

C If necessary load the arrays from the iond_mm.dat and mansig.dat
C files

      CALL ldapmr2(z, aw, abo, ion, sigma, e, iab0, imax)

C initialize other arrays

      if((fgsolr().ne.fgsolr0).or.firstflag)then
         fgsolr0=fgsolr()
         do i=1,imax
            abo(i)=fgabnd(el_name(z(i)))
         enddo
c  insert external abundances (ab_, ab0)
         ab_d = 10.**Ab_
         DO 100 i = 1 , imax
            IF ( z(i).GT.2 ) THEN
               ab(i) = abo(i)*ab_d
            ELSE
               IF ( qHHe ) THEN
                  ab(i) = abo(i)
               ELSE
                  ab(i) = 0.
               ENDIF
            ENDIF
 100     CONTINUE
         IF ( iab0.GT.0 ) ab(iab0) = abo(iab0)*10.**Ab0
      endif

      IF ( Firstflag ) THEN
         IF ( Xi.GT.1.E-20 ) THEN
c  normalized temperature
            t4 = Temp*.0001
c  calculate recombination coefficients - radiative and dielectric
c  from Aldrovandi+Pequignot 1973 Astro, astrophys 25 137
            DO 120 i = 1 , imax , 1
               DO 110 j = 1 , z(i) - 1 , 1
                  e1 = EXP(-ion(i,j,5)/t4)
                  e2 = EXP(-ion(i,j,7)/t4)
                  arec(i,j) = ion(i,j,2)*(t4**(-ion(i,j,3)))
     &                        + ion(i,j,4)*(t4**(-1.5))
     &                        *e1*(1.+ion(i,j,6)*e2)
 110           CONTINUE
c  do radiative recomb to hydrogenic ion NB assumed gaunt factor=1
c  from Gould+Thakur 1970 ann phys 61 351.
               z2 = z(i)**2
               y = 15.8*z2/t4
               arec(i,z(i)) = 1.033E-3*z2*(1.735+LOG(y)+1/(6.*y))
     &                        /SQRT(t4)
 120        CONTINUE
c  first normalise photon spectrum
            denom = 0.
            DO 140 k = 1 , 720 , 1
               denom = denom + (e(k)**(1-Gamma))*(e(k+1)-e(k))
 140        CONTINUE
            denom = 1/denom
c  calculate the photoionization integrals and population balance
            DO 160 i = 1 , imax , 1
               DO 150 j = 1 , z(i) , 1
                  int(i,j) = 0.0
c  do integral in log10 (eV) from lowest energy
c  to the maximum energy- use same energy grid
c  as the cross-sections are defined on
                  DO 145 k = 1 , 720 , 1
                     int(i,j) = int(i,j) + sigma(i,j,k)*(e(k)**(-Gamma))
     &                          *(e(k+1)-e(k))
 145              CONTINUE
                  int(i,j) = int(i,j)*denom
c  const=7.958d-2/2.43d+4
                  ratio(i,j) = DLOG(3.2749D-6*int(i,j)/arec(i,j))
 150           CONTINUE
 160        CONTINUE
c  calculate population numbers
c  set sum=1 as the sum=1+A21+A32A21+....
c  complicated because floating range limit
            DO 220 i = 1 , imax , 1
               xil = DLOG(DBLE(Xi))
               mult = 0.D0
               DO 170 j = z(i) , 1 , -1
                  mul(j) = 0.D0
                  DO 165 n = j , 1 , -1
                     mul(j) = mul(j) + ratio(i,n)
 165              CONTINUE
                  mul(j) = mul(j) + j*xil
                  IF ( mul(j).GT.mult ) mult = mul(j)
 170           CONTINUE
               DO 180 j = z(i) , 1 , -1
                  mul(j) = mul(j) - mult
 180           CONTINUE
               sum = 0.D0
               DO 190 j = z(i) , 1 , -1
                  sum = sum + DEXP(mul(j))
 190           CONTINUE
               sum = sum + DEXP(-mult)
               num(i,1) = -mult - DLOG(sum)
               DO 200 j = 2 , z(i) + 1 , 1
                  num(i,j) = num(i,j-1) + ratio(i,j-1) + xil
 200           CONTINUE
               DO 210 j = 1 , z(i) + 1 , 1
                  num(i,j) = DEXP(num(i,j))
 210           CONTINUE
 220        CONTINUE
         ELSE
            DO 240 i = 1 , imax , 1
               num(i,1) = 1.
               DO 230 j = 2 , z(i) + 1 , 1
                  num(i,j) = 0.
 230           CONTINUE
 240        CONTINUE
         ENDIF
         Firstflag = .FALSE.
      ENDIF
c ------------------------------------------------------------------- c
 
      xp = X0*5.11E5
      sigmaeff = 0.
      IF ( xp.LT.e(721) ) THEN
c  pick closest energy bin
         k = 1
         DO WHILE ( e(k).LT.xp )
            k = k + 1
         ENDDO
         re = (xp-e(k-1))/(e(k)-e(k-1))
         DO 300 i = 1 , imax , 1
            DO 260 j = 1 , z(i) , 1
               sigmaeff = sigmaeff + SNGL(num(i,j))*ab(i)
     &                    *(sigma(i,j,k-1)
     &                    +re*(sigma(i,j,k)-sigma(i,j,k-1)))
 260        CONTINUE
 300     CONTINUE
      ELSE
         re = (xp/e(721))**(-3)
         DO 350 i = 1 , imax , 1
            rp = ab(i)*re
            DO 320 j = 1 , z(i) , 1
               sigmaeff = sigmaeff + SNGL(num(i,j))*sigma(i,j,721)*rp
 320        CONTINUE
 350     CONTINUE
      ENDIF
      sigmaeff = sigmaeff*6.6E-5
 
      CLCWRM3 = sigmaeff

      RETURN
      END
 
c ------------------------------------------------------------------- c
 
      SUBROUTINE ldapmr2(z, aw, abo, ion, sigma, e, iab0, imax)

      REAL aw(20), abo(20), ion(20,30,*), sigma(20,30,*), e(*)
      INTEGER z(20), iab0

      INTEGER i, j, k, imax, ios
      INTEGER ilun10 , ilun11 , lenn
      REAL fgabnd
      CHARACTER*2 el_name
      CHARACTER*128 absfil, wrtstr
      LOGICAL qfirst

      INTEGER lenact
      CHARACTER fgdatd*128
      EXTERNAL lenact, fgdatd

      DATA qfirst /.TRUE./

      IF ( .NOT.qfirst ) RETURN

      qfirst = .FALSE.

c  set up constants and read in data if first time round
 
      iab0 = 0
      CALL GETLUN(ilun10)
      absfil = fgdatd()
      lenn = lenact(absfil)
      absfil = absfil(:lenn)//'iond_mm.dat'
      CALL OPENWR(ilun10,absfil,'old',' ',' ',0,0,ios)
 
c          open(unit=10,file='iond_mm.dat',status='old',iostat=ios)
c          open(unit=11,file='mansig.data',status='old')
c          ilun10=10
c          ilun11=11
 
      i = 1
      READ (ilun10,*,IOSTAT=ios) z(i) , aw(i) , abo(i)
      IF ( z(i).EQ.26 ) iab0 = i
c  NB j=1 is ground state, j=z(i)=fully charged
      DO WHILE ( ios.EQ.0 )
         DO 20 j = 1 , z(i) , 1
            READ (ilun10,*) ion(i,j,1) , ion(i,j,2) , ion(i,j,3) , 
     &                      ion(i,j,4) , ion(i,j,5) , ion(i,j,6) , 
     &                      ion(i,j,7) , ion(i,j,8) , ion(i,j,9) , 
     &                      ion(i,j,10)
c  change units of coefficients
            ion(i,j,2) = ion(i,j,2)*1.E+10
            ion(i,j,4) = ion(i,j,4)*1.E+04
            ion(i,j,5) = ion(i,j,5)*1.E-04
            ion(i,j,7) = ion(i,j,7)*1.E-04
 20      CONTINUE
         i = i + 1
         READ (ilun10,*,IOSTAT=ios) z(i) , aw(i) , abo(i)
         IF ( z(i).EQ.26 ) iab0 = i
      ENDDO
      imax = i - 1

      do i=1,imax
         abo(i)=fgabnd(el_name(z(i)))
      enddo
 
c       ******new bit for xspec
      CLOSE (UNIT=ilun10)
      CALL FRELUN(ilun10)
 
c  check not overflowing array bounds
      IF ( imax.GT.20 .OR. z(imax).GT.30 ) THEN
         CALL xwrite('element array overflowing array bounds', 5)
         WRITE(wrtstr, '(a,i6)') 'number of elements= ' , imax
         CALL xwrite(wrtstr, 5)
         WRITE(wrtstr, '(a,i6)') 'maximum atomic number= ' , z(imax)
         CALL xwrite(wrtstr, 5)
         RETURN
      ENDIF
c  read in photoionization data from Reilman and Manson ApJ supp 40 1979
      CALL GETLUN(ilun11)
      absfil = fgdatd()
      lenn = lenact(absfil)
      absfil = absfil(:lenn)//'mansig.dat'
      CALL OPENWR(ilun11,absfil,'old',' ',' ',0,0,ios)

      DO 50 i = 1 , imax , 1
c  now go over all ion states - values come from Reilman
c  ApJ supp 1979 except hydrogenic which is analytic
         DO 40 j = 1 , z(i) , 1
c  ion identification header
            READ (ilun11,*)
c  read in over all energies
c  the table is valid for 5.eV-19907.7eV
            DO 30 k = 1 , 721 , 1
               READ (ilun11,*) e(k) , sigma(i,j,k)
c  change units of sigma
               sigma(i,j,k) = sigma(i,j,k)/6.6E-27
 30         CONTINUE
 40      CONTINUE
 50   CONTINUE
c       ******new bit for xspec
      CLOSE (UNIT=ilun11)
      CALL FRELUN(ilun11)
 
c           close(ilun10)
c           close(ilun11)
 
c  end of data set up
c ------------------------------------------------------------------- c

      RETURN
      END










c------------------------------------------------------------------- c

      Real function Sppbk(y, Gamma1, Gamma2, Eb, Ec)
c
c Broken power-law with exponential cutoff.
c
      Implicit None
      Real y, Gamma1, Gamma2, Eb, Ec
      Real A, S

c  y = 1/x
      
      A = Eb**(Gamma2 - Gamma1)
      if (1/y.gt.Eb) then
        S = A * y**(Gamma2 - 2)
      else
        S = y**(Gamma1 - 2)
      endif
      if (Ec.eq.0.0) then
        Sppbk = S
      else
        Sppbk = S * Exp(-1 / (Ec * y))
      endif
      return
      end
