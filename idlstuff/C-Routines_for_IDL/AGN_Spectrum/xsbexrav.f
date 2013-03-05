      subroutine xsbexrav(ear,ne,param,ifl,photar) 

      implicit   none
      integer    ne, ifl
      real       param(9),ear(0:ne),photar(ne) 


c 
c     Driver for angle-dependent reflection from an exponentially-cutoff
c     BROKEN power law and neutral medium.
c 
c     number of parameters: 9
c     1: Gamma1, first power law photon index
c     2: E_break, break energy
c     3: Gamma2, second power law photon index
c     4: E_c, the cutoff energy in keV (if E_c=0 there is no cutoff;
c        one needs to change the lower limit for that)
c     5: scaling factor for reflection (1 for isotropic source above disk) 
c     6: cosine of inclination angle 
c     7: abundance of elements heavier than He relative to
c        the solar abundances
c     8: iron abundance relative to the solar iron abundance (7.67)
c     9: redshift 
c 

      INCLUDE 'xspec.inc'

      INTEGER ixtot, isptot, ispref

      integer    n, ne_, nesave, ierr
      real       xn,normfac,SppBk   

      Real Gamma1, Gamma2, Eb, Ec, O2p, z, Ab, AbFe, Xm
      Logical FirstCall
      character contxt*72

      save       FirstCall
      save       ixtot, isptot, ispref, nesave

      data FirstCall/.True./
      data ixtot, isptot, ispref, nesave/4*-1/


      Gamma1 = param(1)
      Eb     = param(2)/511.
      Gamma2 = param(3)
      Ec     = param(4)/511.
      O2p    = param(5)
      Xm     = param(6)
      Ab     = alog10(param(7))
      AbFe   = alog10(param(8))
      z      = param(9)

      ierr = 0

      IF ( firstcall ) THEN
       CALL xwrite('Compton reflection from neutral medium.',5)
       CALL xwrite('See help for details.',5)
       CALL xwrite('If you use results of this model in a paper,',5)
       CALL xwrite
     & ('please refer to Magdziarz & Zdziarski 1995 MNRAS, 273, 837'
     & ,5)
       CALL xwrite('Questions: Andrzej Zdziarski, aaz@camk.edu.pl',5)
       firstcall = .FALSE.
      ENDIF

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
      call GrnTeeBk(MEMR(ixtot), ne_, MEMR(isptot), MEMR(ispref),
     &              Gamma1, Gamma2, Eb, Ec,
     &              O2p, Xm, Ab, AbFe)


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
c ------------------------------------------------------------------- c
 
      subroutine GrnTeeBk
     &   (xtot, jtot, sptot,spref,
     &    Gamma1, Gamma2, Eb, Ec,
     &    scale, xm, ab, ab0)
c     calculates the reflected and total spectrum from pair cascades
c     using Green's functions of Magdziarz & Zdziarski 1995
      implicit   none
      integer    jtot
      Real       Gamma1, Gamma2, Eb, Ec
      Real       xtot(jtot),sptot(jtot),scale,xm,ab,ab0
      integer    j
      real       spref(jtot)  

      if(scale.gt.0.)then
         call 
     &     GrnCepBk(Gamma1, Gamma2, Eb, Ec, xtot,jtot,spref,xm,ab,ab0)
         do j=1,jtot
            sptot(j)=sptot(j)+spref(j)*scale
         enddo 
      elseif(scale.lt.0.)then
         call 
     &     GrnCepBk(Gamma1, Gamma2, Eb, Ec, xtot,jtot,spref,xm,ab,ab0)
         do j=1,jtot
            sptot(j)=spref(j)*abs(scale)
         enddo
      endif
 
      return
      end
 
c ------------------------------------------------------------------ c
 
      subroutine GrnCepBk(Gamma1,Gamma2,Eb,Ec,
     &                    xtot,jtot,spref,xm,ab,ab0)
c     calculates the reflected specutrum according to MZ95
c     ab0 - iron log10 abundance relative to hydrogen
      implicit  none
      integer   jtot
      real      gnr 
      real      xtot(jtot),spref(jtot),xm,ab,ab0
      real      Gamma1,Gamma2,Eb,Ec
      integer   i,j 
      real      xmin,xtransl,xtransh,xrefmax,ymin
      real      pm1(19),pmx(26),apf, ap(3),xmo
      real      xkabs,xnor,xjl,xjd,fjc
      real      y,dym,dy,ddy,y0,sr,srp
      integer   nc,ncp,jmin,jtransl,jtransh,jrefmax
      real      Sigmabfb,SppBk,grxyl,grxyh
 
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
      ap(2)= apf(xmo)
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
         xkabs=Sigmabfb(j,xtot,jtransh,ab,ab0,xnor)
         spref(j)=SppBk(1/xtot(j),Gamma1, Gamma2, Eb, Ec)
     &           *gnr(xkabs,xmo)
 20   continue
c-----------------------
      if(jmin.ge.jtransl) then
         xkabs=Sigmabfb(1,xtot,max(1,jtransh),ab,ab0,xnor)
      endif
      ap(1)=xnor/1.21
      xjl=alog(xtransl)
      xjd=alog(xtransh)-xjl
      do 60 j=jtransl+1,jtransh
         fjc=.5*sin(3.14159*((alog(xtot(j))-xjl)/xjd-.5))
         xkabs=Sigmabfb(j,xtot,jtransh,ab,ab0,xnor)
         y=1/xtot(j)
         spref(j)=SppBk(y,Gamma1, Gamma2, Eb, Ec)*gnr(xkabs,xmo)
         dym=y-ymin
         dy=min(2.,dym)
         ddy=(dy-pmx(26))/(ncp+1)
         y0=y-dy
         sr=SppBk(y0,Gamma1, Gamma2, Eb, Ec)*grxyl(pm1,pmx,dy,y0,ap)
         dy=pmx(26)
         y0=y-dy
         sr=.5*(sr+SppBk(y0,Gamma1, Gamma2, Eb, Ec)
     &        *grxyl(pm1,pmx,dy,y0,ap))
         do 62 i=1,ncp
            dy=dy+ddy
            y0=y-dy
            sr=sr+SppBk(y0,Gamma1, Gamma2, Eb, Ec)
     &           *grxyl(pm1,pmx,dy,y0,ap)
 62      continue
         sr=sr*ddy
         if(dym.gt.2.)then
            ddy=dym-2
            nc=int(ncp*ddy/(dym-pmx(26)))
            ddy=ddy/(nc+1)
            dy=dym
            y0=y-dy
            srp=SppBk(y0,Gamma1, Gamma2, Eb, Ec)
     &          *grxyh(pm1,pmx,dy,y0,ap)
            dy=2
            y0=y-dy
            srp=.5*(srp+SppBk(y0,Gamma1, Gamma2, Eb, Ec)
     &            *grxyh(pm1,pmx,dy,y0,ap))
            do 64 i=1,nc
               dy=dy+ddy
               y0=y-dy
               srp=srp+
     &             SppBk(y0,Gamma1, Gamma2, Eb, Ec)
     &             *grxyh(pm1,pmx,dy,y0,ap)
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
         sr=SppBk(y0,Gamma1, Gamma2, Eb, Ec)*grxyl(pm1,pmx,dy,y0,ap)
         dy=pmx(26)
         y0=y-dy
         sr=.5*(sr+SppBk(y0,Gamma1, Gamma2, Eb, Ec)
     &        *grxyl(pm1,pmx,dy,y0,ap))
         do 32 i=1,ncp
            dy=dy+ddy
            y0=y-dy
            sr=sr+SppBk(y0,Gamma1, Gamma2, Eb, Ec)
     &         *grxyl(pm1,pmx,dy,y0,ap)
 32      continue
         sr=sr*ddy
         if(dym.gt.2.)then
            ddy=dym-2
            nc=int(ncp*ddy/(dym-pmx(26)))
            ddy=ddy/(nc+1)
            dy=dym
            y0=y-dy
            srp=SppBk(y0,Gamma1, Gamma2, Eb, Ec)
     &          *grxyh(pm1,pmx,dy,y0,ap)
            dy=2
            y0=y-dy
            srp=.5*(srp+SppBk(y0,Gamma1, Gamma2, Eb, Ec)
     &            *grxyh(pm1,pmx,dy,y0,ap))
            do 34 i=1,nc
               dy=dy+ddy
               y0=y-dy
               srp=srp+
     &             SppBk(y0,Gamma1, Gamma2, Eb, Ec)
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




 
      REAL FUNCTION SIGMABFB(In,X,N,Ab_,Ab0,Xnor)
c     bound-free cross section for neutral matter with cosmic abundances
c     per hydrogen atom in units of the Thomson cross section
c     ab_ - log10 of abundance of elements heavier than He
c           relative to the Solar abundances
c     ab0 - iron log10 abundance relative to hydrogen
c     x - energy vector of the length n; simgabf is calculated for
c     element in from x

      INCLUDE 'xspec.inc'

      INTEGER i , In , N , nold, ierr
      INTEGER iopac, ixold
      REAL X(*) , abo(17) , ab(17)
      REAL Ab0 , Ab_ , TOTLXSB , Xnor
      REAL ab_in , ab0in
      LOGICAL fl1 , fl2 , firstflag
      REAL fgabnd
      CHARACTER contxt*80
      CHARACTER*4 fgsolr0

      CHARACTER fgsolr*4
      EXTERNAL fgsolr

      SAVE iopac , abo , ab , ixold , nold , fl1 , fl2 , firstflag , 
     &   ab_in , ab0in , fgsolr0

      DATA nold, iopac, ixold/0,2*-1/
 
      DATA ab_in , ab0in/2*999999./
      DATA firstflag , fl1 , fl2/3*.TRUE./
      DATA fgsolr0/'    '/

      ierr = 0

c input standard abundances 
      IF ( (fgsolr().ne.fgsolr0).or.firstflag ) THEN
         firstflag = .FALSE.
         ab_in=999999.
         ab0in=999999.
         fgsolr0=fgsolr()
         abo( 1)=12+log10(fgabnd('H '))
         abo( 2)=12+log10(fgabnd('He'))
         abo( 3)=12+log10(fgabnd('C '))
         abo( 4)=12+log10(fgabnd('N '))
         abo( 5)=12+log10(fgabnd('O '))
         abo( 6)=12+log10(fgabnd('Ne'))
         abo( 7)=12+log10(fgabnd('Na'))
         abo( 8)=12+log10(fgabnd('Mg'))
         abo( 9)=12+log10(fgabnd('Al'))
         abo(10)=12+log10(fgabnd('Si'))
         abo(11)=12+log10(fgabnd('S '))
         abo(12)=12+log10(fgabnd('Cl'))
         abo(13)=12+log10(fgabnd('Ar'))
         abo(14)=12+log10(fgabnd('Ca'))
         abo(15)=12+log10(fgabnd('Cr'))
         abo(16)=12+log10(fgabnd('Fe'))
         abo(17)=12+log10(fgabnd('Ni'))
         ab(1) = abo(1)
         ab(2) = abo(2)
      ENDIF

c     if the energy vector changed from the
c     previous call, recompute the opacity (in totlxs) and reset xold

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
         fl2 = .TRUE.
      ENDIF
      IF ( Ab_.NE.ab_in ) THEN
         DO 100 i = 3 , 17
            ab(i) = abo(i) + Ab_
 100     CONTINUE
         ab(16) = abo(16) + Ab0
         ab_in = Ab_
         ab0in = Ab0
         fl1 = .TRUE.
      ENDIF
      IF ( Ab0.NE.ab0in ) THEN
         ab(16) = abo(16) + Ab0
         ab0in = Ab0
         fl1 = .TRUE.
      ENDIF
c     If either the abundances or the energy vector changed,
c     recompute the opacity
      IF ( fl1 .OR. fl2 ) THEN
         DO 150 i = 1 , N
            MEMR(iopac+i-1) = TOTLXSB(ab,i,X,N,fl1,fl2,Xnor)
 150     CONTINUE
      ENDIF
 
      SIGMABFB = MEMR(iopac+In-1)
 
 999  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         CALL xwrite(contxt, 10)
         SIGMABFB = 0.
      ENDIF

      RETURN
      END
**==TOTLXSB.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
C----------------------------------------------------------------------
C
      FUNCTION TOTLXSB(Ab,In,X,N,Fl1,Fl2,Xnor)
c     *** this version with fully ionized H and He!  AAZ, 18.05.94
c     Optimized version: calculetes cross section for abundances
c     contained in the vector AB at energy in 511.keV units given as
c     an element numbered in of the vector x; n is the length of x.
c     Flag fl1=.true. must indicate changes of AB, and fl2=.true.
c     changes of the energy vector x. For energy >10.keV the law
c     xnor*e^-3 is used, xnor is an output parameter.
c     Note: result in sigma_Thomson units.
C
C     Reference:
C      Monika Balucinska-Church and Dan McCammon
C      "Photoelectric Absorption Cross Sections with Variable Abundances"
C      Ap.J. 400, 699 (1992)
C
C
C     Description:
C     Calculates the effective absorption cross section in units of
c     sigma_Thomson per hydrogen atom at energy E in eV for the
c     abundances of the elements specified in vector AB
C
C     Method:
C     Calls seventeen functions that calculate the mass absorption coeffs
C     in cm**2/g for the following elements: H, He, C, N, O, Ne, Na, Mg,
C     Al, Si, S, Cl, Ar, Ca, Cr, Fe, Ni.
C
C     Deficiencies:
C     Works only in the range of energy from 30 eV to 10,000 eV.
C
C     Bugs:
C     None known -- please report any problems to authors
C
C     Authors:
C     Dan McCammon               (47413::MCCAMMON)
C     Monika Balucinska-Church   (19775::MBC)
C
C     History:
C     8.8.91  - original
C     8.1.92  - modified to remove VAX fortran77 extensions
C     12.3.93 - function HELIUM has been replaced by HELIUM_V2 with
C               Fano line profile
C     21.9.93 - Names of functions simplified: TOTLXS_V2 to TOTLX2,
C               HELIUM_V2 to HEL2.
C     24.1.94 - Names of functions simplified: TOTLX2 to TOTLXS,
C               HEL2 to HELIUM.
C
C     Parameters:
C     E - energy in eV
C     AB - vector of length 17 of Log10 abundances relative to
C          hydrogen = AB(1).  N.B. The order of the elements in AB
C          is as follows: H,He,C,N,O,Ne,Na,Mg,Al,Si,S,Cl,Ar,Ca,Cr,Fe,Ni.
C     TOTLX2 - effective cross section in cm**2/H atom for assumed
C             abundances of cosmic elements.  Note that this will be per
C             hydrogen atom for the relative abundances given in AB.
C
C     Type Definitions:

      INCLUDE 'xspec.inc'

      INTEGER ia, itoth

      INTEGER i, ierr, nsave
      REAL Xnor , xnors
      CHARACTER contxt*80
C
C     Local variables:
      INTEGER k
C          ( index through elements)
      REAL X(*) , ap(17)
      INTEGER In , N
c          x  - energy points
c          in - index for x
c          n  - number of energy points
C          ( mass absorption cross sections in cm**2/g
C           for the cosmic elements -- in the same order as for AB)
C     Import:
      REAL e
C          ( energy in eV)
      REAL Ab(17)
C          ( log10 abundances relative to hydrogen)
 
C     Export:
      REAL TOTLXSB
C          ( effective cross section in cm**2/H atom)
C     External functions: ( mass absorption coefficients for the elements)
      REAL CARBON0
      REAL NITRO0
      REAL OXYGEN0
      REAL NEON0
      REAL SODIUM0
      REAL MAGNES0
      REAL ALUM0
      REAL SILICN0
      REAL SULFUR0
      REAL CHLORN0
      REAL ARGON0
      REAL CALC0
      REAL CHROM0
      REAL IRON0
      REAL NICKEL0
 
 
C     Local constants:
      REAL aw(17)
C          ( atomic weights of the elements)
      REAL av
C          ( Avogadro's number) <- multiplied by sigma_Thomson
      INTEGER numb
C          ( number of elements)
      INTEGER inh
c          ( limiting bin for asymptotic absorption )
 
      LOGICAL fl0 , Fl1 , Fl2
      SAVE ia , itoth, aw , av , numb , ap , fl0 , xnors, nsave
 
      DATA aw/1.00797 , 4.0026 , 12.01115 , 14.0067 , 15.9994 , 20.183 , 
     &     22.9898 , 24.312 , 26.9815 , 28.086 , 32.064 , 35.453 , 
     &     39.94 , 40.08 , 51.996 , 55.847 , 58.71/
 
c      DATA AV /6.022045E23/
      DATA av/.400614/
 
      DATA numb/17/
 
      DATA fl0/.TRUE./

      DATA ia, itoth, nsave/3*-1/
 
C     Start:

      ierr = 0
      IF ( N .NE. nsave ) THEN
         CALL udmget(N*17, 6, ia, ierr)
         contxt = 'Failed to get memory for a'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL udmget(N, 6, itoth, ierr)
         contxt = 'Failed to get memory for toth'
         IF ( ierr .NE. 0 ) GOTO 999
         nsave = N
      ENDIF
 
c     normalizing factor for the asymptotic relation
      IF ( Fl1 ) THEN
C        mass absorption coefficients for the elements at 10.keV
         IF ( fl0 ) THEN
            e = 10000.
c            ap( 1)=HYDRO0(E)
c            ap( 2)=HELIUM0(E)
            ap(1) = 0
            ap(2) = 0
            ap(3) = CARBON0(e)
            ap(4) = NITRO0(e)
            ap(5) = OXYGEN0(e)
            ap(6) = NEON0(e)
            ap(7) = SODIUM0(e)
            ap(8) = MAGNES0(e)
            ap(9) = ALUM0(e)
            ap(10) = SILICN0(e)
            ap(11) = SULFUR0(e)
            ap(12) = CHLORN0(e)
            ap(13) = ARGON0(e)
            ap(14) = CALC0(e)
            ap(15) = CHROM0(e)
            ap(16) = IRON0(e)
            ap(17) = NICKEL0(e)
            fl0 = .FALSE.
         ENDIF
         xnors = 0.
         DO 50 k = 1 , numb
            xnors = xnors + aw(k)*ap(k)*(10.**(Ab(k)-Ab(1)))
 50      CONTINUE
         xnors = xnors/av*(0.01947)**3
         Fl1 = .FALSE.
      ENDIF
      Xnor = xnors
 
      IF ( Fl2 ) THEN
         DO 100 i = 1 , N
c           above 10 keV, use e**-3 law
            IF ( X(i).LT.0.01947 ) THEN
               e = X(i)*511000
c               A( 1,i)=HYDRO0(E)
c               A( 2,i)=HELIUM0(E)
               MEMR(ia+0+17*(i-1)) = 0
               MEMR(ia+1+17*(i-1)) = 0
               MEMR(ia+2+17*(i-1)) = CARBON0(e)
               MEMR(ia+3+17*(i-1)) = NITRO0(e)
               MEMR(ia+4+17*(i-1)) = OXYGEN0(e)
               MEMR(ia+5+17*(i-1)) = NEON0(e)
               MEMR(ia+6+17*(i-1)) = SODIUM0(e)
               MEMR(ia+7+17*(i-1)) = MAGNES0(e)
               MEMR(ia+8+17*(i-1)) = ALUM0(e)
               MEMR(ia+9+17*(i-1)) = SILICN0(e)
               MEMR(ia+10+17*(i-1)) = SULFUR0(e)
               MEMR(ia+11+17*(i-1)) = CHLORN0(e)
               MEMR(ia+12+17*(i-1)) = ARGON0(e)
               MEMR(ia+13+17*(i-1)) = CALC0(e)
               MEMR(ia+14+17*(i-1)) = CHROM0(e)
               MEMR(ia+15+17*(i-1)) = IRON0(e)
               MEMR(ia+16+17*(i-1)) = NICKEL0(e)
               inh = i
            ELSE
               MEMR(itoth+i-1) = X(i)**(-3)
            ENDIF
 100     CONTINUE
         inh = inh + 1
         Fl2 = .FALSE.
      ENDIF
 
      IF ( X(In).LT..01947 ) THEN
         TOTLXSB = 0.
C        loop through elements
         DO 150 k = 1 , numb
            TOTLXSB = TOTLXSB + aw(k)*MEMR(ia+k-1+17*(In-1))
     &                               *(10.**(Ab(k)-Ab(1)))
 150     CONTINUE
         TOTLXSB = TOTLXSB/av
      ELSE
         TOTLXSB = Xnor*MEMR(itoth+In-1)
      ENDIF

 999  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         CALL xwrite(contxt, 10)
         TOTLXSB = 0.
      ENDIF

 
      RETURN
      END



      real function apf (x)
      implicit none
      real x, xinv, xinv2,xinv4
      real c0,c1,c2,c3,c4,c4inv
      data c0/0.802/c1/-1.019/c2/2.528/c3/-3.198/c4/1.457/c4inv/0.030/
* code outside will check for x == 0

      xinv = (c4inv/x)
      xinv2 = xinv*xinv
      xinv4 = xinv2*xinv2
           
      apf = c0 + x*(c1 + x*(c2 + x*(c3 + x*c4))) + xinv4
      return
      end








