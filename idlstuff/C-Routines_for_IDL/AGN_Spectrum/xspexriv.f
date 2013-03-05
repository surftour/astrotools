**==xspexriv.spg  processed by SPAG 4.50J  at 11:24 on 27 Feb 1996
      SUBROUTINE XSPEXRIV(Ear,Ne,Param,Ifl,Photar)
 
      IMPLICIT NONE
      INTEGER Ne , Ifl
      REAL Ear(0:Ne) , Param(*) , Photar(Ne)
 
c     Driver for angle-dependent reflection from an exponentially-cutoff
c     power law and ionized medium.
c     It needs to be compiled together with xspexrav.f (neutral reflection).
c
c     See Magdziarz & Zdziarski, 1995, MNRAS.
c     See Zdziarski et al., 1995, ApJL 438, L63 for description of
c     calculation of ionization (based on Done et al. 1992, ApJ 395, 275).
c     The abundances are defined by the command abund
c
c     The output spectrum is the sum of the cutoff power law and the
c     reflection component.
c     The reflection component alone can be obtained
c     for scale (see below) = rel_refl < 0. Then the actual
c     reflection normalization is |scale|. Note that you need to
c     change then the limits of rel_refl. The range of rel_refl in that case
c     should exclude zero (as then the direct component appears).
c
c     If E_c=0, there is no cutoff in the power law.
c
c     This version allows for changes of the vector 'ear' between subsequent
c       calls.
c
c
c     number of model parameters:9
c     1: Gamma, power law photon index, N_E prop. to E^{-Gamma}
c     2: E_c, the cutoff energy in keV (if E_c=0 there is no cutoff;
c        one needs to change the lower limit for that)
c     3: scale, scaling factor for reflection; if <0, no direct component
c        (scale=1 for isotropic source above disk)
c     4: redshift, z
c     5: abundance of elements heavier than He relative to
c        the solar abundances
c     6: iron abundance relative to the solar abundances
c     7: cosine of inclination angle
c     8: disk temperature in K
c     9: disk ionization parameter = L/nR^2
c     algorithm:
c          a[E*(1+z)]=E^{-Gamma}*exp(-E/E_c)+scale*reflection
c     Normalization is the photon flux at 1 keV (photons keV^-1 cm^-2 s^-1)
c     of the cutoff power law only (without reflection)
c     and in the earth frame.
c
c
      INCLUDE 'xspec.inc'

      INTEGER isptot, ispref, isppow, ix, iphotaux      

      INTEGER i , j , jmax, ierr, nesave
      LOGICAL firstcall
      CHARACTER contxt*127

      SAVE isptot, ispref, isppow, ix, iphotaux      
      SAVE firstcall, nesave

      DATA isptot, ispref, isppow, ix, iphotaux, nesave/6*-1/
      DATA firstcall/.TRUE./

      ierr = 0
 
      IF ( firstcall ) THEN

         CALL xwrite('Compton reflection from ionized medium.',5)
         CALL xwrite('See help for details.',5)
         CALL xwrite('If you use results of this model in a paper,',5)
         CALL xwrite
     & ('please refer to Magdziarz & Zdziarski 1995 MNRAS, 273, 837',5)
         CALL xwrite('Questions: Andrzej Zdziarski, aaz@camk.edu.pl',5)

         firstcall = .FALSE.
      ENDIF

      IF ( Ne .NE. nesave ) THEN
         CALL udmget(Ne+1, 6, isptot, ierr)
         contxt = 'Failed to get memory for sptot'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL udmget(Ne+1, 6, ispref, ierr)
         contxt = 'Failed to get memory for spref'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL udmget(Ne+1, 6, isppow, ierr)
         contxt = 'Failed to get memory for sppow'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL udmget(Ne+1, 6, ix, ierr)
         contxt = 'Failed to get memory for x'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL udmget(Ne+1, 6, iphotaux, ierr)
         contxt = 'Failed to get memory for photaux'
         IF ( ierr .NE. 0 ) GOTO 999
         nesave = Ne
      ENDIF

 
c  x is in the source frame

      DO 100 i = 0 , Ne
         MEMR(ix+i) = (Ear(i)/511.)*(1.+Param(4))
 100  CONTINUE

      jmax = Ne + 1

c     x is the energy array (units m_e c^2)
c     sppow is the exp. cutoff power law alone (E F_E)
c     spref is the reflected spectrum array (E F_E)
c     sptot is the total spectrum array (E F_E), = sppow if no reflection
c     all dimensions = jmax
      CALL POWREFJ(Param,MEMR(ix),jmax,MEMR(isptot),MEMR(isppow),
     &             MEMR(ispref))
      DO 200 j = 0 , jmax-1
         MEMR(iphotaux+j) = MEMR(isptot+j)/Ear(j)**2
 200  CONTINUE
      DO 300 i = 1 , Ne
         Photar(i) = .5*(MEMR(iphotaux+i)+MEMR(iphotaux+i-1))
     &                 *(Ear(i)-Ear(i-1))
 300  CONTINUE
 
 999  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         CALL xwrite(contxt, 10)
      ENDIF
 
      RETURN
      END
**==powrefj.spg  processed by SPAG 4.50J  at 11:24 on 27 Feb 1996
 
 
      SUBROUTINE POWREFJ(Param,X,Jmax,Sptot,Sppow,Spref)
c
      IMPLICIT NONE
      REAL gamma , ec , scale , z_red , normfac , xm , xn
      REAL temp , xi
      REAL ab_ , ab0 , SPP0 , SPPC
      REAL Param(*) , Sptot(*) , Sppow(*) , Spref(*) , X(*)
      INTEGER Jmax , j
      EXTERNAL SPP0 , SPPC
 
c--------
C  INPUT:
c the power law photon index:
      gamma = Param(1)
c the break energy in units of m_e c^2
      ec = Param(2)/511
c the scaling factor for the reflected spectrum;
c 1 corresponds to seeing equal
c contributions from the reflected and direct spectra
      scale = Param(3)
c redshift
      z_red = Param(4)
c abundance
      ab_ = LOG10(Param(5))
c relative Iron abundance
      ab0 = LOG10(Param(6))
c Cosine of inclination angle
      xm = Param(7)
c reflecting matter temperature
      temp = Param(8)
c ionization parameter
      xi = Param(9)
c--------

c normalization is at 1 keV in the Earth frame of the power law only
      xn = (1+z_red)/511
 
      IF ( ec.EQ.0. ) THEN
         normfac = 1/SPP0(1/xn,gamma,ec)
         IF ( scale.NE.0. ) THEN
            CALL GREENCJ(X,Jmax,Spref,xm,ab_,ab0,SPP0,gamma,ec,temp,xi)
            DO 120 j = 1 , Jmax
               Spref(j) = Spref(j)*ABS(scale)*normfac
 120        CONTINUE
         ELSE
            DO 140 j = 1 , Jmax
               Spref(j) = 0
 140        CONTINUE
         ENDIF
         IF ( scale.GE.0 ) THEN
            DO 160 j = 1 , Jmax
               Sppow(j) = SPP0(1/X(j),gamma,ec)*normfac
               Sptot(j) = Sppow(j) + Spref(j)
 160        CONTINUE
         ELSE
            DO 180 j = 1 , Jmax
               Sptot(j) = Spref(j)
 180        CONTINUE
         ENDIF
      ELSE
         normfac = 1/SPPC(1/xn,gamma,ec)
         IF ( scale.NE.0. ) THEN
            CALL GREENCJ(X,Jmax,Spref,xm,ab_,ab0,SPPC,gamma,ec,temp,xi)
            DO 200 j = 1 , Jmax
               Spref(j) = Spref(j)*ABS(scale)*normfac
 200        CONTINUE
         ELSE
            DO 220 j = 1 , Jmax
               Spref(j) = 0
 220        CONTINUE
         ENDIF
         IF ( scale.GE.0 ) THEN
            DO 240 j = 1 , Jmax
               Sppow(j) = SPPC(1/X(j),gamma,ec)*normfac
               Sptot(j) = Sppow(j) + Spref(j)
 240        CONTINUE
         ELSE
            DO 260 j = 1 , Jmax
               Sptot(j) = Spref(j)
 260        CONTINUE
         ENDIF
      ENDIF
 
      RETURN
      END
**==greencj.spg  processed by SPAG 4.50J  at 11:24 on 27 Feb 1996
 
c ------------------------------------------------------------------ c
 
      SUBROUTINE GREENCJ(X,Jmax,Spref,Xm,Ab_,Ab0,SPP,Gamma,Ec,Temp,Xi)
c     calculates the reflected spectrum
c     ab0 - iron log10 abundance relative to hydrogen
      IMPLICIT NONE
      REAL Ab0 , Ab_ , ddy , dy , dym , Ec , fjc , Gamma , GNRJ , 
     &     GRXYH , GRXYL , SIGMABFJ , SPP , sr , srp , Temp , Xi , xjd , 
     &     xjl , xkabs , Xm
      REAL xmin , xmo , xnor , xrefmax , xtransh , xtransl , y , y0 , 
     &     ymin
      INTEGER i , j , Jmax , jmin , jrefmax , jtransh , jtransl , nc , 
     &        ncp
      REAL X(*) , Spref(*)
      REAL pm1(19) , pmx(26) , apf, ap(3)
      EXTERNAL SPP
 
c  precision factor for Green's function integration
      ncp = 100
c  ranges for different methods
      xmin = 2.E-4
      xtransl = 0.01957
      xtransh = 0.02348
      ymin = .03333
c  angle dependent parameters
      xmo = MAX(.05,MIN(.95,Xm))
      CALL PM1Y(xmo,pm1)
      CALL PMXY(xmo,pmx)
      ap(2) = apf(xmo)
      ap(3) = 0.381
      xrefmax = 1/(ymin+pmx(26))
 
      jmin = 0
      DO 100 j = 1 , Jmax
         IF ( X(j).GT.xmin ) GOTO 200
         jmin = j
 100  CONTINUE
 200  jtransl = 0
      DO 300 j = MAX(1,jmin) , Jmax
         IF ( X(j).GT.xtransl ) GOTO 400
         jtransl = j
 300  CONTINUE
 400  jtransh = 0
      DO 500 j = MAX(1,jtransl) , Jmax
         IF ( X(j).GT.xtransh ) GOTO 600
         jtransh = j
 500  CONTINUE
 600  jrefmax = 0
      DO 700 j = MAX(1,jtransh) , Jmax
         IF ( X(j).GT.xrefmax ) GOTO 800
         jrefmax = j
 700  CONTINUE
 
 800  DO 900 j = 1 , jmin
         Spref(j) = 0
 900  CONTINUE
c-----------------------
      DO 1000 j = jmin + 1 , jtransl
         xkabs = SIGMABFJ(j,X,jtransh,Gamma,Temp,Xi,Ab_,Ab0,xnor)
         Spref(j) = SPP(1/X(j),Gamma,Ec)*GNRJ(xkabs,xmo)
 1000 CONTINUE
c-----------------------
      IF ( jmin.GE.jtransl ) xkabs = SIGMABFJ(1,X,MAX(1,jtransh),Gamma,
     &                               Temp,Xi,Ab_,Ab0,xnor)
      ap(1) = xnor/1.21
      xjl = ALOG(xtransl)
      xjd = ALOG(xtransh) - xjl
      DO 1100 j = jtransl + 1 , jtransh
         fjc = .5*SIN(3.14159*((ALOG(X(j))-xjl)/xjd-.5))
         xkabs = SIGMABFJ(j,X,jtransh,Gamma,Temp,Xi,Ab_,Ab0,xnor)
         y = 1/X(j)
         Spref(j) = SPP(y,Gamma,Ec)*GNRJ(xkabs,xmo)
         dym = y - ymin
         dy = MIN(2.0,dym)
         ddy = (dy-pmx(26))/(ncp+1)
         y0 = y - dy
         sr = SPP(y0,Gamma,Ec)*GRXYL(pm1,pmx,dy,y0,ap)
         dy = pmx(26)
         y0 = y - dy
         sr = .5*(sr+SPP(y0,Gamma,Ec)*GRXYL(pm1,pmx,dy,y0,ap))
         DO 1050 i = 1 , ncp
            dy = dy + ddy
            y0 = y - dy
            sr = sr + SPP(y0,Gamma,Ec)*GRXYL(pm1,pmx,dy,y0,ap)
 1050    CONTINUE
         sr = sr*ddy
         IF ( dym.GT.2. ) THEN
            ddy = dym - 2
            nc = INT(ncp*ddy/(dym-pmx(26)))
            ddy = ddy/(nc+1)
            dy = dym
            y0 = y - dy
            srp = SPP(y0,Gamma,Ec)*GRXYH(pm1,pmx,dy,y0,ap)
            dy = 2
            y0 = y - dy
            srp = .5*(srp+SPP(y0,Gamma,Ec)*GRXYH(pm1,pmx,dy,y0,ap))
            DO 1060 i = 1 , nc
               dy = dy + ddy
               y0 = y - dy
               srp = srp + SPP(y0,Gamma,Ec)*GRXYH(pm1,pmx,dy,y0,ap)
 1060       CONTINUE
            sr = sr + srp*ddy
         ENDIF
         Spref(j) = (.5-fjc)*Spref(j) + (.5+fjc)*sr*X(j)
 1100 CONTINUE
c-----------------------
      DO 1200 j = jtransh + 1 , jrefmax
         y = 1/X(j)
         dym = y - ymin
         dy = MIN(2.0,dym)
         ddy = (dy-pmx(26))/(ncp+1)
         y0 = y - dy
         sr = SPP(y0,Gamma,Ec)*GRXYL(pm1,pmx,dy,y0,ap)
         dy = pmx(26)
         y0 = y - dy
         sr = .5*(sr+SPP(y0,Gamma,Ec)*GRXYL(pm1,pmx,dy,y0,ap))
         DO 1150 i = 1 , ncp
            dy = dy + ddy
            y0 = y - dy
            sr = sr + SPP(y0,Gamma,Ec)*GRXYL(pm1,pmx,dy,y0,ap)
 1150    CONTINUE
         sr = sr*ddy
         IF ( dym.GT.2. ) THEN
            ddy = dym - 2
            nc = INT(ncp*ddy/(dym-pmx(26)))
            ddy = ddy/(nc+1)
            dy = dym
            y0 = y - dy
            srp = SPP(y0,Gamma,Ec)*GRXYH(pm1,pmx,dy,y0,ap)
            dy = 2
            y0 = y - dy
            srp = .5*(srp+SPP(y0,Gamma,Ec)*GRXYH(pm1,pmx,dy,y0,ap))
            DO 1160 i = 1 , nc
               dy = dy + ddy
               y0 = y - dy
               srp = srp + SPP(y0,Gamma,Ec)*GRXYH(pm1,pmx,dy,y0,ap)
 1160       CONTINUE
            sr = sr + srp*ddy
         ENDIF
         Spref(j) = sr*X(j)
 1200 CONTINUE
c----------------------
      DO 1300 j = jrefmax + 1 , Jmax
         Spref(j) = 0
 1300 CONTINUE
 
      RETURN
      END
**==gnrj.spg  processed by SPAG 4.50J  at 11:24 on 27 Feb 1996
 
c ------------------------------------------------------------------ c
c  Reflection factor for nonrelativistic Compton reflection,
c  reflection orders >1 approximated with isotropic phase function,
c  and approximation of H-function by Basko ApJ 223,268,(1978),
c  valid up to 14keV;
c  sig - sigma_bf per hydrogen atom in sigma_Thomson units
c  xm  - cosine of the angle between normal to the slab and the
c        observation direction
 
      REAL FUNCTION GNRJ(Sig,Xm)
      IMPLICIT NONE
      REAL Sig , Xm
      REAL al0 , al1 , al2 , xl , x2 , x4 , SQ3
      PARAMETER (SQ3=1.732051)
 
      al0 = 1/(1+Sig/1.21)
      al1 = SQRT(1-al0)
      al2 = SQ3*al1
      xl = LOG(1+1/Xm)
      x2 = Xm*Xm
      x4 = x2*x2
 
      GNRJ = .5*al0*Xm*(.375*((3-2*x2+3*x4)*xl+(3*x2-1)*(.5-Xm))
     &       +((1+(1/al1-1)*(1-LOG(1+al2)/al2))*(1+Xm*SQ3)/(1+Xm*al2)-1)
     &       *xl)
 
      RETURN
      END
**==sigmabfj.spg  processed by SPAG 4.50J  at 11:24 on 27 Feb 1996
 
c ------------------------------------------------------------------- c
 
      FUNCTION SIGMABFJ(In,X,N,Gamma0,Temp0,Xi0,Ab_0,Ab0,Asc0)
c     bound-free cross section for ionized matter with cosmic abundances
c     and variable iron abundance (in units of solar iron abundance)
c     per hydrogen atom in units of the Thomson cross section
c     e is in units of 511 keV

      INCLUDE 'xspec.inc'

      REAL ab , Ab0 , ab_ , Ab_0 , asc , Asc0 , gamma , Gamma0 , 
     &     SIGMABFJ , temp , Temp0 , WARMABJ , xas , xi , Xi0
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
         asc = 1.503E2*WARMABJ(xas,gamma,temp,xi,ab_,ab,firstflag,
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
            MEMR(iopac+i-1) = 1.503E2*WARMABJ(X(i),gamma,temp,xi,ab_,ab,
     &                firstflag,.false.)
 100     CONTINUE
         fl1 = .FALSE.
      ENDIF
 
      Asc0 = asc
      SIGMABFJ = MEMR(iopac+In-1)

 999  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         CALL xwrite(contxt, 10)
         SIGMABFJ = 0.
      ENDIF
 
      RETURN
      END

**==warmabj.spg  processed by SPAG 4.50J  at 11:24 on 27 Feb 1996

c -------------------------------------------------------------------

      FUNCTION WARMABJ(X0,Gamma,Temp,Xi,Ab_,Ab0,Firstflag,qHHe)

      REAL WARMABJ , X0 , Gamma , Temp , Xi , Ab_ , Ab0
      LOGICAL Firstflag, qHHe

      INCLUDE 'xspec.inc'

      REAL clcwrm
      REAL aw(20), abo(20), e(800)
      INTEGER z(20), iab0, imax
      INTEGER isigma, iion, iint, iarec, inum, iratio, ierr
      CHARACTER*72 contxt

      COMMON /aprmdt/ isigma, iion, z, aw, abo, e, iab0, imax

      SAVE iint, iarec, inum, iratio

      DATA isigma, iion, iint, iarec, inum, iratio /6*-1/

      ierr = 0

C Get the memory for the arrays required by CLCWRM

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

      WARMABJ = CLCWRM(X0,Gamma,Temp,Xi,Ab_,Ab0,Firstflag,qHHe,
     &                 MEMR(isigma),MEMR(iion),MEMR(iint),MEMR(iarec),
     &                 MEMD(inum),MEMD(iratio), z, aw, abo, e, iab0,
     &                 imax)

 999  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         CALL xwrite(contxt, 10)
         WARMABJ = 0.
      ENDIF

      RETURN
      END
 
c -------------------------------------------------------------------
 
      FUNCTION CLCWRM(X0,Gamma,Temp,Xi,Ab_,Ab0,Firstflag,qHHe,sigma,
     &                ion,int,arec,num,ratio,z,aw,abo,e,iab0,imax)

c 	calculate cross-sections using Daltabuit and Cox, use cosmic
c       abundances, and radiative recombination rates from
c       Aldrovandi and Pequignot 1973 and explicit formula
c       for recombination to hydrogenlike atoms
c       result in units of 10^-22
      IMPLICIT NONE
      REAL CLCWRM , X0 , Gamma , Temp , Xi , Ab_ , Ab0
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

      CALL ldapmr(z, aw, abo, ion, sigma, e, iab0, imax)

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
 
      CLCWRM = sigmaeff

      RETURN
      END
 
c ------------------------------------------------------------------- c
 
      SUBROUTINE ldapmr(z, aw, abo, ion, sigma, e, iab0, imax)

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
      IF ( ios .NE. 0 ) THEN
         WRITE(wrtstr,'(a,i5)') 'Failed to open iond_mm.dat : ios = ', 
     &                          ios
         CALL xwrite(wrtstr, 10)
         RETURN
      ENDIF
 
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
      IF ( ios .NE. 0 ) THEN
         WRITE(wrtstr,'(a,i5)') 'Failed to open mansig.dat : ios = ', 
     &                          ios
         CALL xwrite(wrtstr, 10)
         RETURN
      ENDIF

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

c ------------------------------------------------------------------- c

      function el_name(z)

      CHARACTER*2 el_name,elt_long(28)
      INTEGER z

      DATA elt_long
     &        /'H ','He',
     &         'Li','Be','B ','C ','N ','O ','F ','Ne',
     &         'Na','Mg','Al','Si','P ','S ','Cl','Ar',
     &         'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni'/ 

      el_name=elt_long(z)

      RETURN
      END

c ------------------------------------------------------------------- c

