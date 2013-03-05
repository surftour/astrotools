**==xspexrav.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
      SUBROUTINE XSPEXRAV(Ear,Ne,Param,Ifl,Photar)
 
      IMPLICIT NONE
      INTEGER Ne , Ifl
      REAL Ear(0:Ne) , Param(*) , Photar(Ne)
 
c     Driver for angle-dependent reflection from an exponentially-cutoff
c     power law and neutral medium.
c     See Magdziarz & Zdziarski, 1995, MNRAS.
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
c     Version with variable iron abundance and new opacities of Balucinska
c     & McCammon (1992, and 1994, private communication). As expected in AGNs,
c     H and He are assumed to be fully ionized.
c
c     This version allows for changes of the vector 'ear' between subsequent
c       calls.
c
c
c     number of model parameters: 7
c     1: Gamma, power law photon index, N_E prop. to E^{-Gamma}
c     2: E_c, the cutoff energy in keV (if E_c=0 there is no cutoff;
c        one needs to change the lower limit for that)
c     3: scale, scaling factor for reflection; if <0, no direct component
c        (scale=1 for isotropic source above disk)
c     4: redshift, z
c     5: abundance of elements heavier than He relative to
c        the solar abundances
c     6: iron abundance relative to the solar iron abundance
c     7: cosine of inclination angle
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
         CALL xwrite('Compton reflection from neutral medium.',5)
         CALL xwrite('See help for details.',5)
         CALL xwrite('If you use results of this model in a paper,',5)
         CALL xwrite
     &  ('please refer to Magdziarz & Zdziarski 1995 MNRAS, 273, 837'
     &              ,5)
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

c   x is in the source frame

      DO 100 i = 0 , Ne
         MEMR(ix+i) = (Ear(i)/511.)*(1.+Param(4))
 100  CONTINUE

      jmax = Ne + 1

c     x is the energy array (units m_e c^2)
c     sppow is the exp. cutoff power law alone (E F_E)
c     spref is the reflected spectrum array (E F_E)
c     sptot is the total spectrum array (E F_E), = sppow if no reflection
c     all dimensions = jmax
      CALL POWREFA(Param,MEMR(ix),jmax,MEMR(isptot),MEMR(isppow),
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
**==powrefa.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
 
      SUBROUTINE POWREFA(Param,X,Jmax,Sptot,Sppow,Spref)
c
      IMPLICIT NONE
      REAL gamma , ec , scale , z_red , normfac , xm , xn
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
c Iron abundance
      ab0 = LOG10(Param(6))
c Cosine of inclination angle
      xm = Param(7)
c--------
c normfac = normalization is at 1 keV in the Earth frame of the cutoff power law
      xn = (1+z_red)/511
 
      IF ( ec.EQ.0. ) THEN
         normfac = 1/SPP0(1/xn,gamma,ec)
         IF ( scale.NE.0. ) THEN
            CALL GREENCA(x,Jmax,Spref,xm,ab_,ab0,SPP0,gamma,ec)
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
               Sppow(j) = SPP0(1/x(j),gamma,ec)*normfac
               Sptot(j) = Sppow(j) + Spref(j)
 160        CONTINUE
         ELSE
            DO 180 j = 1 , Jmax
               Sptot(j) = Spref(j)
 180        CONTINUE
         ENDIF
      ELSE
c ec.NE.0
         normfac = 1/SPPC(1/xn,gamma,ec)      
         IF ( scale.NE.0. ) THEN
            CALL GREENCA(x,Jmax,Spref,xm,ab_,ab0,SPPC,gamma,ec)
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
               Sppow(j) = SPPC(1/x(j),gamma,ec)*normfac
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
**==spp0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
c ------------------------------------------------------------------ c
 
      REAL FUNCTION SPP0(Y,Gamma,Ec)
      IMPLICIT NONE
      REAL Y , Gamma , Ec
      SPP0 = Y**(Gamma-2)
      RETURN
      END
**==sppc.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
      REAL FUNCTION SPPC(Y,Gamma,Ec)
      IMPLICIT NONE
      REAL Y , Gamma , Ec
      SPPC = Y**(Gamma-2)*EXP(-1/(Ec*Y))
      RETURN
      END
**==greenca.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
c ------------------------------------------------------------------ c
 
      SUBROUTINE GREENCA(X,Jmax,Spref,Xm,Ab_,Ab0,SPP,Gamma,Ec)
c     calculates the reflected spectrum
c     ab0 - iron log10 abundance relative to hydrogen
      IMPLICIT NONE
      REAL Ab0 , Ab_ , ddy , dy , dym , Ec , fjc , Gamma , GNR , GRXYH , 
     &     GRXYL , SIGMABFA , SPP , sr , srp , xjd , xjl , xkabs , Xm , 
     &     xmin , xmo
      REAL xnor , xrefmax , xtransh , xtransl , y , y0 , ymin
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
         xkabs = SIGMABFA(j,X,jtransh,Ab_,Ab0,xnor)
         Spref(j) = SPP(1/X(j),Gamma,Ec)*GNR(xkabs,xmo)
 1000 CONTINUE
c-----------------------
      IF ( jmin.GE.jtransl ) xkabs = SIGMABFA(1,X,MAX(1,jtransh),Ab_,
     &                               Ab0,xnor)
      ap(1) = xnor/1.21
      xjl = ALOG(xtransl)
      xjd = ALOG(xtransh) - xjl
      DO 1100 j = jtransl + 1 , jtransh
         fjc = .5*SIN(3.14159*((ALOG(X(j))-xjl)/xjd-.5))
         xkabs = SIGMABFA(j,X,jtransh,Ab_,Ab0,xnor)
         y = 1/X(j)
         Spref(j) = SPP(y,Gamma,Ec)*GNR(xkabs,xmo)
         dym = y - ymin
         dy = MIN(2.,dym)
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
         dy = MIN(2.,dym)
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
**==gnr.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
c ------------------------------------------------------------------ c
c  Reflection factor for nonrelativistic Compton reflection,
c  reflection orders >1 approximated with isotropic phase function,
c  and approximation of H-function by Basko ApJ 223,268,(1978),
c  valid up to 14keV;
c  sig - sigma_bf per hydrogen atom in sigma_Thomson units
c  xm  - cosine of the angle between normal to the slab and the
c        observation direction
 
      REAL FUNCTION GNR(Sig,Xm)
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
 
      GNR = .5*al0*Xm*(.375*((3-2*x2+3*x4)*xl+(3*x2-1)*(.5-Xm))
     &      +((1+(1/al1-1)*(1-LOG(1+al2)/al2))*(1+Xm*SQ3)/(1+Xm*al2)-1)
     &      *xl)
 
      RETURN
      END
**==grxyh.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
c ------------------------------------------------------------------ c
c  Green's functions for Compton reflection with asymtotic
c  absorption included; output is (dy+y0)*G(mu,dy,y0)*W(mu,dy,y0),
c  for .0316<y0<100.; in two parts according to dy=2.
c  grxyh - subroutine for dy>2
c  grxyl - subroutine for pmx(26)<dy<2
c  pm1   - vector of angle-dependent coefficients,
c          set-up by pm1y subroutine
c  pmx   - vector of angle-dependent coefficients,
c          set-up by pmxy subroutine
c  dy    - wavelength shift dy=1/x-1/x0
c  y0    - incident photon wavelength y0=1/x0
c  ap    - absorption parameters
 
      REAL FUNCTION GRXYH(Pm1,Pmx,Dy,Y0,Ap)
      IMPLICIT NONE
      REAL Pm1(19) , Pmx(26) , Dy , Y0 , Ap(3)
      REAL pyx(10) , p1 , p4 , gr , AIC
 
      p1 = LOG(Y0)
      p4 = 1/Y0
      pyx(5) = .230 + .033*p1 + .763*Y0**(-.507) + .007*Y0**(-1.553)
      pyx(6) = 75.180*p4
      pyx(7) = -.100 + .115*p1 + .100*Y0**1.623
      IF ( Y0.LT.1. ) THEN
         pyx(8) = Pmx(19) + Pmx(20)*p1
      ELSE
         pyx(8) = Pmx(21) + Pmx(22)*p1 + Pmx(23)*(Pmx(24)+p4)**Pmx(25)
      ENDIF
      p1 = 1/(Y0+Dy)
      gr = pyx(5)*(1+pyx(6)*p1)**(pyx(7)*(1+pyx(8)*p1))
      p1 = 1 + Dy
      gr = gr*Pm1(18)*Dy**(-1.5)*p1*(1+Pm1(14)/p1)
     &     **(Pm1(15)+Pm1(16)*p1**Pm1(17))
      gr = gr*EXP(AIC(Y0,(Dy-.15),Ap(1)))
 
      GRXYH = gr
 
      RETURN
      END
**==grxyl.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
      REAL FUNCTION GRXYL(Pm1,Pmx,Dy,Y0,Ap)
      IMPLICIT NONE
      REAL Pm1(19) , Pmx(26) , Dy , Y0 , Ap(3)
      REAL pyx(10) , p1 , p2 , p3 , p4 , gr , AIC , xla0 , xla
 
      p1 = LOG(MIN(Y0,10.))
      p2 = p1*p1
      p3 = p2*p1
      pyx(1) = p1*(Pmx(1)+Pmx(2)*p1+Pmx(3)*p2+Pmx(4)*p3)
      pyx(2) = p1*(Pmx(5)+Pmx(6)*p1+Pmx(7)*p2+Pmx(8)*p3)
      pyx(3) = p1*(Pmx(9)+Pmx(10)*p1+Pmx(11)*p2+Pmx(12)*p3)
      pyx(4) = p1*(Pmx(13)+Pmx(14)*p1+Pmx(15)*p2+Pmx(16)*p3+Pmx(17)
     &         *p3*p1)
      pyx(9) = 1 + Pmx(18)*(.326*((Y0*.1)**1.692-1)+p1-2.303)
      IF ( Y0.LT.10. ) THEN
         p1 = Y0 + 2
         p1 = LOG(p1/(Y0+Dy))/LOG(p1/Y0)
         p2 = p1*p1
         p3 = p2*p1
         gr = EXP(pyx(1)+pyx(2)*p1+pyx(3)*p2+pyx(4)*p3)
      ELSE
         p4 = 1/(10+Dy)
         p1 = LOG(12*p4)*5.4848
         p2 = p1*p1
         p3 = p2*p1
         p4 = pyx(9)*(Y0+Dy)*p4
         gr = p4*EXP(pyx(1)+pyx(2)*p1+pyx(3)*p2+pyx(4)*p3)
      ENDIF
      p1 = Pm1(1) + Pm1(2)*Dy
      p2 = (Pm1(3)+Pm1(4)*MAX(0.,(Pm1(5)-Dy))**Pm1(6)+Pm1(7)*(Pm1(8)-Dy)
     &     **Pm1(9))**Pm1(10)/(1+EXP(Pm1(11)*(Pm1(12)-Dy**Pm1(13))))
      gr = gr*MIN(p1,p2)
      xla0 = 1/(1+Ap(1)*Y0**3)
      xla = 1/(1+Ap(1)*(Y0+Dy)**3)
      gr = gr*xla0*(Ap(2)+(1-Ap(2))*(1+(EXP(Ap(3)*Dy)-2)*SQRT(1-xla))
     &     *EXP(AIC(Y0,Dy,Ap(1))))

      GRXYL = gr
 
      RETURN
      END
**==aic.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
      REAL FUNCTION AIC(Y0,Dy,Apc)
      IMPLICIT NONE
      REAL Y0 , Dy , Apc
      REAL yp , a0 , a1 , c0 , c1 , c2 , d0l , d0h , d1l , d1h
      REAL aih , ail
 
      yp = Y0 + Dy
      a0 = Apc
      a1 = Apc**.3333
      c0 = .5/a1
      c1 = -3.4641*c0
      c2 = 0.5774
      d0l = a0*Y0**3
      d0h = a0*yp**3
      d1l = a1*Y0
      d1h = a1*yp
      ail = Y0*(3-LOG(1+d0l)) + c0*LOG((1-d1l+d1l**2)/(1+d1l)**2)
     &      + c1*ATAN(c2*(2*d1l-1))
      aih = yp*(3-LOG(1+d0h)) + c0*LOG((1-d1h+d1h**2)/(1+d1h)**2)
     &      + c1*ATAN(c2*(2*d1h-1))
 
      AIC = aih - ail
 
      RETURN
      END
**==pmxy.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
c ------------------------------------------------------------------ c
c  Angle dependent coefficients for Green's function for Compton
c  reflection;
c  xm  - cosine of the angle between normal to the slab and the
c        observation direction
c  pmx - output coefficients vector
 
      SUBROUTINE PMXY(Xm,Pmx)
      IMPLICIT NONE
      REAL Xm , Pmx(26)
      REAL p2 , p3 , p4
 
      p2 = Xm*Xm
      p3 = p2*Xm
      p4 = p3*Xm
      Pmx(1) = .8827 - .5885*Xm + .7988*p2 - .4117*p3
      Pmx(2) = .0517 + .1076*Xm - .1691*p2 + .0678*p3
      Pmx(3) = .0014 + .0043*Xm
      Pmx(4) = -.0003 - .0027*Xm + .0021*p2
      Pmx(5) = -1.3259 + 1.6214*Xm - 3.7007*p2 + 2.1722*p3 + 
     &         EXP(20.3*(Xm-1.063))
      Pmx(6) = .0790 - .4029*Xm + .6619*p2 + .2210*p3 - 
     &         EXP(17.7*(Xm-1.041))
      Pmx(7) = .0751 - .1848*Xm + .4068*p2 - .4126*p3 + 
     &         EXP(9.8*(Xm-1.185))
      Pmx(8) = -.0020 - .0394*Xm + .1004*p2 - .0597*p3
      Pmx(9) = 3.2953 - 3.6996*Xm + 7.9837*p2 - 4.6855*p3
      Pmx(10) = -.5278 + .9857*Xm - 1.7454*p2 - .3464*p3 + 
     &          EXP(27.6*(Xm-0.959))
      Pmx(11) = -.1919 + .5798*Xm - 1.2879*p2 + 1.4885*p3 - 
     &          EXP(18.3*(Xm-0.986))
      Pmx(12) = .0200 + .0832*Xm - .0333*p2 - .2370*p3 + 
     &          EXP(16.6*(Xm-1.086))
      Pmx(13) = -2.2779 + 2.4031*Xm - 4.0733*p2 + 1.9499*p3 - 
     &          EXP(19.6*(Xm-0.968))
      Pmx(14) = .4790 - 1.0166*Xm + 3.1727*p2 - 4.0108*p3 + 3.0545*p4 - 
     &          EXP(30.4*(Xm-0.957))
      Pmx(15) = .1122 - .3580*Xm + .4985*p2 - .3750*p3 - .5349*p4 + 
     &          EXP(31.2*(Xm-0.972))
      Pmx(16) = -.0160 - .0471*Xm - .0536*p2 + .2006*p3 + .2929*p4 - 
     &          EXP(30.6*(Xm-1.009))
      Pmx(17) = .0005 + .0021*Xm - .0447*p2 + .1749*p3 - .2303*p4 + 
     &          EXP(16.9*(Xm-1.130))
      Pmx(18) = .006 + .089*Xm - .102*p2 + .056*p3
      Pmx(19) = .618 - .830*Xm
      Pmx(20) = .128 - .132*Xm
      Pmx(21) = .632 - .875*Xm
      Pmx(22) = -.672 - .286*Xm + .717*p2 - .481*p3
      Pmx(23) = .0126 - .0160*Xm + .0077*p2
      Pmx(24) = .0111 + .0030*Xm - .0014*p2
      Pmx(25) = -2.437 - .328*Xm - .260*p2 + .279*p3
      Pmx(26) = 1 - SQRT(1-(Xm-.05)**2)
 
      RETURN
      END
**==pm1y.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
c ------------------------------------------------------------------ c
c  Angle dependent coefficients for Green's function for Compton
c  reflection of incident photon at energy x0=1;
c  xm  - cosine of the angle between normal to the slab and the
c        observation direction
c  pm1 - output coefficients vector
 
      SUBROUTINE PM1Y(Xm,Pm1)
      IMPLICIT NONE
      REAL Xm , Pm1(19)
      REAL p2 , p3 , p4
 
      p2 = .911 - .549*(1-Xm)**1.471
      p3 = .254 - .041*Xm**(-3.798)
      Pm1(1) = p2 - 2*p3
      Pm1(2) = p3
      Pm1(3) = .161 + .439*Xm**.189 - .791*Xm**7.789
      Pm1(4) = MAX(0.,(-.871+1.740*Xm**.176-1.088*Xm**.907))
      Pm1(5) = .934 + .054*Xm**(-.666)
      p4 = MIN(.524,Xm)
      Pm1(6) = 4.647 - 11.574*p4 + 11.046*p4*p4
      Pm1(7) = .012 + .199*Xm**.939 + .861*Xm**8.999
      Pm1(8) = 2.150
      p2 = Xm*Xm
      p3 = p2*Xm
      Pm1(9) = -1.281 + 2.374*Xm - 4.332*p2 + 2.630*p3
      Pm1(10) = 2.882 + .035*Xm**(-1.065)
      Pm1(11) = 30.28 + 69.29*Xm
      Pm1(12) = 1.037*(1-SQRT(1-.874*p2))**.123
      Pm1(13) = .123
      p4 = p3*Xm
      Pm1(14) = 56.50 - 27.54*Xm + 671.8*p2 - 1245.*p3 + 708.4*p4
      Pm1(15) = -.897 - .826*Xm + 2.273*p2 - .984*p3
      Pm1(16) = 1.240 - 1.297*Xm
      Pm1(17) = -.490 + 1.181*Xm - 1.038*p2
      Pm1(18) = 1.858*(Xm**2.033+0.745*Xm**1.065)
      Pm1(19) = 1 - SQRT(1-(Xm-.05)**2)
 
      RETURN
      END
**==sigmabfa.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
c ------------------------------------------------------------------- c
 
      REAL FUNCTION SIGMABFA(In,X,N,Ab_,Ab0,Xnor)
c     bound-free cross section for neutral matter with cosmic abundances
c     per hydrogen atom in units of the Thomson cross section
c     ab_ - log10 of abundance of elements heavier than He
c           relative to the Solar abundances
c     ab0 - iron log10 abundance relative to hydrogen
c     x - energy vector of the length n; simgabf is calculated for
c     element in from x

      INCLUDE '../../inc/xspec.inc'

      INTEGER i , In , N , nold, ierr
      INTEGER iopac, ixold
      REAL X(*) , abo(17) , ab(17)
      REAL Ab0 , Ab_ , TOTLXSA , Xnor, xnorsav
      REAL ab_in , ab0in
      LOGICAL fl1 , fl2 , firstflag
      REAL fgabnd
      CHARACTER contxt*80
      CHARACTER*4 fgsolr0

      CHARACTER fgsolr*4
      EXTERNAL fgsolr

      SAVE iopac , abo , ab , ixold , nold , fl1 , fl2 , firstflag , 
     &   ab_in , ab0in , fgsolr0, xnorsav

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
            MEMR(iopac+i-1) = TOTLXSA(ab,i,X,N,fl1,fl2,Xnor)
 150     CONTINUE
         xnorsav = Xnor
      ENDIF
 
      Xnor = xnorsav
      SIGMABFA = MEMR(iopac+In-1)
 
 999  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         CALL xwrite(contxt, 10)
         SIGMABFA = 0.
      ENDIF

      RETURN
      END
**==totlxsa.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
C----------------------------------------------------------------------
C
      FUNCTION TOTLXSA(Ab,In,X,N,Fl1,Fl2,Xnor)
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

      INCLUDE '../../inc/xspec.inc'

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
      REAL TOTLXSA
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
         TOTLXSA = 0.
C        loop through elements
         DO 150 k = 1 , numb
            TOTLXSA = TOTLXSA + aw(k)*MEMR(ia+k-1+17*(In-1))
     &                               *(10.**(Ab(k)-Ab(1)))
 150     CONTINUE
         TOTLXSA = TOTLXSA/av
      ELSE
         TOTLXSA = Xnor*MEMR(itoth+In-1)
      ENDIF

 999  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         CALL xwrite(contxt, 10)
         TOTLXSA = 0.
      ENDIF

 
      RETURN
      END
**==alum0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
 
 
C-------------------------------------------------------------------------
C
C  File XSCTNS.FOR
C
C        This set of subroutines calculates the photoelectric absorption
C   cross sections for the elements H, He, C, N, O, Ne, Na, Mg, Al, Si,
C   S, Cl, A, Ca, Cr, Fe, and Ni.  The result is in cm**2/g, given the
C   photon energy in eV.  These functions are valid only over the energy
C   range 30 - 10,000 eV, but do NOT check that the input energy is
C   within the valid range.  These functions are called by TOTLX2 to
C   calculate the total effective cross section, given a set of relative
C   abundances.  They can also be used by themselves.
C
C   history : 26 May, 1993 - original HELIUM function replaced by HEL2 ;
C              HELIUM not deleted from the program (19775::MBC)
C
C           : 21 Sep, 1993 - Function HELIUM now deleted, all remnants
C             of VAX Fortran 77 extensions removed, file renamed
C             XSCTN2 (19775::MBC)
C
C           : 22 Sep, 1993 - bug in FANO routine (called by HEL2) corrected
C             that resulted in incorrect line shape and energy (19775::MBC)
C
C           : 24 Jan, 1994 - name HEL2 changed to HELIUM
C   References:
C
C      Monika Balucinska-Church and Dan McCammon
C      "Photoelectric Absorption Cross Sections with Variable Abunances"
C      Ap.J. 400, 699 (1992)
C
C      All data are from:
C      B. L. Henke, P. Lee, T. J. Tanaka, R. L. Shimabukuro and B. K.
C      Fujikawa, 1982, Atomic Data and Nuclear Data Tables, vol 27, p 1.
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C  Real Function: ALUM
C  Source: Atomic & Nuclear Data Tables, January 1982
C
C  Description:
C      Calculates mass absorption coefficient (mu/rho) for aluminum.
C  History:  updated below L-edge (72.78 eV) - Aug 30, 1991 (MBC)
C
C  Usage:  FUNCTION ALUM(E)
C      E = Energy in eV.
C      ALUM returns mu/rho in cm**2/gm.
C
C  COMMON Blocks:
C      none
C  IMPLICIT
C      none
C  Subroutines called by ALUM:
C      none
C
C--------------------------------------------------------------------------
      FUNCTION ALUM0(E)
      IMPLICIT NONE
 
      REAL E , elog , x , ALUM0
 
      elog = ALOG(E)
      IF ( E.LT.72.78 ) THEN
         x = 26.90487 + (3.-9.135221)*elog + 1.175546*elog*elog
 
      ELSEIF ( E.LT.1559.9 ) THEN
         x = -38.1232 + 29.5161*elog - 4.45416*elog*elog + 
     &       0.226204*elog*elog*elog
      ELSE
         x = 14.6897 + 4.22743*elog - 0.344185*elog*elog + 
     &       8.18542E-3*elog*elog*elog
      ENDIF
 
      ALUM0 = EXP(x)/(E*E*E)
 
      RETURN
      END
**==argon0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C------------------------------------------------------------------------
C------------------------------------------------------------------------
C------------------------------------------------------------------------
C
C  REAL FUNCTION: ARGON
C  Source: Atomic & Nuclear Data Tables, January 1982
C  Description: ARGON calculates the mass absorbtion coefficient of argon.
C
C ****          works well for energies above 40 eV !!!
C
C  History: updated below L-edge (245 eV) - Aug 30, 1991 (MBC)
C
C  Usage: FUNCTION ARGON(E)
C      E = Energy in eV
C      ARGON returns the mass absorption cross section in cm**2/g
C
C  COMMON BLOCKS:
C      NONE
C
C  IMPLICIT
C      NONE
C  SUBROUTINES CALLED BY ARGON:
C      NONE
C
C---------------------------------------------------------------------------
      FUNCTION ARGON0(E)
      IMPLICIT NONE
 
      REAL E , elog , x , ARGON0
 
      elog = ALOG(E)
      IF ( E.LT.245.0 ) THEN
         x = -330.3509 + (267.7433+3.)*elog - 78.90498*elog*elog + 
     &       10.35983*(elog**3) - 0.5140201*(elog**4)
 
      ELSEIF ( E.LT.3202.9 ) THEN
         x = -5.71870 + (8.85812*elog) + (-0.307357*elog*elog)
     &       + (0.00169351*(elog**3)) + (-0.0138134*(elog**4))
     &       + (0.00120451*(elog**5))
      ELSE
         x = 19.1905 + (2.74276*elog) + (-0.164603*elog*elog)
     &       + (0.00165895*elog*elog*elog)
      ENDIF
 
      ARGON0 = EXP(x)/(E*E*E)
 
      RETURN
      END
**==calc0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C--------------------------------------------------------------------------
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Real Function: CALC
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for calcium
C
C     History: original Aug 6, 1991 (MBC)
C
C     Usage:   FUNCTION CALC(E)
C              E = Energy in eV
C              CALC returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT
C            none
C-------------------------------------------------------------------------
C
      FUNCTION CALC0(E)
      IMPLICIT NONE
C
C
      REAL E , elog , x , CALC0
C
C
      elog = ALOG(E)
 
 
      IF ( E.LT.349.31 ) THEN
         x = -873.972 + (865.5231+3.)*elog - 339.678*elog*elog + 
     &       66.83369*(elog**3) - 6.590398*(elog**4)
     &       + 0.2601044*(elog**5)
 
      ELSEIF ( E.LT.4038.1 ) THEN
         x = -3449.707 + (2433.409+3.)*elog - 682.0668*elog*elog + 
     &       95.3563*(elog**3) - 6.655018*(elog**4)
     &       + 0.1854492*(elog**5)
 
      ELSE
         x = 18.89376 + (3.-0.2903538)*elog - 0.1377201*elog*elog
 
      ENDIF
      CALC0 = EXP(x)/(E*E*E)
 
 
      RETURN
      END
**==carbon0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
C
C  REAL FUNCTION: CARBON
C  Source: Atomic & Nuclear Data Tables, Jan. 1982
C
C  Description: CARBON calculates the mass absorbtion cross-section of carbon.
C
C  USAGE: FUNCTION CARBON(E)
C      E = Energy in EV
C      CARBON returns the mass absorption cross-section in cm**2/g
C
C  COMMON BLOCKS:
C      NONE
C
C  IMPLICIT
C      NONE
C
C  SUBROUTINES CALLED BY CARBON:
C      NONE
C
C-------------------------------------------------------------------------
      FUNCTION CARBON0(E)
      IMPLICIT NONE
C
C
      REAL E , elog , x , CARBON0
C
C
      elog = ALOG(E)
      IF ( E.LT.284.0 ) THEN
         x = 8.74161 + (7.13348*elog) + (-1.14604*elog*elog)
     &       + (0.0677044*elog*elog*elog)
      ELSE
         x = 3.81334 + (8.93626*elog) + (-1.06905*elog*elog)
     &       + (0.0422195*elog*elog*elog)
      ENDIF
 
      CARBON0 = EXP(x)/(E*E*E)
      RETURN
      END
**==chlorn0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
C     Real Function: CHLORN
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for chlorine
C
C     History: original Aug 6, 1991 (MBC)
C
C     Usage:   FUNCTION CHLORN(E)
C              E = Energy in eV
C              CHLORN returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT
C            none
C-------------------------------------------------------------------------
C
      FUNCTION CHLORN0(E)
      IMPLICIT NONE
C
C
      REAL E , elog , x , CHLORN0
C
C
      elog = ALOG(E)
 
 
      IF ( E.LT.202.0 ) THEN
         x = 6253.247 + (3.-8225.248)*elog + 4491.675*elog*elog - 
     &       1302.145*(elog**3) + 211.4881*(elog**4)
     &       - 18.25547*(elog**5) + 0.6545154*(elog**6)
 
      ELSEIF ( E.LT.2819.6 ) THEN
         x = -233.0502 + (143.9776+3.)*elog - 31.12463*elog*elog + 
     &       2.938618*(elog**3) - 0.104096*(elog**4)
 
      ELSE
         x = -23.74675 + (14.50997+3.)*elog - 1.857953*elog*elog + 
     &       6.6208832E-2*(elog**3)
 
      ENDIF
      CHLORN0 = EXP(x)/(E*E*E)
 
 
      RETURN
      END
**==chrom0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
C     Real Function: CHROM
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for chromium
C
C     History: original Aug 5, 1991 (MBC)
C
C     Usage:   FUNCTION CHROM(E)
C              E = Energy in eV
C              CHROM returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT
C            none
C-------------------------------------------------------------------------
C
      FUNCTION CHROM0(E)
      IMPLICIT NONE
C
C
      REAL E , elog , x , CHROM0
C
C
      elog = ALOG(E)
 
 
      IF ( E.LT.598.0 ) THEN
         x = -0.4919405 + (12.66939+3.)*elog - 5.199775*elog*elog + 
     &       1.086566*(elog**3) - 0.1196001*(elog**4)
     &       + 5.2152011E-3*(elog**5)
 
      ELSEIF ( E.LT.691.0 ) THEN
         x = 27.29282 + (3.-2.703336)*elog
 
      ELSEIF ( E.LT.5988.8 ) THEN
         x = -15.2525 + (13.23729+3.)*elog - 1.966778*elog*elog + 
     &       8.062207E-2*(elog**3)
 
      ELSE
         x = 8.307041 + (2.008987+3.)*elog - 0.2580816*elog*elog
 
      ENDIF
      CHROM0 = EXP(x)/(E*E*E)
 
 
      RETURN
      END
**==helium0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C------------------------------------------------------------------------------
C-----------------------------------------------------------------------------
C
C     Real Funcion : HELIUM
C     Source : Marr, G. V., and West, J. B., Atomic and Nuclear Data Tables,
C                (1976) 18, 497.
C             Oza, D. H., (1986), Phys. Rev. A, 33,  824.
C             Fernley, J. A., Taylor, K. T., and Seaton, M. J., (1987),
C                J. Phys. B., 20, 6457.
C
C     Description :
C     calculates mass absorption coefficient (mu/rho) in cm2/g for neutral
C     helium for the given energy in eV.
C     Cross sections come from experimental data compiled by Marr and
C     West (Atomic Data and Nuclear Data Tables (1976) 18, 497).
C     The four strongest autoionization resonances are taken into account;
C     numbers come from Oza (Phys Rev A (1986), 33, 824), and Fernley et al.
C     (J. Phys B (1987) 20, 6457).
C
C     Deficiencies :
C     works in the energy range from 30 eV to 10,000 eV
 
C     Bugs :
C     if any are found please report to the authors
C
C     History :
C     this subroutine replaces the previous version of HELIUM which
C     calculated mass absoprtion coefficients based on Henke's data
C     (Henke, B. L., et al., (1982), Atomic and Nuclear Data Tables, 27, 1).
C     This version of HELIUM returns mass  absorption coefficients which
C     are in better agreement with the best experiments as well as
C     theoretical models (see Chen, W. F., Cooper, G., and Brion, C. E.,
C     (1991), Phys. Rev. A, 44, 186).  This fortran-77 version of the
C     subroutine is based on Jelinsky's program written in C
C     (EUVE Archive)
C
C     History :
C     Jan 4, 1993 : original (19775::MBC)
C     Feb 23, 1993 : comments added and modified to remove VAX
C                    fortran 77 extensions
C
C     Sep 21, 1993 : further remnants of VAX fortran 77 extensions
C                    have been removed (19775::MBC)
C
C     Sep 22, 1993 : bug in the FANO routine has been removed (19775::MBC)
C
C     Usage : FUNCTION HELIUM(E)
C            E = Energy in eV
C
C     Common Blocks :
C           none
C
C     Implicit :
C           none
C
C     Functions called by HELIUM
C           FANO
C
C------------------------------------------------------------------------------
      FUNCTION HELIUM0(E)
      IMPLICIT NONE
 
C    Type definitions :
C     IMPLICIT NONE
C    Global variables :
C    Structure definitions :
C    Function declarations :
      REAL FANOO
C    Local constants :
      INTEGER IP
C         ( index through loop )
      PARAMETER (IP=8)
      INTEGER IF
C         ( index through loop )
      PARAMETER (IF=4)
      REAL AV
C         ( Avogadro's number )
      PARAMETER (AV=6.022045E23)
      REAL AW
C         ( atomic weight of HYDROgen )
      PARAMETER (AW=4.0026E0)
C    Local variables :
      REAL lambda
C          ( wavelength in Angstroms)
      REAL x
      REAL y
      REAL sigma
C          ( cross section in cm2/atom)
      INTEGER i
C          ( index trough loop)
C     Import :
      REAL E
C          ( energy in eV)
C     Export :
      REAL HELIUM0
C          ( cross section in cm**2/g)
C    Local data :
      REAL c1(IP)
      REAL c2(IP)
      REAL q(IF)
      REAL nu(IF)
      REAL gamma(IF)
 
C          ( polynomial coefficients for Marr and West data)
      DATA c1/ - 2.953607E1 , 7.083061E0 , 8.678646E-1 , -1.221932E0 , 
     &     4.052997E-2 , 1.317109E-1 , -3.265795E-2 , 2.500933E-3/
 
C          ( polynomial coefficients for Marr and West data )
      DATA c2/ - 2.465188E1 , 4.354679E0 , -3.553024E0 , 5.573040E0 , 
     &     -5.872938E0 , 3.720797E0 , -1.226919E0 , 1.576657E-1/
 
C          ( parameters Q for resonances (Fernley et al. 1987) )
      DATA q/2.81E0 , 2.51E0 , 2.45E0 , 2.44E0/
 
C          ( parameters NU for resonances (Oza 1986) )
      DATA nu/1.610E0 , 2.795E0 , 3.817E0 , 4.824E0/
 
C          ( parameters GAMMA for resonances (Oza 1986) )
      DATA gamma/2.64061E-3 , 6.20116E-4 , 2.56061E-4 , 1.320159E-4/
 
 
C     Start :
 
      lambda = 12398.54E0/E
      x = ALOG10(lambda)
 
      IF ( lambda.GT.503.97E0 ) THEN
         HELIUM0 = 0.E0
         RETURN
 
      ELSEIF ( lambda.LT.46.E0 ) THEN
         y = 0.E0
         DO 50 i = 1 , IP
            y = y + c2(i)*(x**(i-1))
 50      CONTINUE
 
      ELSE
         y = 0.E0
         DO 100 i = 1 , IP
            y = y + c1(i)*(x**(i-1))
 100     CONTINUE
 
         DO 150 i = 1 , IF
            y = y + ALOG10(FANOO(q(i),nu(i),gamma(i),lambda))
 150     CONTINUE
 
      ENDIF
 
      sigma = 10.E0**y
      HELIUM0 = sigma*AV/AW
 
      END
**==fanoo.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C
C----------------------------------------------------------------------------
C----------------------------------------------------------------------------
C    Real Function FANO
C
C    Source :
C             Fernley, J. A., Taylor, K. T., and Seaton, M. J., (1987),
C                J. Phys. B., 20, 6457.
C             Oza, D. H., (1986), Phys. Rev. A, 33,  824.
C
C    Description :
C     returns a Fano line profile for use by HELIUM. The form of
C     this line shape is taken from Fernley, Taylor and Seaton (above).
C
C    Deficiencies :
C     works in the energy range from 30 eV to 10,000 eV
C
C    Bugs :
C    if any are found please report to the authors
C
C    22 Sep, 1993 - bug fixed in calculations of EPSI - devided by 1.807317
C    not added (19775::MBC)
C
C    History :
C     Jan 4, 1993 : original, based on Jelinsky's program written in
C     C (EUVE Archive) - (19775::MBC)
C
C     Feb 23, 1993 : comments added and modified to remove VAX
C     fortran 77 extensions (19775::MC)
C
C     Sep 21, 1993 : further remnants of VAX fortran 77 extensions
C     removed (19775::MBC)
C
C
C    Common Blocks :
C           none
C    Implicit :
C           none
C
C    Functions called by FANO
C           none
C
C----------------------------------------------------------------------------
C
      FUNCTION FANOO(A,B,C,Lambda)
      IMPLICIT NONE
 
C    Type definitions :
C     IMPLICIT NONE
C    Global variables :
C    Structure definitions :
C    Function declarations :
C    Local constants :
C    Local variables :
      REAL eps
C          ( energy in Rydbergs )
      REAL epsi
      REAL x
C          ( log_10 of wavelength in Angstroms )
C     Import :
      REAL A
C          ( Q coefficient (Fernley et al. 1987) )
      REAL B
C          ( NU coefficient (Oza 1986) )
      REAL C
C          ( GAMMA coefficient (Oza 1986) )
      REAL Lambda
C          ( wavelength in Angstroms )
C     Export :
      REAL FANOO
C    Start :
 
      eps = 911.2671E0/Lambda
      epsi = 3.0E0 - 1.E0/(B*B) + 1.807317
      x = 2.0*(eps-epsi)/C
      FANOO = (x-A)*(x-A)/(1.0E0+x*x)
 
      END
**==hydro0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C-----------------------------------------------------------------------
C       REAL FUNCTION HYDRO(E)
C
C       DATE:3/6/84
C       AUTHOR: ANGIE BETKER
C       SOURCE: ATOMIC AND NUCLEAR DATA TABLES, JANUARY 1982
C
C        History: modified: 6/5/87 - J. Bloch - Created F77 Vax/VMS version.
C                 updated Aug 30, 1991 (MBC)
C
C
C       DESCRIPTION: HYDRO CALCULATES MU/RHO FOR HYDROGEN IN CM**2/GM
C
C       COMMON BLOCKS:
C               NONE
C
C       SUBROUTINES CALLED:
C               NONE
C
C       IMPLICIT
C               NONE
C
C---------------------------------------------------------------------------
C
      FUNCTION HYDRO0(E)
      IMPLICIT NONE
C
C
      REAL E , elog , x , HYDRO0
C
C
C        IF (E.LT.180) HYDRO0 = (10 ** 10.11) * (E ** -3.03)
C        IF ((E.GE.180).AND.(E.LT.750)) HYDRO0=(10**10.54)*(E**-3.23)
C        IF ((E.GE.750).AND.(E.LT.6000)) HYDRO0=(10**10.94)*(E**-3.37)
C        IF (E.GE.6000) HYDRO0 = (10 ** 10.42) * (E ** -3.24)
 
      elog = ALOG(E)
 
      x = 21.46941 + (3.-2.060152)*elog - 0.1492932*elog*elog + 
     &    5.4634294E-3*(elog**3)
 
      HYDRO0 = EXP(x)/(E*E*E)
 
      RETURN
      END
**==iron0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C
C----------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Real Function: IRON
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for iron
C     History: original Aug 6, 1991 (MBC)
C
C     Usage:   FUNCTION IRON(E)
C              E = Energy in eV
C              IRON returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT:
C            none
C
C-------------------------------------------------------------------------
C
      FUNCTION IRON0(E)
      IMPLICIT NONE
C
C
      REAL E , elog , x , IRON0
C
C
      elog = ALOG(E)
 
      IF ( E.LT.707.4 ) THEN
         x = -15.07332 + (18.94335+3.)*elog - 4.862457*elog*elog + 
     &       0.5573765*(elog**3) - 3.0065542E-2*(elog**4)
     &       + 4.9834867E-4*(elog**5)
 
      ELSEIF ( E.LT.7111.2 ) THEN
         x = -253.0979 + (135.4238+3.)*elog - 25.47119*elog*elog + 
     &       2.08867*(elog**3) - 6.4264648E-2*(elog**4)
 
      ELSE
         x = -1.037655 + (4.022304+3.)*elog - 0.3638919*elog*elog
 
      ENDIF
      IRON0 = EXP(x)/(E*E*E)
 
 
      RETURN
      END
**==magnes0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C-------------------------------------------------------------------------
C-------------------------------------------------------------------------
C     Real Function: MAGNES
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for magnesium
C     History: original Aug 6, 1991 (MBC)
C
C     Usage:   FUNCTION MAGNES(E)
C              E = Energy in eV
C              MAGNES returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT:
C            none
C-------------------------------------------------------------------------
C
      FUNCTION MAGNES0(E)
      IMPLICIT NONE
C
C
      REAL E , elog , x , MAGNES0
C
C
      elog = ALOG(E)
 
 
      IF ( E.LT.49.45 ) THEN
         x = 7.107172 + (0.7359418+3.)*elog
 
      ELSEIF ( E.LT.1303.4 ) THEN
         x = -81.32915 + (62.2775+3.)*elog - 15.00826*elog*elog + 
     &       1.558686*(elog**3) - 6.1339621E-2*(elog**4)
 
      ELSE
         x = -9.161526 + (10.07448+3.)*elog - 1.435878*elog*elog + 
     &       5.2728362E-2*(elog**3)
 
      ENDIF
      MAGNES0 = EXP(x)/(E*E*E)
 
 
      RETURN
      END
**==neon0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C
C------------------------------------------------------------------------
C----------------------------------------------------------------------
C  REAL FUNCTION: NEON
C  Source: Atomic and Nuclear Data Tables, Jan. 1982
C
C  Description:  NEON calculates the mass absorption coefficient
C      for NEON0 gas.
C
C  Usage:  REAL FUNCTION NEON(E)
C      E = Energy in eV
C      NEON returns the mass absorption cross section in cm**2/g.
C
C  COMMON BLOCKS:
C      NONE
C
C  IMPLICIT
C      NONE
C
C  SUBROUTINES CALLED BY NEON:
C      NONE
C
C------------------------------------------------------------------------------
      FUNCTION NEON0(E)
      IMPLICIT NONE
C
C
      REAL E , elog , x , NEON0
C
C
      elog = ALOG(E)
      IF ( E.LT.867. ) THEN
         x = -3.04041 + (13.0071*elog) + (-1.93205*elog*elog)
     &       + (0.0977639*elog*elog*elog)
      ELSE
         x = 17.6007 + (3.29278*elog) + (-0.263065*elog*elog)
     &       + (5.68290E-3*elog*elog*elog)
      ENDIF
      NEON0 = EXP(x)/(E*E*E)
      RETURN
      END
**==nickel0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
C
C  REAL FUNCTION: NICKEL
C  Source: Atomic & Nuclear Data Tables, January 1982
C
C  Description: NICKEL calculates the mass absorbtion coefficient for nickel.
C  History: updated below L-edge (853.6 eV) - Aug 30, 1991 (MBC)
C
C  Usage: REAL FUNCTION NICKEL(E)
C      E = Energy in eV
C      NICKEL returns the mass absorption cross section in cm**2/g.
C
C  COMMON BLOCKS:
C      NONE
C  IMPLICIT
C      NONE
C
C  SUBROUTINES CALLED BY NICKEL:
C      NONE
C
C
C---------------------------------------------------------------------------
      FUNCTION NICKEL0(E)
      IMPLICIT NONE
 
      REAL E , elog , x , NICKEL0
 
      elog = ALOG(E)
 
      IF ( E.LT.853.6 ) THEN
         x = -7.919931 + (11.06475+3.)*elog - 1.935318*elog*elog + 
     &       9.3929626E-2*(elog**3)
 
      ELSEIF ( E.LT.8331.6 ) THEN
         x = 3.71129 + (8.45098*elog) + (-0.896656*elog*elog)
     &       + (0.0324889*elog*elog*elog)
      ELSE
         x = 28.4989 + (0.485797*elog)
      ENDIF
 
      NICKEL0 = EXP(x)/(E*E*E)
 
      RETURN
      END
**==nitro0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
C
C  REAL FUNCTION: NITRO
C  Source: Atomic & Nuclear Data Tables, January 1982
C
C  Description:  NITRO calculates the mass absorption cross-section of
C        nitrogen as a function of energy.
C
C  Usage: REAL FUNCTION NITRO(E)
C      E = Energy in eV
C      NITRO returns the mass absorption cross-section in cm**2/g.
C
C  COMMON BLOCKS:
C      NONE
C
C  IMPLICIT
C      NONE
C
C  SUBROUTINES CALLED BY NITRO:
C      NONE
C
C--------------------------------------------------------------------------
      FUNCTION NITRO0(E)
      IMPLICIT NONE
C
C
      REAL E , elog , x , NITRO0
C
C
      elog = ALOG(E)
      IF ( E.LT.401. ) THEN
         x = 9.24058 + (7.02985*elog) + (-1.08849*elog*elog)
     &       + (0.0611007*elog*elog*elog)
      ELSE
         x = -13.0353 + (15.4851*elog) + (-1.89502*elog*elog)
     &       + (0.0769412*elog*elog*elog)
      ENDIF
      NITRO0 = EXP(x)/(E*E*E)
      RETURN
      END
**==oxygen0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C
C------------------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C  REAL FUNCTION: OXYGEN
C
C  Source:        Atomic & Nuclear Data Tables, January 1982
C
C  Description:  OXY Calculates the mass absorption cross-section of oxygen
C      as a function of energy.
C  History: updated above K-edge (531.7 eV) - Aug 30, 1991 (MBC)
C
C  Usage: REAL FUNCTION OXYGEN(E)
C      E = Energy in eV
C      OXYGEN returns the mass absorption cross-section in cm**2/g.
C
C  COMMON BLOCKS:
C      NONE
C
C  IMPLICIT
C      NONE
C  SUBROUTINES CALLED BY OXYGEN:
C      NONE
C
C------------------------------------------------------------------------
      FUNCTION OXYGEN0(E)
      IMPLICIT NONE
 
      REAL E , elog , x , OXYGEN0
 
      elog = ALOG(E)
      IF ( E.LT.531.7 ) THEN
         x = 2.57264 + (10.9321*elog) + (-1.79383*elog*elog)
     &       + (0.102619*elog*elog*elog)
      ELSE
         x = 16.53869 + (0.6428144+3.)*elog - 0.3177744*elog*elog + 
     &       7.9471897E-3*(elog**3)
 
      ENDIF
 
      OXYGEN0 = EXP(x)/(E*E*E)
 
      RETURN
      END
**==silicn0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
C       FUNCTION SILICN
C
C       Source: Atomic and Nuclear Data Tables, January 1982
C
C       Description: SILICN calculates the mass absorption cross section
C                for silicon in cm**2/g.
C       History: updated Aug 30, 1991 (MBC)
C                updated March 4, 1992 (MBC)
C
C       COMMON BLOCKS:
C               NONE
C
C       IMPLICIT
C               none
C
C       SUBROUTINES CALLED:
C               NONE
C
C---------------------------------------------------------------------------
      FUNCTION SILICN0(E)
      IMPLICIT NONE
 
      REAL E , elog , x , SILICN0
C
      elog = ALOG(E)
 
      IF ( E.LT.100.6 ) THEN
         x = -3.066295 + (7.006248+3.)*elog - 0.9627411*elog*elog
      ELSEIF ( E.LT.1840.0 ) THEN
         x = -182.7217 + (125.061+3.)*elog - 29.47269*elog*elog + 
     &       3.03284*(elog**3) - 0.1173096*(elog**4)
      ELSE
         x = -33.39074 + (18.42992+3.)*elog - 2.385117*elog*elog + 
     &       8.887583E-2*(elog**3)
 
      ENDIF
 
      SILICN0 = EXP(x)/(E*E*E)
      RETURN
      END
**==sodium0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C-------------------------------------------------------------------------
C------------------------------------------------------------------------
C     Real Function: SODIUM
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for sodium
C
C     History: original Aug 6, 1991 (MBC)
C
C     Usage:   FUNCTION SODIUM(E)
C              E = Energy in eV
C              SODIUM returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT:
C            none
C-------------------------------------------------------------------------
C
      FUNCTION SODIUM0(E)
      IMPLICIT NONE
C
C     IMPLICIT NONE
C
      REAL E , elog , x , SODIUM0
C
C
      elog = ALOG(E)
 
      IF ( E.LT.1071.7 ) THEN
         x = -2737.598 + (2798.704+3.)*elog - 1009.892*elog*elog + 
     &       87.16455*(elog**3) + 43.20644*(elog**4)
     &       - 15.27259*(elog**5) + 2.180531*(elog**6)
     &       - 0.1526546*(elog**7) + 4.3137977E-3*(elog**8)
 
      ELSE
         x = 1.534019 + (6.261744+3.)*elog - 0.9914126*elog*elog + 
     &       3.5278253E-2*(elog**3)
 
 
      ENDIF
      SODIUM0 = EXP(x)/(E*E*E)
 
 
      RETURN
      END
**==sulfur0.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     Real Function: SULFUR
C     Source: Atomic & Nuclear Data Tables, January 1982
C
C     Description:
C     Calculates mass absorption coefficient (mu/rho) for sulfur
C
C     History: original Aug 6, 1991 (MBC)
C
C     Usage:   FUNCTION SULFUR(E)
C              E = Energy in eV
C              SULFUR returns mu/rho in cm**2/g
C
C     COMMON Blocks:
C            none
C
C     IMPLICIT
C            none
C-------------------------------------------------------------------------
C
      FUNCTION SULFUR0(E)
      IMPLICIT NONE
C
C
      REAL E , elog , x , SULFUR0
C
C
      elog = ALOG(E)
 
 
      IF ( E.LT.165.0 ) THEN
         x = 598.2911 + (3.-678.2265)*elog + 308.1133*elog*elog - 
     &       68.99324*(elog**3) + 7.62458*(elog**4)
     &       - 0.3335031*(elog**5)
 
      ELSEIF ( E.LT.2470.5 ) THEN
         x = 3994.831 + (3.-3693.886)*elog + 1417.287*elog*elog - 
     &       287.9909*(elog**3) + 32.70061*(elog**4)
     &       - 1.968987*(elog**5) + 4.9149349E-2*(elog**6)
 
      ELSE
         x = -22.49628 + (14.24599+3.)*elog - 1.848444*elog*elog + 
     &       6.6506132E-2*(elog**3)
 
      ENDIF
      SULFUR0 = EXP(x)/(E*E*E)
 
 
      RETURN
      END
C
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
 
 
