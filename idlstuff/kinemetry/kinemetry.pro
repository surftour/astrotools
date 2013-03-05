;#######################################################################
;
; Copyright (C) 2005, Davor Krajnovic
;
; If you have found this software useful for your
; research, I would appreciate an acknowledgment to use of the
; `Kinemetry: a generalisation of photometry to the higher moments
; of the line-of-sight velocity distribution' by Krajnovic et al. (2006)'
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;#######################################################################
;NAME:
;   KINEMETRY
;
;PURPOSE:
;   Perform harmonic expansion of 2D maps of observed kinematic
;   moments (velocity, velocity dispersion, h3, h4...) along the best
;   fitting ellipses (either fixed or free to change along the radii).
;
;EXPLANATION:
;   This program is a generalisation of the ellipse fitting method used
;   in photometry (e.g. Jedrzejewsky 1987) to the odd moments of the 
;   line-of-sight velocity distribution (LOSVD). The even moments are 
;   treated as in in photometry, while there is a modification for the 
;   odd moments. The method assumes that along the best fitting ellipse
;   the odd moment is well approximated by a simple cosine term, as in 
;   the tilted ring method by Begeman (1987). We use interpolation to 
;   sample the kinematic moment along the ellipses, as done in all 
;   ellipse-fitting implementations for photometry, in the few-pixels
;   regime.
;
;   For a given radius a best fit ellipse is described by flattening
;   'q' (= 1 - ellipticity) and its position angle 'PA'. The ellipse
;   parameters are found similarly to Jedrzejewsky (1987) photometry
;   approach, but in the case of odd moments by minimising 'a1', 'a3',
;   and 'b3' coefficients (or coeffs. next to sin(x), sin(3x), cos(3x))
;   of the Fourier expansion of a kinematic profile extracted along an 
;   ellipse. This is possible because a small error in q produces a 
;   non-zero b3 coefficient, and a small error in PA produces non-zero 
;   a1, a3, and b3 coefficients. Errors in the position of the centre 
;   produces nonzero even terms (a2,b2), which are ignored for the moment.
;
;   The determination of q and PA is done in two steps. First a
;   grid of q and PA is specified (not very dense but covering the
;   whole parameters space) and for each combination of the
;   ellipse parameters a kinematic profile is extracted. Using a
;   least-squares fit (singular value decomposition) the Fourier
;   coefficients are determined. The combination (q,PA), for which
;   a1^2+a3^2+b3^2 is the smallest, is used as the input for the
;   second step, which consists in a non-linear search for best
;   parameters performed by MPFIT.
;
;   After determination of the best fitting ellipse parameters
;   (described by q and PA for a given radius), a new Fourier expansion
;   is performed. The results are the Fourier coefficients and 
;   reconstructed kinematic moment maps.
;
;CALLING SEQUENCE:
;   KINEMETRY, xbin, ybin, moment, rad, pa, q, cf, $
;     [ NTRM=NTRM, ERROR=error, SCALE=scale, NAME=name, $
;	PAQ=paq, NPA=npa, NQ=nq, RANGEQ=rangeq, ALL=all, EVEN=even,$
;	VSYS=VSYS, VELCIRC=velCirc, VELKIN=velkin, ER_CF=er_cf, $
;	ER_PA=er_pa, ER_q=er_q, XELLIP=xellip, YELLIP=yellip, $
;	RING=ring, RADIUS = radius, PLOT=plot, VERBOSE=verbose ]
;
;INPUTS:
;   XBIN   - 1D array with X coordinates describing the map
;   YBIN   - 1D array with Y coordinates describing the map
;   MOMENT - 1D array with kin.moment (e.g. velocity) values 
;            at XBIN, YBIN positions
;
;OPTIONAL INPUTS KEYWORDS:   
;   /NTRM  - scalar specifying the number of terms for the harmonic analysis 
;	     of the profile extracted along the best fitting ellipse. Default
;	     value is 6 odd terms, which means the following terms will be 
;	     used in the expansion: a1, b1, a3, b3, a5 and a5, or the terms
;	     corresponding to sin(x),cos(x),sin(3x),cos(3x),sin(5x),cos(5x).
;   ERROR  - 1D array with errors to VELBIN values. If this
;            keyword is specified then the program will calculate
;	     formal (1sigma) errors of the coefficients which are
;	     returned in ER_PA, ER_Q (determined by MPFIT), ER_CF 
;	     (determined by SVD) variables. 
;   SCALE  - scalar specifying the pixel scale on the map. If not set, the
;            SAURON pixel scale (0.8 arcsec) is assumed.
;   NAME   - name of the object (used by VERBOSE keyword and for internal
;            plotting)
;   PAQ    - 2 element or 2*NRAD element vector specifying position angle (PA)
;	     and flattening (q) of the ellipses in that order (kept constant).
;	     It is possible to specify a set of PA and q values (that 
;	     correspond to given radii (see RADIUS keyword)), for which one
;	     wants to get the Fourier coefficients. In this case PAQ should
;	     be set as follows: PAQ=[PA1,Q1, PA2,Q2...., PAnrad,Qnrad]
;   NPA    - scalar specifying the number of PA used to crudely estimate
;            the parameters of the best fit ellipse before entering
;            MPFIT. Default value is 21. To speed up the process and for
;	     quick tests it is useful to use a small number (e.g 5). 
;            Complicated maps may require a bigger number (e.g. 41).
;   NQ     - scalar specifying the number of q used to crudely estimate
;            the parameters of the best fit ellipse before entering
;            MPFIT. Default value is 21. To speed up the process and for
;	     quick tests it is useful to use a small number (e.g 5). 
;            Complicated maps may require a bigger number (e.g. 41).
;   RANGEQ - 2 element vector specifying the min and max value for 
;	     flattening Q. Default values are 0.2 and 1.0.
;   /ALL   - If this keyword is set then the harmonic analysis of the rings 
;	     will include both even and odd terms. If this keyword is set,
;	     and NTRM = n then the following terms are used in the expansion:
;	     a1, b2, a2, b2, a3, b3,...., an, bn (or coeffs nex to: sin(x),
;	     cos(x),sin(2x),cos(2x),sin(3x),cos(3x),...,sin(nx),cos(nx))
;   /EVEN  - set this keyword to do kinemetry on even kinematic moments. 
;	     In this case, kinemetry reduces to photometry and the best
;	     fitting ellipse is obtained by minimising a1, b1, a2, b2
;	     terms. When this keyword is set, keyword /ALL is automatically
;	     set and NTRM should be increased (e.g. NTRM=10 will use the 
;	     following terms in the expansion: a1, b2, a2, b2, a3, b3, a4, b4
;	     (or coeffs. next to sin(x),cos(x),sin(2x), cos(2x),sin(3x),
;	     cos(3x),sin(4x),cos(4x)))
;   /VSYS  - if this keyword is set the zeroth term (a0) is not extracted.
;	     (for odd moments).This is useful for determining rotation curves.
;	     One can first run kinemetry without setting this keyword to find 
;	     the systemic velocity (given as cf[*,0]). Then subtract
;	     the systemic velocity form the velocity map and re-run kinemetry 
;	     with /vsys set. In this case the zeroth terms will be zero. For 
;	     completeness, it is also possible to input VSYS, e.g. VSYS=10. 
; 	     The zeroth term will not be calculated, but it will be set to 10 
;	     in output. 
;   RING   - scalar specifying desired radius of the first ring. Set this
;	     keyword to a value at which the extraction should begin. This 
;	     is useful in case of ring-like structures sometimes observed in 
;	     HI data.
;   RADIUS - 1D array with values specifying the lenght of the semi-major axis
;	     at which the data (kin.profile) should be extracted from the map 
;	     for the kinemetric analisys. Te values should be in pixel space 
;            (not in physical units such as arcsec).  If this keyword is set, 
;	     the values are coopied into the output variable: RAD.
;   /PLOT  - If this keyword is set, diagnostic plots are shown for each
;            radii: 
;		- the best ellipse (overploted on kin.moment map), 
;		- position of minimum (green diamond) in PA-Q diagram 
;		  (the parameters of the best fit ellipse) determined by
;		  MPFIT, where the initial (input to MPFIT) values of PA
;		  and q are presented by dots and the colours
;	          show the Chi2 square contours (linearly interpolated 
;		  between the PA,Q points),
;		- fit to kin.profile (white is DATA, red is the FIT, where
;		  FIT is given by a0+b1*cos(x) for odd, and a0 for even 
;		  moments), 
;		- residuals (DATA - FIT), and overplotted higher order
;		  terms (green: a1,a3 and b3, red: a1,a3,b3,a5 and b5; 
;		  for the /EVEN case - green: a1,b1,a2,b2, red:a1,b1,
;		  a2,b2,a4,b4)
;   /VERBOSE - set this keyword to print status of the fit on screen
;
;OUTPUT:
;   RAD    - 1D array with radii at which kin.profiles were extracted
;   PA     - 1D array with position angle of the best fitting ellipses,
;            PA is first determined on an interval PA=[-90,90], where
;            PA=0 along positive x-axis. Above x-axis PA > 0 and below
;            x-axis Pa < 0. PA does not differentiate between receding
;            and approaching sides of (velocity) maps. This is
;            transformed to the usual East of North system, where the 
;	     East is the negative x-axis, and the North is the positive 
;	     y-axis. For odd kin.moments PA is measured from the North 
;	     to the receding (positive) side of the galaxy (which is 
;	     detected by checking the sign of the cos(theta) term. For
;	     the even terms it is left degenerate to 180 degrees rotation.
;   Q      - 1D array with flattening of the best fitting ellipses
;            (q=1-ellipticity), defined on interval q=[0.2,1]
;   CF     - 2D array containing coefficients of the Fourier expansion
;            for each radii cf=[Nradii, Ncoeff]. For example: 
;	     a0=cf[*,0], a1=cf[*,1], b1=cf[*,2]....
;
;OPTIONAL OUTPUT KEYWORDS:
;   VELKIN - 1D array of reconstructed kin.moment using NTRM harmonic
;            terms at positions XBIN,YBIN, obtained by linear interpolation
;	     from points given in XELLIP and YELLIP keywords
;   VELCIRC - 1D array containg 'circular velocity' or a0 + b1*cos(theta)
;            at positions XBIN, YBIN (velcirc = a0, in case of EVEN moments), 
;	     obtained by linear interpolation from points given in XELLIP 
;	     and YELLIP keywords;   
;   ER_PA  - 1D array of 1sigma errors to the ellipse position angle
;   ER_Q   - 1D array of 1sigma errors to the ellipse axial ratio
;   ER_CF  - 2D array containing 1 sigma errors to the coefficients 
;	     of the Fourier expansion for each radii
;   XELLIP - 1D array with X coordintes of the best fitting ellipses
;   YELLIP - 1D array with Y coordintes of the best fitting ellipses
;
;RESTRICTIONS:
;   Speed and robustness of the program depends on the number of NQ
;   and NPA that define the initial grid. Small grid is fast, but
;   not precise in complicated cases.
;
;NOTE: determination of the centre and penalisation towards extraction 
;      on the circles is not yet included in this distribution. 
;
;REQUIRED ROUTINES:
;   MPFIT: by C.B. Markwardt from http://astrog.physics.wisc.edu/~craigm/idl/
;   RANGE: by M. Cappellari (included in kinemetry distribution)
;   SAURON_COLORMAP: by Michele Cappellari & Eric Emsellem, 2001
;			(included in kinemetry distribution)
;EXAMPLE:
;    1) 
;    Run kinemetry on a velocity map defined with coordinate arrays:
;    Xbin, Ybin, and a velocity array: Velbin. Desired outputs are: 
;    position angle and flattening of the ellipses, harmonic terms:
;    a0, a1, b1, a3, b3, a5, b5 and a reconstructed map of the circular
;    velocity: 
; 
;    KINEMETRY, Xbin, Ybin, Velbin, rad, pa, q, cf, NTRM=6, VELCIRC=velcirc,$
;	        /plot, /verbose 
;
;    2)
;    Run kinemetry on a velocity map starting at radius=5". Desired outputs 
;    are position angle and flattening of the ellipses, harmonic terms:
;    a0, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a reconstructed map 
;    of the circular velocity and a map reconstructed with all terms: 
;
;    KINEMETRY, Xbin, Ybin, Velbin, rad, pa, q, cf, NTRM=10, /ALL, $
;    VELCIRC=velcirc, VELKIN=velkin, RING=5, /plot, /verbose 
;
;    3) 
;    Run kinemetry on a velocity dispersion map (given by array Sigma). 
;    Desired outputs are position angle and flattening of the ellipses,
;    and harmonic terms: a0, a1, b1, a2, b2, a3, b3, a4, b4:
;
;    KINEMETRY, Xbin, Ybin, Sigma, rad, pa, q, cf, NTRM=10, \EVEN, $
;               /plot, /verbose
;
;REVISION HISTORY:
;   V1.0 - Written by Davor Krajnovic and Michele Cappellari (March, 2005)
;   V2.0 - First released version, Davor Krajnovic, Oxford, 7.12.2005. 
;   V2.1 - Add keyword RADIUS, DK, Oxford, 27.01.2006. 
;   v2.2 - Changed definition of PAQ keyword, DK, Oxford, 02.02.2006.
;   v2.3 - Corrected bug: keyowrd \even was not pass correctly to MPFIT. 
;	   Thanks to Roland Jesseit for prompting. Minor homogonisation 
;          of the code. DK, Oxford, 06.02.2006
;#######################################################################
PRO sauron_colormap
;
; The official SAURON colormap: type "sauron_colormap" to load the colormap.
;
; Michele Cappellari & Eric Emsellem, Leiden, 10 July 2001
;
; Start with these 7 equally spaced coordinates, then add 4 additional points
; x = findgen(7)*255/6. + 1
; 1.0  43.5  86.0  128.5  171.0  213.5  256.0
;
x = [1.0, 43.5, 86.0, 86.0+20, 128.5-10, 128.5, 128.5+10, 171.0-20, 171.0, 213.5, 256.0]

red =   [0.01, 0.0, 0.4,  0.5, 0.3, 0.0, 0.7, 1.0, 1.0,  1.0, 0.9]
green = [0.01, 0.0, 0.85, 1.0, 1.0, 0.9, 1.0, 1.0, 0.85, 0.0, 0.9]
blue =  [0.01, 1.0, 1.0,  1.0, 0.7, 0.0, 0.0, 0.0, 0.0,  0.0, 0.9]

xnew = FINDGEN(256)+1

redV = INTERPOL(red, x, xnew)
greenV = INTERPOL(green, x, xnew)
blueV = INTERPOL(blue, x, xnew)

TVLCT, redV*255, greenV*255, blueV*255 ; load the SAURON colormap

END
;----------------------------------------------------------------------------
FUNCTION range, x1, x2, n
;
; RANGE(x1,x2) = x1,x1+1,...,x2. In this case x1, x2 should be integers.
; RANGE(x1,x2,n) = x1,x1+dx,...,x2 with N integer the result has length N.
; The result will have the type of x1, but if three parameters are used
; the result will be at least of type float.
;
; Michele cappellari, Leiden, 16 October 2001
;
compile_opt idl2
on_error, 2

t = SIZE(x1,/TYPE)
CASE n_params() of
    2: IF x1 LT x2 THEN $
            v = x1 + INDGEN(x2-x1+1, TYPE=t) $
        ELSE $
            v = x1 - INDGEN(x1-x2+1, TYPE=t)
    3: BEGIN
            v = x1 + (x2-x1)/(n-1.0)*INDGEN(n, TYPE=t)
    END
    ELSE: message, '2 or 3 parameters are needed'
ENDCASE

RETURN, v
END
;----------------------------------------------------------------------
PRO kinem_trigrid_irregular, xbin, ybin, moment, xnew, ynew, momNew, _EXTRA=extra
compile_opt idl2
;
; Given a set of irregular coordinates and corresponding values
; output coordinates, which are also irregular, it gives an interpolated values 
; at the output positions.
; The interpolation is done using TRIGRID and any option can be passed
; to this routine via the _EXTRA mechanism.
;
; V1.0: Michele Cappellari, Leiden, 17 December 2004
; V1.1: small cosmetic changes, Davor Krajnovic, Edinburgh, 24.11.2005.
;----------
nbin = N_ELEMENTS(xBin)
IF nbin NE N_ELEMENTS(ybin) or nbin NE N_ELEMENTS(moment) THEN $
    message, 'XBIN, YBIN and VELBIN must have the same size'
nNew = N_ELEMENTS(xNew)
IF nNew NE n_elements(yNew) THEN message, 'XNEW and YNEW must have the same size'

TRIANGULATE, xbin, ybin, tr
momNew = FLTARR(nNew, /NOZERO)
FOR j=0,nNew-1 DO momNew[j] = (TRIGRID(xBin, yBin, moment, tr, $
    XOUT=xNew[j]+[0,1], YOUT=yNew[j]+[0,1], _EXTRA=extra))[0]

END
;----------------------------------------------------------------------
FUNCTION kinem_SVD_Solve, a, b, ERRORS=errors
compile_opt idl2
;
; Solve a linear system using the method
; of the Singular Values Decomposition.
; IF error keyword set, a and b arrays are 
; considered to be properly divided by 
; Gaussian errors. Uncertainties of the returned
; parameters are calculated by computing the 
; covariance matrix and taking the SQRT of 
; the diagonal elements (following svdvar routine
; from Numerical Recepies, Press et al. 1992.). 
;
; This function returns coefficients of the linear 
; system. If /ERRORS is used then first N_ELEMENTS(w)
; values are the coefficients, while the second 
; N_ELEMENTS(w) values are the corresponding 
; uncertainties. 
;
; Written by Davor Krajnovic, Oxford, 25.05.2005, 
; (expanding on original routine by M.Cappellari)
; 

sa = SIZE(a)
sb = SIZE(b)
IF (sa[1] NE sb[1] or sa[0] NE 2 or sb[0] NE 1) THEN $
    message, 'SVD_Solve Error: incompatible dimensions'

SVDC, a, w, u, v, /COLUMN, /DOUBLE
small = WHERE(w LT MAX(w)*1d-13,m)
IF m NE 0 THEN w[small] = 0d
weights=SVSOL(u, w, v, b, /COLUMN, /DOUBLE)


IF keyword_set(errors) THEN BEGIN
	np = N_ELEMENTS(weights)
	wi=dblarr(np) & er = dblarr(np)
	cvm=dblarr(np,np)	; covariance matrix

	FOR i=0, np-1 DO BEGIN 
		FOR j=0, i DO BEGIN
			sum=0.
			FOR k=0, N_ELEMENTS(weights)-1 DO BEGIN
				IF w[k] ne 0. THEN wi[k] = 1./(w[k]*w[k])
				sum = sum + v[i,k]*v[j,k]*wi[k]
			ENDFOR
			cvm[i,j]=sum
			cvm[j,i]=sum
		ENDFOR
	ENDFOR

	FOR i=0, N_elements(weights)-1 DO er[i]=SQRT(cvm[i,i])

	param=[weights, er]
	RETURN, param
ENDIF

RETURN, weights

END
;----------------------------------------------------------------------
FUNCTION kinem_fit_trig_series, x, y, yfit, NTERMS=nterms, COEFF=coeff, $
			        ERR=err, ALL=all, EVEN=even, VSYS=vsys
compile_opt idl2
;
; v1.0 M. Cappellari, Leiden, 15 December 2004
; v1.1 D. Krajnovic, Paranal, 09.03.2004, adapted to be used for
;      reconstruction of kin.moment using arbitrary number of terms
; v1.2 expanded for error treatment, DK, Oxford, 26.05.2005
; v1.3 add vsys keyword and clean comments, DK, Edinburgh, 22.05.2005.
; v1.4 add even keyword, DK, Oxfrod, 06.02.2006.
;
;
; The harmonic series below can include 
;   1) odd J values:
;       1, [sin(x), cos(x)], [sin(3x), cos(3x)],...
;      The even trigonometric terms are zero in the point-symmetric case.
;      NTERMS=4 for all odd terms up to [sin(3x), cos(3x)] and
;      NTERMS=6 for all odd terms up to [sin(5x), cos(5x)]
;   2) even J values
;      1, [sin(x), cos(x)], [sin(2x), cos(2x)],...
;      This is used for a complete harmonic analysis, or to extract
;      EVEN terms. 
;      NTERMS=10 for all terms up to [sin(5x), cos(5x)]
;   3) odd or even, but without zero-th term (or systemic velocity)
;----------------------------------------------------------

arr = DBLARR(n_elements(x),nterms+1)
;
; if keyword vsys set, do not extract zero-th term
;
IF KEYWORD_SET(vsys) THEN arr[*,0] = 0d $
		     ELSE arr[*,0] = 1d


IF keyword_set(all) or keyword_set(even) THEN $ 
	FOR j=1,nterms,2 DO arr[*,[j,j+1]] = [sin(((j+1)/2)*x),cos(((j+1)/2)*x)] $
		    ELSE $
	FOR j=1,nterms,2 DO arr[*,[j,j+1]] = [sin(j*x),cos(j*x)]


;
; divide arrays with corresponding errors
;
IF keyword_set(err) THEN BEGIN
	y1 = y/err
	brr=arr*0.
	FOR j=0, Nterms DO brr[*,j] = arr[*,j]/err
	IF N_ELEMENTS(coeff) EQ 0 THEN coeff = KINEM_SVD_SOLVE(brr,y1, /ERRORS)
	koeff = coeff[0:nterms]
ENDIF ELSE BEGIN
	IF N_ELEMENTS(coeff) EQ 0 THEN coeff = KINEM_SVD_SOLVE(arr,y1)
	koeff = coeff
ENDELSE
yfit = arr # koeff

RETURN, coeff
END
;----------------------------------------------------------------------
FUNCTION kinem_fitfunc_ellipse, p, $
	COEFF=coeff, NTERMS=nterms, XBAR=xbar, YBAR=ybar, MOMENT=moment, $
	RAD=r, XELL=xEll, YELL=yEll, MOMELL=momEll, THETA=theta, $
	MOMFIT=momFit, MPFIT=mpf, ELEM=elem, W=w, ERROR=error, $
	ER_MOMELL=er_momEll, ALL=all, EVEN=even, VSYS=vsys
compile_opt idl2

;
; construction of elliptical coordinates on which kin.moment is
; interpolated; expansion of kin.profile in harmonic series;
; used by both 'brute force' and MPFIT minimisation
ang = p[0]/!RADEG
mi=(360./(180./10))
theta = range(0.0,2.0*!pi, mi>10*r<100) 
	
x = r*COS(theta)
y = r*SIN(theta)*p[1]
xEll = x*COS(ang) - y*SIN(ang)
yEll = x*SIN(ang) + y*COS(ang)

KINEM_TRIGRID_IRREGULAR, xBar, yBar, moment, xEll, yEll, momEll, MISSING=12345678
w = WHERE(momEll NE 12345678, elem)
IF elem eq 0 THEN RETURN, 1e30
	
er_momEll=dblarr(N_ELEMENTS(w))
xxe=xell[w] & yye=yell[w]

;
; Find the closest bin for the interpolated position 
; (elliptical coordinates) and asign the error of 
; the bin to that position 
;
FOR i=0, N_ELEMENTS(w)-1 DO BEGIN
	dist2 = (xbar-xxe[i])^2 + (ybar-yye[i])^2
	wd = WHERE(dist2 eq MIN(dist2))	
	er_momEll[i]=error[wd[0]]
ENDFOR

;
; Calculation of coeffs. with options for ALL, EVEN or VSYS
;
coeff = KINEM_FIT_TRIG_SERIES(theta[w], momEll[w], momFit, NTERMS=nterms,$
				      ERR=er_momEll, ALL=all, EVEN=even, VSYS=vsys)

		
IF KEYWORD_SET(even) THEN BEGIN
	; 
	; Following eq.(1) of Jedrzejewski (1987), it tries to 
	; minimize a1,a2,a2,b2, which indicate incorrect PA and
	; flattening of the trial ellipse (terms are defined
	; differently than in 'odd' case	
	;
	IF KEYWORD_SET(mpf) THEN RETURN, coeff[[1,2,3,4]] $
			     ELSE return, TOTAL(coeff[[1,2,3,4]]^2)


ENDIF ELSE BEGIN
	;
	; Following eq.(1) of Jedrzejewski (1987), but for odd kinematic
	; moments, it tries to minimize a1,a3,b3, which indicate
	; incorrect PA and flattening of the trial ellipse
	;
	IF KEYWORD_SET(mpf) THEN RETURN, coeff[[1,3,4]] $
			    ELSE RETURN, TOTAL(coeff[[1,3,4]]^2)

ENDELSE

END
;----------------------------------------------------------------------
;        MAIN PROGRAM
;----------------------------------------------------------------------
PRO kinemetry, xbin, ybin, moment, rad, pa, q, cf, $
	NTRM=NTRM, ERROR=error, SCALE=scale, NAME=name, $
	PAQ=paq, NPA=npa, NQ=nq, RANGEQ=rangeq, ALL=all, EVEN=even,$
	VSYS=VSYS, VELCIRC=velCirc, VELKIN=velkin, ER_CF=er_cf, $
	ER_PA=er_pa, ER_q=er_q, XELLIP=xellip, YELLIP=yellip, $
	RING=ring, RADIUS=radius, PLOT=plot, VERBOSE=verbose

compile_opt idl2

IF N_ELEMENTS(scale) EQ 0 THEN scale = 0.8 ; pixel size of SAURON maps
IF N_ELEMENTS(ntrm) EQ 0 THEN NTRM=6	   ; number of terms in F.expansion

;
; if errors not specified assume a unit error for all bins
;
IF N_ELEMENTS(error) EQ 0 THEN error = 1. + 0*moment
IF KEYWORD_SET(even) THEN even=1 ELSE even = 0
IF KEYWORD_SET(all) THEN all=1 ELSE all = 0

;
; change to pixel scale
;
xbar=xbin/scale
ybar=ybin/scale

;
; setting radii
;
IF KEYWORD_SET(RADIUS) THEN BEGIN
	rad = radius 
	nrad = N_ELEMENTS(rad)
ENDIF ELSE BEGIN
	nrad = 100 ; Large value: it will be truncated anyway
	pix = FINDGEN(nrad)
	rad = pix + (1.1^pix)
ENDELSE

;
; The central pixel is left unchanged in the reconstruction.
; Shifting of the first radius in case of a central hole.
;
IF KEYWORD_SET(RING) THEN BEGIN
	rad = ring/scale + rad 
	xellip = 0.
	yellip = 0.
	vv = 0.
	vrec = 0.
ENDIF ELSE BEGIN
	mini=MIN(sqrt(xbar^2+ybar^2))
	ww = WHERE(sqrt(xbar^2+ybar^2) eq mini)
	xellip = xbin[ww]
	yellip = ybin[ww]
	vv = moment[ww]
	vrec = moment[ww]
ENDELSE

;
; Initialised vectors of results
;
pa = FLTARR(nrad)
q = pa
cf = FLTARR(nrad,ntrm+1)
er_cf = cf
er_pa = pa
er_q = q

;
; Initialises parameters for MPFIT
;
parinfo = REPLICATE({step:1.0,limits:[0d,0d],limited:[1,1]}, 2)
parinfo[0].limits = [-95d,95d]  ; PA limits in degrees
IF keyword_set(RANGEQ) THEN parinfo[1].limits = [RANGEQ[0],RANGEQ[1]] $
                       ELSE parinfo[1].limits = [0.2d,1d]   ; q limits
parinfo[0].step = 0.5d  ; Step in degrees (of the order of the expected accuracy)
parinfo[1].step = 0.01d ; q

;
; Set grids for global minimization
;
IF N_ELEMENTS(npa) EQ 0 THEN npa = 21
IF N_ELEMENTS(nq) EQ 0 THEN nq = 21
IF N_ELEMENTS(rangeQ) EQ 0 THEN rangeQ=[0.2,1.0] 

pa_grid = RANGE(-90.0,90.0,npa) # REPLICATE(1,nq)
q_grid = REPLICATE(1,npa) # RANGE(rangeQ[0],rangeQ[1],nq)
chi2_grid = FLTARR(nPa,nq)
k = 0
;
; loop over radii
;
FOR i=0, nrad-1 DO BEGIN
    ;
    ; check if PA and q are set constant for the whole map (PAQ)
    ;
    IF N_ELEMENTS(PAQ) EQ 0 THEN BEGIN
        ;
        ; Perform brute-force *global* minimization of chi2 on a regular grid. 
        ; This is needed as the problem can have multiple minima.
        ;
        FOR j=0,npa*nq-1 DO $
	      chi2_grid[j] = kinem_fitfunc_ellipse([pa_grid[j], q_grid[j]], $
        	XBAR=xbar, YBAR=ybar, MOMENT=moment, NTERMS=4, RAD=rad[i], $
		ERROR=error, EVEN=even); $
        tmp = MIN(chi2_grid,jmin)
        ;
        ; Perform least-squares minimization of the a1,a3,b3 coefficients
        ; starting from the best values of the global minimization.
        ;
        par = [pa_grid[jmin], q_grid[jmin]]
        fa = {MOMENT:moment, XBAR:xbar, YBAR:ybar, NTERMS:4, RAD:rad[i], MPF:1, $
	      ERROR:error, EVEN:even}
        sol = MPFIT('kinem_fitfunc_ellipse', par, FUNCTARGS=fa, FTOL=1e-3, $
             PARINFO=parinfo, NFEV=ncalls, BESTNORM=fmin, PERROR=err_sol, /QUIET)
        PA_min = sol[0]
        q_min = sol[1]

	;
	; errors to ellipse parameters
	;
	er_PA_min = err_sol[0]	
	er_q_min  = err_sol[1] 


        IF KEYWORD_SET(verbose) THEN print, 'done ', i, '-th radius: ', rad[i]*scale

    ENDIF ELSE BEGIN ; PA and Q are set to fixed values: skip all minimization
	
	;
	; check if PAQ is an array of (nrad*2) values or has just 2 values
	;
	IF N_ELEMENTS(PAQ) gt 2 THEN BEGIN
        	pa_min = paq[k]
        	q_min = paq[k+1]
		k = k + 2
	ENDIF ELSE BEGIN
        	pa_min = paq[0]
        	q_min = paq[1]
	ENDELSE
	er_PA_min = 0.
	er_q_min  = 0. 

    ENDELSE

    ;
    ; Final harmonic expansion along the best fitting ellipse
    ;
    tmp = kinem_fitfunc_ellipse([pa_min, q_min], $
        XBAR=xbar, YBAR=ybar, MOMENT=moment, NTERMS=ntrm, RAD=rad[i], COEFF=coeff, $
        XELL=xEll, YELL=yEll, MOMELL=momEll, MOMFIT=recon, THETA=theta, W=w, $
	ERROR=error, ER_MOMELL=er_momEll, ALL=all, EVEN=even, VSYS=vsys)

    ;	
    ; errors of harmonic coeffs
    ;
    er_coeff = coeff[ntrm+1 : 2*ntrm+1]
    coeff = coeff[0:ntrm]

    ;
    ; Stops the fit when there are less than 3/4 of the pixels sampled
    ; along the best fitting ellipse. This condition can be relaxed when 
    ; PAQ set.
    ;
    IF N_ELEMENTS(w) LT N_ELEMENTS(xEll)*0.75 THEN break


    ;
    ; Assigning of outputs
    ;
    pa[i] = pa_min
    q[i] = q_min
    cf[i,*] = coeff

    er_cf[i,*] = er_coeff
    er_pa[i] = er_pa_min
    er_q[i] = er_q_min

    IF KEYWORD_SET(vsys) THEN BEGIN
			IF vsys ne 1 THEN cf[*,0] = Vsys ELSE cf[*,0] = 0.
			er_cf[*,0] = 0.
    ENDIF

    xellip = [xellip,xell[w]]
    yellip = [yellip,yell[w]]
    IF KEYWORD_SET(even) THEN vv = [vv, coeff[0] + 0.*coeff[2]*cos(theta[w])] $
			 ELSE vv = [vv, coeff[0] + coeff[2]*cos(theta[w])]
    vrec = [vrec, recon]

    ;
    ; Optional plotting
    ;
    IF KEYWORD_SET(plot) THEN BEGIN
        !P.MULTI=[0,1,4]
        sauron_colormap

        IF KEYWORD_SET(EVEN) THEN BEGIN
	        mn = MIN(SIGRANGE(moment),MAX=mx)
		levels = mn + (mx-mn) * FINDGEN(64)/(64-1)
		CONTOUR, moment>mn<mx, xbin, ybin, /FILL, /ISO, LEVELS=levels, $
		    /XSTYLE, /YSTYLE, /IRREGULAR, XTITLE='arcsec', YTITLE='arcsec'
        ENDIF ELSE BEGIN
	        mn = MIN(SIGRANGE(moment-median(moment)),MAX=mx)
                mx = mx < abs(mn)
		levels = -mx + (mx-(-mx)) * FINDGEN(64)/(64-1)
		CONTOUR, (moment-median(moment))>(-mx)<mx, xbin, ybin, /FILL, /ISO, LEVELS=levels, $
		    /XSTYLE, /YSTYLE, /IRREGULAR, XTITLE='arcsec', YTITLE='arcsec'
	ENDELSE
        PLOTS, xEll*scale, yEll*scale, COLOR=0
        OPLOT, [-Rad[i]*scale,Rad[i]*scale]*COS(PA[i]/!RADEG), [-rad[i]*scale,rad[i]*scale]*SIN(pa[i]/!RADEG), COLOR=0

        LOADCT, 4, /silent
        CONTOUR, ALOG(chi2_grid), pa_grid, q_grid, NLEVELS=15, XTITLE='PA (deg)', YTITLE='q', /FILL
        PLOTS, pa_min, q_min, PSYM=4, COLOR=125, THICK=3
        PLOTS, pa_grid, q_grid, PSYM=3, THICK=2

        LOADCT, 12, /silent
        !Y.MARGIN=[0,0] ; show plots with shared X axis
        !Y.OMARGIN=[5,10] ; allow for space for the axis labels

	IF KEYWORD_SET(EVEN) THEN BEGIN
	  IF N_ELEMENTS(error) gt 1 THEN $
	        PLOTERROR, theta[w]*!RADEG, momEll[w], er_momEll, PSYM=-1, TITLE=rad[i]*scale, XTICKFORMAT='(A1)', YTITLE='V' $
	  ELSE $
        	PLOT, theta[w]*!RADEG, momEll[w], PSYM=-1, TITLE=rad[i], XTICKFORMAT='(A1)', YTITLE='V'
          OPLOT, theta*!RADEG, coeff[0]+coeff[2]*COS(theta), PSYM=-4, COLOR=200

          PLOT, theta*!RADEG, SIGRANGE(momEll-coeff[0]), PSYM=-4, XTITLE='!7h!X!N', YTITLE='residuals', YRANGE=[-20,20]
          OPLOT, theta*!RADEG, coeff[1]*SIN(theta) + coeff[2]*COS(theta) + coeff[3]*SIN(2*theta) + coeff[4]*COS(2*theta), PSYM=-1, COLOR=40
	  IF NTRM ge 8 THEN $
	          OPLOT, theta*!RADEG, coeff[1]*SIN(theta) + coeff[2]*COS(theta) + coeff[3]*SIN(2*theta) + coeff[4]*COS(2*theta) $
                       + coeff[7]*sin(4*theta) + coeff[8]*cos(4*theta), psym=-1, COLOR=200	

	ENDIF ELSE BEGIN

	  IF KEYWORD_SET(all) THEN BEGIN
		a1 = coeff[1] & b1 = coeff[2] & a3 = coeff[5] & b3 = coeff[6]
	  ENDIF ELSE BEGIN
		a1 = coeff[1] & b1 = coeff[2] & a3 = coeff[3] & b3 = coeff[4]
	  ENDELSE		

	  IF N_ELEMENTS(error) gt 1 THEN $
	        PLOTERROR, theta[w]*!RADEG, momEll[w], er_momEll, PSYM=-1, TITLE=rad[i]*scale, XTICKFORMAT='(A1)', YTITLE='V' $
	  ELSE $
        	PLOT, theta[w]*!RADEG, momEll[w], PSYM=-1, TITLE=rad[i], XTICKFORMAT='(A1)', YTITLE='V'
          OPLOT, theta*!RADEG, coeff[0]+b1*COS(theta), PSYM=-4, COLOR=200

          PLOT, theta*!RADEG, SIGRANGE(momEll-coeff[0]-b1*COS(theta)), PSYM=-4, XTITLE='!7h!X!N', YTITLE='residuals', YRANGE=[-20,20]
          OPLOT, theta*!RADEG, a1*SIN(theta) + a3*SIN(3*theta) + b3*COS(3*theta), PSYM=-1, COLOR=40
	  IF NTRM ge 6 THEN BEGIN
		IF KEYWORD_SET(all) and NTRM ge 9 THEN BEGIN
			a5 = coeff[9] & b5 = coeff[10]
		ENDIF ELSE BEGIN
	 		a5 = coeff[5] & b5 = coeff[6]
	  	ENDELSE		
          	OPLOT, theta*!RADEG, a1*SIN(theta) + a3*SIN(3*theta) + b3*COS(3*theta) $
                     + a5*sin(5*theta) + b5*cos(5*theta), psym=-1, COLOR=200	
	  ENDIF
        ENDELSE

        !Y.MARGIN=[4,2] ; back to default values
        !Y.OMARGIN=[0,0]
        !P.MULTI=0
    ENDIF
ENDFOR  ; end of loop over radii



;
; Final outputs (back to physical scale (arcsec))
;
xellip=xellip*scale
yellip=yellip*scale
rad = rad[0:i-1]*scale
pa = pa[0:i-1]
q = q[0:i-1]
cf = cf[0:i-1,*]

IF N_ELEMENTS(error) gt 1 THEN BEGIN
	er_cf = er_cf[0:i-1,*]
	er_pa = er_pa[0:i-1]
	er_q  = er_q[0:i-1]
ENDIF

;
; PA correction
;
IF N_ELEMENTS(PAQ) EQ 0 THEN BEGIN
	PA=270 + PA 		     	 ; PA measured East of North
	wb=WHERE(cf[*,2] lt 0., red) 	 ; PA measured from receding side of galaxy (positive vel)
	IF red ne 0 AND even eq 0 THEN $ ; no 180 deg correction for even moments
		   PA[wb]=pa[wb]-180
ENDIF
;
; reconstruction of kinematic moment map
;
KINEM_TRIGRID_IRREGULAR, xellip, yellip, vv, xbin, ybin, velCirc, MISSING=12345
KINEM_TRIGRID_IRREGULAR, xellip, yellip, vrec, xbin, ybin, velKin

;
; Reconstruction cosmetic lines: un-commenting these lines will assign 
; original values to bins that are outside analysed region. *Not* recommended
; except for presentation purposes, since there is no science
; associated and may introduce confusion when interpreting results.
;
;wf = WHERE(velCirc eq 12345, elem)
;IF elem ne 0 THEN velCirc[wf] = velBin[wf] ; Leave unchanged values not fitted
;IF elem ne 0 THEN velKin[wf] = velBin[wf]


END
;----------------------------------------------------------------------
