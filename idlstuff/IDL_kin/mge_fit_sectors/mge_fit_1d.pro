;######################################################################
;
; Copyright (C) 1999-2001, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; For details on the method see:
;   Cappellari M., 2002, MNRAS, 333, 400
;
; Updated versions of the software are available from my web page
; http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
; If you have found this software useful for your
; research, I would appreciate an acknowledgment to use of
; `the MGE fitting software developed by Cappellari (2002)'.
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;######################################################################
;+
; NAME:
;     MGE_FIT_1D
;
; AUTHOR:
;       Michele Cappellari, Astrophysics Sub-department, University of Oxford, UK
;
; PURPOSE:
;       Perform a Multi Gaussian Expansion fit to a one-dimensional profile.
;
; EXPLANATION:
;       Further information on MGE_FIT_1D is available in
;       Cappellari M., 2002, MNRAS, 333, 400
;
; CALLING SEQUENCE:
;       MGE_FIT_1D, X, Y, NGAUSS=ngauss, OUTER_SLOPE=outer_slope, RBOUNDS=rbounds
;
; INPUTS:
;       X = Vector of the logarithmically sampled (Important !)
;           abscissa for the profile that has to be fitted.
;       Y = Ordinate of the profile evaluated at the abscissa X.
;
; OUTPUTS:
;       No output parameters. The results are printed on the screen
;       and passed with the optional keyword SOL
;
; OPTIONAL INPUT KEYWORDS:
;       NGAUSS - Number of Gaussian on want to fit. Typical values are 10-20.
;       OUTER_SLOPE - This scalar forces the surface brightness profile of
;               the MGE model to decrease at least as fast as R^(-OUTER_SLOPE)
;               at the largest measured radius (Default: OUTER_SLOPE=2).
;       RBOUNDS - Two elements vector giving the minimum and maximum sigma
;               allowed in the MGE fit. This is in the same units of X.
;
; EXAMPLE:
;       The sequence of commands below was used to generate the
;       one-dimensional MGE fit of Fig.3 in Cappellari (2002).
;
;           n = 200
;           rminLog = ALOG10(0.01)
;           rmaxLog = ALOG10(300)
;           R = 10^((rmaxLog-rminLog)*FINDGEN(n)/(n-1)+rminLog) ; Radius
;           rho = (1 + R)^(-4) ; The profile is logarithmically sampled!
;           mge_fit_1d, R, rho, NGAUSS=16, SOL=sol
;      
;        In the common case in which rho represents an intrinsic density
;        in Msun/pc^3 and R is in pc, the output keyword sol[0,*] will 
;        contain the surface density of the projected Gaussians in Msun/pc^2. 
;        This is already the required input format for the JAM modelling 
;        routines here http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
; OPTIONAL OUTPUT KEYWORDS:
;       SOL - Output keyword containing a 2xNgauss array with the
;               best-fitting solution:
;               1) sol[0,*] = TotalCounts, of the Gaussians components.
;                   The relation TotalCounts = Height*SQRT(2*!PI)*Sigma
;                   can be used compute the Gaussian central surface
;                   brightness (Height).
;                   IMPORTANT: TotalCounts is defined as the integral under a
;                   1D Gaussian curve and not the one of a 2D Gaussian surface.
;               2) sol[1,*] = sigma, is the dispersion of the best-fitting
;                   Gaussians in units of X.
; PROCEDURES USED:
;       The following procedures are contained in the main MGE_FIT_1D program.
;           MGE_FIT_1D_PRINT   -- Show intermediate results while fitting
;           FITFUNC_MGE_1D     -- This the function that is optimized during the fit.
;
;       Other IDL routines needed:
;           BVLS  -- Michele Cappellari porting of Lawson & Hanson generalized NNLS
;                    http://www.strw.leidenuniv.nl/~mcappell/idl/
;           MGE_MPFIT -- Craig Markwardt porting of Levenberg-Marquardt MINPACK-1
;
; MODIFICATION HISTORY:
;       V1.0 Michele Cappellari, Padova, February 2000
;       V1.1 Minor revisions, MC, Leiden, June 2001
;       V1.2 Updated documentation, MC, 18 October 2001
;       V1.3 Added compilation options and explicit stepsize in parinfo
;            structure of MPFIT, MC, Leiden 20 May 2002
;       V1.31 Added ERRMSG keyword to MPFIT call. Leiden, 13 October 2002, MC
;       V1.32: Use N_ELEMENTS instead of KEYWORD_SET to test
;           non-logical keywords. Leiden, 9 May 2003, MC
;       V1.33: Use updated calling sequence for BVLS. Leiden, 20 March 2004, MC
;       V1.4: Force the surface brightness of the MGE model to decrease at
;           least as R^-2 at the largest measured radius. Leiden, 8 May 2004, MC
;       V1.41: Make sure this routine uses the Nov/2000 version of Craig Markwardt
;           MPFIT which was renamed MGE_MPFIT to prevent potential conflicts with
;           more recent versions of the same routine. Vicenza, 23 August 2004, MC.
;       V1.42: Allow forcing the outer slope of the surface brightness profile of
;           the MGE model to decrease at least as R^-n at the largest measured
;           radius (cfr. version 1.4). Leiden, 23 October 2004, MC
;       V1.43: Changed axes labels in plots. Leiden, 18 October 2005, MC
;-
;----------------------------------------------------------------------------
PRO mge_fit_1d_print, MYFUNCT, logsigma, iter, fnorm, $
    FUNCTARGS=fcnargs, PARINFO=parinfo, QUIET=quiet, DOF=dof
;
; This is a plotting routine that is called every NPRINT iteration of MPFIT
;
COMPILE_OPT IDL2, HIDDEN

COMMON mge_fit_1d, x, y, gauss, soluz, yfit, err

PRINT, 'Iteration: ', STRTRIM(iter,2), ';  chi2: ', STRTRIM(fnorm,2)
ngauss = N_ELEMENTS(logsigma)

!P.MULTI=[0,1,2]
!Y.MARGIN=[0,0] ; show plots with shared X axis
!Y.OMARGIN=[5,5] ; allow for space for the axis labels

PLOT, x, y, /XLOG, /YLOG, YTITLE='counts', $
    /XSTYLE, YSTYLE=2, XTICKFORMAT='(A1)', TICKLEN=0.04
OPLOT, x, yfit, LINE=2
FOR j=0,ngauss-1 DO OPLOT, x, gauss[*,j]*soluz[j]

PLOT, x, err*100, /XLOG, /YSTYLE, YTITLE='error (%)', XTITLE='arcsec', $
    /XSTYLE, XTICKNAME=!x.tickname, YRANGE=[-20,19.9], TICKLEN=0.04, YMINOR=5
FOR j=0,ngauss-1 DO PLOTS, 10^logsigma[[j,j]], [-0.02,-0.01], THICK=2

!P.MULTI=0
!Y.MARGIN=[4,2] ; back to default values
!Y.OMARGIN=[0,0]

END
;----------------------------------------------------------------------------
FUNCTION fitfunc_mge_1d, logsigma
;
COMPILE_OPT IDL2, HIDDEN

COMMON mge_fit_1d

ngauss = N_ELEMENTS(logsigma)
npoints = N_ELEMENTS(x)
A = DBLARR(npoints,ngauss)
gauss = A
sigma2 = 10d^(2d*logsigma)

FOR j=0,ngauss-1 DO BEGIN
    gauss[*,j] = EXP(-x^2/(2d*sigma2[j])) / SQRT(2d*!DPI*sigma2[j])
    A[*,j] = gauss[*,j]/y
ENDFOR
B = REPLICATE(1d,npoints)

BND = DBLARR(2,ngauss)
BND[0,*] = 0d
BND[1,*] = (MACHAR()).XMAX
BVLS, A, B, BND, soluz, IERR=ierr, ITMAX=ngauss*30
IF ierr NE 0 THEN MESSAGE, 'BVLS Error n. ' + STRTRIM(ierr,2) + ': try to decrease NGAUSS'

yfit = gauss # soluz
err = 1d - yfit/y

RETURN, FLOAT(err)
END
;----------------------------------------------------------------------------
PRO mge_fit_1d, x1, y1, NGAUSS=ngauss, RBOUNDS=rbounds, $
    SOL=sol, OUTER_SLOPE=outer_slope
;
; This is the main routine that has to be called from outside
;
COMPILE_OPT IDL2
ON_ERROR, 2

COMMON mge_fit_1d

x = x1 ; load parameters into the COMMON block
y = y1

IF N_ELEMENTS(outer_slope) EQ 0 THEN outer_slope = 2.0 $
    ELSE outer_slope = outer_slope > 1.0 < 4.0
IF N_ELEMENTS(ngauss) EQ 0 THEN ngauss=12
IF N_ELEMENTS(rbounds) EQ 0 THEN BEGIN
    rmin = MIN(x, MAX=rmax)
    rminLog = ALOG10(rmin)
    rmaxLog = ALOG10(rmax/sqrt(outer_slope))
ENDIF ELSE BEGIN
    rminlog = ALOG10(rbounds[0])
    rmaxlog = ALOG10(rbounds[1])
ENDELSE
logsigma = ((rmaxLog-rminLog)*(0.5+FINDGEN(ngauss))/ngauss + rminLog)
parinfo = REPLICATE({limited:[1,1], limits:[rminLog,rmaxLog], step:0.01}, ngauss)

t = SYSTIME(1)
p = mge_mpfit('fitfunc_mge_1d', logsigma, ITERPROC='mge_fit_1d_print', $
    NPRINT=5, PARINFO=parinfo, BESTNORM=bestnorm, ERRMSG=errmsg)
IF errmsg NE '' THEN MESSAGE, errmsg
;
; Print the results for the nonzero Gaussians
;
w = WHERE(soluz NE 0, m)
logSigma = p[w]
soluz = soluz[w]
j = SORT(logSigma) ; Sort by increasing sigma

sol = FLTARR(2,m)
sol[0,*] = soluz[j]
sol[1,*] = 10^logSigma[j]

PRINT, '############################################'
PRINT, ' Computation time: ', STRTRIM(SYSTIME(1)-t,2), ' seconds'
PRINT, 'Nonzero Gaussians: ', STRTRIM(m,2)
PRINT, ' Unused Gaussians: ', STRTRIM(ngauss-m,2)
PRINT, ' Chi2: ', STRTRIM(bestnorm,2)
PRINT, ' STDEV: ', STRTRIM(STDDEV(err),2)
PRINT, ' MEANABSDEV: ', STRTRIM(MEANABSDEV(err),2)
PRINT, ' Total_Counts     Sigma'
PRINT, '############################################'
PRINT, sol
PRINT, '############################################'

END
;----------------------------------------------------------------------------
