;######################################################################
;
; Copyright (C) 1999-2005, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; For details on the method see:
;   Cappellari M., 2002, MNRAS, 333, 400
;
; Updated versions of the software are available from my web page
; http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
; If you have found this software useful for your
; research, I would appreciate an acknowledgment to
; `use of the MGE fitting software developed by Cappellari (2002)'.
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
;       TEST_MGE_FIT
;
; AUTHOR:
;       Michele Cappellari, Astrophysics Sub-department, University of Oxford, UK
;
; PURPOSE:
;       Exercise all the routines in the MGE_FIT_SECTORS package.
;       This procedure is intended to be used as a template to be
;       customized for each particular MGE fitting problem.
;
; EXPLANATION:
;       Further information is available in
;       Cappellari M., 2002, MNRAS, 333, 400
;
; CALLING SEQUENCE:
;       TEST_MGE_FIT
;
; EXAMPLE:
;       The following command will take various minutes to complete,
;       while plotting intermediate results on the screen
;
;           test_mge_fit
;
; PROCEDURES USED:
;       The following procedures are contained in the TEST_MGE_FIT program.
;           FIT_M32     -- Perform an MGE fit of the galaxy M32
;           FIT_NGC4342 -- Perform an MGE fit of the galaxy NGC 4342
;           FIT_NGC4473 -- Perform an MGE fit of the galaxy NGC 4473
;           FIT_DOUBLE_POWERLAW_1D -- Fit a 1D MGE model
;           MGE_FIT_1D_HERNQUIST_MODEL -- Compute the circular velocity of an MGE 1D fit
;           FIT_NGC5831_TWIST -- Perform a twisted MGE fit of the galaxy NGC 5831
;
;       Other routines needed from the MGE_FIT_SECTORS package by
;       Michele Cappellari (http://www.strw.leidenuniv.nl/~mcappell/idl)
;           SECTORS_PHOTOMETRY       -- perform photometry along sectors
;           MGE_FIT_SECTORS          -- do the actual MGE fit
;           MGE_PRINT_CONTOURS       -- plot the results
;           SECTORS_PHOTOMETRY_TWIST -- perform photometry along sectors with point symmetry
;           MGE_FIT_SECTORS_TWIST    -- do the actual MGE fit with possible isophote twist
;           MGE_PRINT_CONTOURS_TWIST -- plot the results with possible isophote twist
;           WFPC2_MGE_FIT            -- do an MGE fit of a WFPC2 image (PC1 + MOSAIC)
;           FIND_GALAXY              -- find galaxy center and position angle
;           MGE_FIT_1D               -- perform a 1D MGE fit
;
;       Other IDL routines needed:
;           BVLS  -- Michele Cappellari porting of Lawson & Hanson generalized NNLS
;                    http://www-astro.physics.ox.ac.uk/~mxc/idl/
;           MPFIT -- Craig Markwardt porting of Levenberg-Marquardt MINPACK-1
;                    http://astrog.physics.wisc.edu/~craigm/idl/
;           JAM package -- The Jeans Anisotropic MGE (JAM) package is required only
;                    if one wants to run the optional example MGE_FIT_1D_HERNQUIST_MODEL
;                    http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
;       Astronomy User's Library routines needed (http://idlastro.gsfc.nasa.gov):
;           FITS_READ
;           DIST_CIRCLE
;
; MODIFICATION HISTORY:
;       V1.0: Michele Cappellari, Leiden, January 2000
;       V2.0: Updated all examples, MC, Leiden July 2001
;       V2.1: Added 1D MGE fit example, MC, Leiden, 10 May 2003
;       V2.11: Removed the un-necessary parameter EPS from the routine
;           SECTORS_PHOTOMETRY_TWIST. MC, Leiden, 1 September 2004
;       V2.12: Replaced LOGRANGE keyword with the new MAGRANGE.
;           MC, Leiden, 1 May 2005
;       V2.2: Included a new example MGE_FIT_1D_HERNQUIST_MODEL.
;           MC, Oxford, 30 November 2008
;-
;----------------------------------------------------------------------------
PRO fit_m32
;
; This procedure reproduces Figures 6-7 in Cappellari (2002).
;
; This example illustrates how to fit multiple images together.
; We model an HST/WFPC2/F814W image of M32 and an I-band ground-based one.
;
scale1 = 0.0455     ; (arcsec/pixel) PC1. This is used as scale and flux reference!
scale2 = 0.549      ; (arcsec/pixel) ING (Peletier 1993)
scaleRatio = scale1/scale2
fluxRatio = 0.9579  ; = flux1/flux2

; Perform photometry on the HST/WFPC2/F814W (I-band) image
; The geometric parameters below were obtained using my FIND_GALAXY program
;
eps1 = 0.2
ang1 = 104.6
xc1 = 314
yc1 = 377
fits_read, 'm32_f814w_pc.fits', img1, h ; F814W
skylev = 4.48
img1 = img1 - skylev
sectors_photometry, img1, eps1, ang1, xc1, yc1, radius1, angle1, cnt1

; Perform photometry on Peletier (1993) ING/I-band image
; The geometric parameters below were obtained using my FIND_GALAXY program
;
eps2 = 0.2
ang2 =  67.1
xc2 = 184
yc2 = 376
fits_read, 'm32_i.fits', img2, h ; I-band
skylev = 32.0 ; sky determined to make outer profile asymptotically power-law
img2 = img2 - skylev
sectors_photometry, img2, eps2, ang2, xc2, yc2, radius2, angle2, cnt2, MINLEVEL=5.0
radius2 = radius2/scaleRatio
cnt2 = cnt2*fluxRatio

; Exclude pixels at small radii (<3") in Peletier's image to avoid
; PSF effects and merges the profiles of the different images.
; The HST image is used as flux and spatial scale reference,
; the ground-based data are simply scaled to the HST units.
;
w = WHERE(radius2 GT 3./scale1)
radius = [radius1,radius2[w]]
angle = [angle1,angle2[w]]
counts = [cnt1,cnt2[w]]

sigmaPSF = 0.8  ; WFPC2/PC1/F814W image (here we simply use one Gaussian PSF)
ngauss = 11
eps = 0.2       ; This is just an initial estimate for the ellipticity

; Do the actual MGE fit. The profiles are saved in idl.ps
;
MGE_fit_sectors, radius, angle, counts, eps, $
    NGAUSS=ngauss, SIGMAPSF=sigmaPSF, SCALE=scale1, SOL=sol, QBOUNDS=[0.3,0.85]

; Plot contours of the HST image
;
MGE_print_contours, img1, ang1, xc1, yc1, sol, $
    FILE='m32_hst.ps', BINNING=4, SCALE=scale1, SIGMAPSF=sigmaPSF

; Scale the solution parameters to the ground-based image
;
ngauss = N_ELEMENTS(sol)/3
j = INDGEN(ngauss)*3
sol[j] = sol[j]*scaleRatio^2/fluxRatio
sol[j+1] = sol[j+1]*scaleRatio ; sigma
sigmaPSF = 1.93/2.35 ; seeing FWHM = 1.0 arcsec for the ground based image
MGE_print_contours, img2, ang2, xc2, yc2, sol, $
    FILE='m32_ground.ps', SCALE=scale2, BINNING=4, MAGRANGE=9.5, SIGMAPSF=sigmaPSF

END
;----------------------------------------------------------------------------
PRO fit_ngc4342
;
; This procedure reproduces Figures 8-9 in Cappellari (2002)
;
; This example illustrates a simple MGE fit to one single HST/WFPC2/F814W image.
;
fits_read, 'ngc4342_f814w_pc.fits', img, h
skylev = 0.55
img = img - skylev
scale = 0.0455

ngauss = 12

; Here we use an accurate four gaussians MGE PSF for
; the HST/WFPC2/F814W filter, taken from Table 3 of
; Cappellari et al. (2002, ApJ, 578, 787)

sigmaPSF = [0.494, 1.44, 4.71, 13.4]      ; In PC1 pixels
normPSF = [0.294, 0.559, 0.0813, 0.0657]  ; total(normPSF)=1

; Here we use FIND_GALAXY directly inside the procedure. Usually you may want
; to experiment with different values of the FRACTION keyword, before adopting
; given values of Eps, Ang, Xc, Yc.

find_galaxy, img, majorAxis, eps, ang, xc, yc, /PLOT

; Perform galaxy photometry

sectors_photometry, img, eps, ang, xc, yc, radius, angle, counts, MINLEVEL=0.2

; Do the actual MGE fit

MGE_fit_sectors, radius, angle, counts, eps, $
    NGAUSS=ngauss, SIGMAPSF=sigmaPSF, NORMPSF=normPSF, SOL=sol, SCALE=scale;, /BULGE

; Print the data-model contours comparison of the whole image

MGE_print_contours, img, ang, xc, yc, sol, $
    FILE='ngc4342.ps', SCALE=scale, MAGRANGE=9, $
    SIGMAPSF=sigmaPSF, NORMPSF=normPSF, BINNING=7

; Print the data-model contours comparison of the central regions

s = SIZE(img)
img = img[xc-s[1]/9:xc+s[1]/9,yc-s[2]/9:yc+s[2]/9]
MGE_print_contours, img, ang, s[1]/9, s[2]/9, sol, $
    FILE='ngc4342_nuclear.ps', SCALE=scale, MAGRANGE=9, $
    SIGMAPSF=sigmaPSF, NORMPSF=normPSF

END
;----------------------------------------------------------------------------
PRO fit_ngc4473
;
; This procedure fits an MGE model to the WFPC2 photometry of NGC4473
;
; These parameters are for the PC1 trimmed image as given by find_galaxy
;
skylevel = 4.17
sigmaPSF = 0.8 ; pixels
eps = 0.42
ang = -79.4
xc = 378
yc = 376
ngauss = 16
minlevel = 0.1

WFPC2_MGE_fit, eps, ang, xc, yc, NGAUSS=ngauss, SIGMAPSF=sigmaPSF, $
    SKYLEVEL=skylevel, MINLEVEL=minlevel, FILENAME='ngc4473_f814w.fits'

END
;----------------------------------------------------------------------------
PRO fit_double_powerlaw_1d
;
; This example reproduces Figure 3 in Cappellari (2002)
; See the next procedure mge_fit_1d_hernquist_model 
; for an example using physical units.

n = 200
rminLog = ALOG10(0.01)
rmaxLog = ALOG10(300)
x = 10d^((rmaxLog-rminLog)*FINDGEN(n)/(n-1)+rminLog)
y = (1d + x)^(-4) ; The profile should be logarithmically sampled!
mge_fit_1d, x, y, NGAUSS=16

END
;----------------------------------------------------------------------------
PRO mge_fit_1d_hernquist_model
;
; This example shows how to use the MGE_FIT_1D to fit an intrinsic density
; profile, for subsequent use in the Jeans Anisotropic MGE (JAM) package. 
; Here we use the density of a Hernquist model 
; (Hernquist 1990, ApJ, 356, 359; hereafter H90).
; 
; To run this example the JAM package must also be downloaded from here 
; http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
; Michele Cappellari, Oxford, 28 November 2008

M = 1d11         ; Total mass of the H90 model in Solar Masses
a = 1d3          ; Break radius of the H90 profile in pc
distance = 16.5d ; Assume Virgo distance in Mpc
pc = distance*!dpi/0.648d ; Constant factor to convert arcsec --> pc

n = 100 ; Number of values to sample the profile for the fit
r = range(a/100,a*10,n,/LOG) ; logarithically spaced radii in pc
rho = M*a/(2d*!dpi*r*(r+a)^3) ; Density in Msun/pc^3 (H90 equation 2)

mge_fit_1d, r, rho, NGAUSS=16, SOL=sol

surf = sol[0,*]     ; Surface density in Msun/pc^2
sigma = sol[1,*]/pc ; Gaussian dispersions in arcsec
qobs = surf*0+1     ; Assume a spherical geometry
inc = 90d
mbh = 0d
rad = range(1d,100,100) ; desired output radii in arcsec

mge_circular_velocity, surf, sigma, qobs, inc, mbh, distance, rad, vcirc

plot, rad, vcirc, XTITLE='R (arcsec)', YTITLE='Vc (km/s)'

; Compare with the analytic H90 result for the circular velocity
;
G = 0.00430237d ; (km/s)^2 pc/Msun [6.674e-11 SI units]
vc = sqrt(G*M*r)/(r+a) ; (H90 equation 16)
loadct, 12
oplot, r/pc, vc, COLOR=200

END
;----------------------------------------------------------------------------
PRO fit_ngc5831_twist
;
; This procedure reproduces Figures 11-12 in Cappellari (2002)
;
; These parameters are given by find_galaxy for the mosaic image
;
skylevel = 13.0
sigmaPSF = 0.4 ; pixels
eps = 0.28
ang = 141.0 ; major axis in the inner regions (gives a starting guess for the PA)
xc = 969
yc = 974
ngauss = 11
minlevel = 1.0
scale = 0.100

filename = 'ngc5831_f702w_mosaic.fits'
fits_read, filename, img, h
badpixelsMosaic = where(img le 0)
img = img - skyLevel
imgClean = img ; save a non-masked image for the final contour plot

; Mask a galaxy and two bright stars
;
s = SIZE(img)
dist_circle, r, s[1:2], 1408.09, 357.749
badpixelsMosaic = [badpixelsMosaic,WHERE(r LT 200)]
dist_circle, r, s[1:2], 539.243, 161.264
badpixelsMosaic = [badpixelsMosaic,WHERE(r LT 50)]
dist_circle, r, s[1:2], 393.028, 713.681
badpixelsMosaic = [badpixelsMosaic,WHERE(r LT 50)]

; Perform galaxy photometry
;
sectors_photometry_twist, img, ang, xc, yc, radius, angle, counts, $
    BADPIXELS=badpixelsMosaic, MINLEVEL=minlevel

; Do the actual MGE fit
;
MGE_fit_sectors_twist, radius, angle, counts, eps, $
    NGAUSS=ngauss, SIGMAPSF=sigmaPSF, SOL=sol, SCALE=scale, $
    QBOUNDS=[0.2,1], PABOUNDS=[-40.0,40.0]

; Print the data-model contours comparison
;
MGE_print_contours_twist, imgClean, ang, xc, yc, sol, $
    FILE='NGC5831_Mosaic.ps', SCALE=scale, BINNING=6, MAGRANGE=10, SIGMAPSF=sigmaPSF

END
;----------------------------------------------------------------------------
PRO test_MGE_fit
;
; This is the main routine to call in succession the MGE fits to
; M32, NGC4342, NGC4473, power-law and NGC5831, and measure the execution time.
; A run of this program takes: 
; - 691 s on a Pentium III, 1.0GHz PC, with IDL 5.4.
; - 270 s on a Pentium M, 1.6GHz PC, with IDL 6.1.
; - 150 s on a Core2 Duo, 2.2GHz PC, with IDL 7.0.
; It was tested with IDL 5.6-7.0 under both Windows and Linux.
;
t = SYSTIME(1)
fit_m32
fit_ngc4342
fit_ngc4473
fit_double_powerlaw_1d
fit_ngc5831_twist
PRINT, 'Total computation time:', SYSTIME(1) - t, ' Seconds'

END
;----------------------------------------------------------------------------
