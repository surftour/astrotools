;######################################################################
;
; Copyright (C) 1999-2004, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; Updated versions of the software are available from my web page
; http://www-astro.physics.ox.ac.uk/~mxc/idl/
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
;       FIND_GALAXY
;
; AUTHOR:
;       Michele Cappellari, Astrophysics Sub-department, University of Oxford, UK
;
; PURPOSE:
;       Find the largest region of connected pixels (after smoothing)
;       lying above a given intensity level of the image.
;       This is useful to automatically identify the location and orientation of
;       a galaxy in an image, assuming it is the largest positive fluctuation.
;       The convention used by this routine are the same as for the rest
;       of the MGE_FIT_SECTORS package.
;
; EXPLANATION:
;       This procedure uses the weighted first and second moments of the intensity
;       distribution for the computation of the galaxy center and position angle.
;       Further information on FIND_GALAXY is available in
;       Cappellari M., 2002, MNRAS, 333, 400
;
; CALLING SEQUENCE:
;       FIND_GALAXY, Img, Majoraxis, Eps, Ang, Xpeak, Ypeak, Xmed, Ymed, $
;           FRACTION=fraction, INDEX=ind, /PLOT, NBLOB=nblob, /QUIET, X=x, Y=y
;
; INPUTS:
;       Img = The galaxy images as a 2D array.
;
; OUTPUTS:
;       Eps = The galaxy "average" ellipticity Eps = 1-b/a = 1-q'.
;           Photometry will be measured along elliptical annuli with
;           constant axial ellipticity Eps. The four quantities
;           (Eps,Ang,Xc,Yc) can be measured with the routine FIND_GALAXY.
;       Ang = Position angle measured from the image Y axis,
;           counter-clockwise to the galaxy major axis.
;           All angles are measured with respect to this direction
;       Xpeak = X coordinate of the galaxy center in pixels. To be precise this
;           coordinate represents the coordinate of the brightest pixels within
;           a 40x40 pixels box centered on the galaxy weighted average center.
;       Ypeak = Y coordinate of the galaxy center in pixels.
;       Xmed = X coordinate of luminosity weighted galaxy center.
;       Ymed = Y coordinate of luminosity weighted galaxy center.
;
; OPTIONAL INPUT KEYWORDS:
;       FRACTION - This corresponds (approximately) to the fraction
;           [0 < FRACTION < 1] of the image pixels that one wants to
;           consider to estimate galaxy parameters (default 0.1 = 10%)
;       PLOT - display an image in the current graphic window showing
;           the pixels used in the computation of the moments.
;       NBLOB - If NBLOB=1 (default) the procedure selects the largest feature
;           in the image, if NBLOB=2 the second largest is selected, and so
;           on for increasing values of NBLOB. This is useful when the
;           galaxy is not the largest feature in the image.
;       QUIET - do not print numerical values on the screen.
;
; OPTIONAL OUTPUT KEYWORDS:
;       INDEX - img[index] contain the pixels of the galaxy from which the
;           inertia ellipsoid and the ellipticity and PA are derived.
;       X - The x-coordinate of the pixels with index given above.
;       Y - The y-coordinate of the pixels with index given above.
;
; EXAMPLE:
;       This command locates the position and orientation of a galaxy
;       in the image IMG and produces a plot showing the obtained results
;       and the region of the image used in the computation.
;
;           find_galaxy, img, /PLOT
;
;       This command only uses the 2% of the image pixels to estimate the
;       intensity weighted moments and show the results.
;
;           find_galaxy, img, FRACTION=0.02, /PLOT
;
; PROCEDURES USED:
;       The following procedures are contained in the main FIND_GALAXY program.
;           DISPLAY        -- Show images even larger than the screen by resampling
;           SECOND_MOMENTS -- Compute second moments on the region selected
;
;       Astronomy User's Library routines needed (http://idlastro.gsfc.nasa.gov/):
;           TVELLIPSE  -- Draw an ellipse on the current display
;
; MODIFICATION HISTORY:
;       V1.0: Written by Michele Cappellari, Padova, 30 Giugno 1999
;       V1.01: Michele Cappellari, ESO Garching, 27 september 1999
;       V1.1: Made a more robust automatic level selection, MC, Leiden, August 2001
;       V1.11: Added compilation options, MC, Leiden 20 May 2002
;       V1.12: Load proper color table for plot. MC, Leiden, 9 May 2003
;       V1.2: Do not use a widget to display the image. Just resample the
;           image if needed. Added /QUIET keyword and (xmed,ymed) output.
;           After suggestion by R. McDermid. MC, Leiden, 29 July 2003
;       V1.21: Make the procedure work also with very small images,
;           where it does not make sense to extract a subimage.
;           MC, Leiden, 1 August 2003
;       V1.22: Added NBLOB keyword. Useful when the galaxy is not the
;           largest feature in the image. MC, Leiden, 9 October 2004
;       V1.23: Gives an error message if IMG is not a two-dimensional array.
;           MC, Leiden, 11 May 2006
;       V1.24: Included optional output keywords INDEX, X and Y.
;           MC, Munich, 14 December 2006
;-
;----------------------------------------------------------------------------
PRO display, img, ratio
COMPILE_OPT IDL2, HIDDEN

swx = !D.X_VSIZE    ; Window size in X in device units
swy = !D.Y_VSIZE    ; Size in Y
s = SIZE(img)
six = FLOAT(s[1])   ; Image size in X in pixels
siy = FLOAT(s[2])   ; Size in Y

ratio = siy/swy > six/swx
TVSCL, POLY_2D(img, [[0,0],[ratio,0]],[[0,ratio],[0,0]], 0, six/ratio, siy/ratio)

END
;--------------------------------------------------------------------
PRO second_moments, img, ind, majoraxis, eps, theta, xpeak, ypeak, xmed, ymed, x, y
COMPILE_OPT IDL2, HIDDEN
;
; Restrict the computation of the first and second moments to
; the region containing the galaxy, defined by vector IND.

img1 = img[ind]
s = SIZE(img)
x = FLOAT(ind MOD s[1])
y = FLOAT(ind/s[1])

; Compute coefficients of the moment of inertia tensor.
;
i = TOTAL(img1)
xmed = TOTAL(img1*x)/i
ymed = TOTAL(img1*y)/i
x2 = TOTAL(img1*x^2)/i - xmed^2
y2 = TOTAL(img1*y^2)/i - ymed^2
xy = TOTAL(img1*x*y)/i - xmed*ymed

; Diagonalize the moment of inertia tensor.
; theta is the angle, measured counter-clockwise,
; from the image Y axis to the galaxy major axis.
;
theta = ATAN(2.0*xy,x2-y2)/2.0 * !RADEG + 90.0
a2 = (x2+y2)/2.0 + SQRT(((x2-y2)/2.0)^2 + xy^2)
b2 = (x2+y2)/2.0 - SQRT(((x2-y2)/2.0)^2 + xy^2)
eps = 1.0 - SQRT(b2/a2)
majoraxis = SQRT(a2)

; If the image has many pixels then compute the coordinates of the
; highest pixel value inside a 40x40 pixels region centered on the
; first intensity moments (Xmed,Ymed), otherwise just return the
; coordinates of the highest pixel value in the whole image.
;
n = 20
xmed1 = ROUND(xmed)
ymed1 = ROUND(ymed)
IF xmed1-n GE 0 AND xmed1+n LT s[1] AND $ ; Check if subimage fits...
   ymed1-n GE 0 AND ymed1+n LT s[2] THEN BEGIN
    tmp = MAX(img[xmed1-n:xmed1+n,ymed1-n:ymed1+n],j)
    xpeak = (j MOD (2*n+1)) + xmed1-n
    ypeak = j/(2*n+1) + ymed1-n
ENDIF ELSE BEGIN    ; ...otherwise use full image
    tmp = MAX(img,j)
    xpeak = j MOD s[1]
    ypeak = j/s[1]
ENDELSE

END
;-------------------------------------------------------------------------
PRO find_galaxy, img, majoraxis, eps, theta, xpeak, ypeak, xmed, ymed, $
     FRACTION=fraction, PLOT=PLOT, QUIET=quiet, NBLOB=nblob, INDEX=ind, X=x, Y=y
COMPILE_OPT IDL2
ON_ERROR, 2
;
; With nblob=1 find the largest connected region in the image,
; with nblob=2 find the second in size and so on...

if n_elements(nblob) eq 0 then nblob = 1
s = SIZE(img)
if s[0] ne 2 then message, 'IMG must be a two-dimensional array'
a = MEDIAN(img,5)
j = SORT(a)
IF N_ELEMENTS(fraction) EQ 0 THEN fraction = 0.1
j = WHERE(a GT a[j[s[4]*(1.0-fraction)]])

a = BYTARR(s[1],s[2])
a[j] = 1
a = LABEL_REGION(a,/EIGHT)             ; Get blob indices.
h = HISTOGRAM(a, REVERSE_INDICES=r)    ; Get population and members of each blob.
gal = (SORT(h))[N_ELEMENTS(h)-1-nblob] ; Find the nblob-th largest region (after the sky)
ind = r[r[gal]:r[gal+1]-1]             ; Find subscripts of members of region "gal"

second_moments, img, ind, majoraxis, eps, theta, xpeak, ypeak, xmed, ymed, x, y

IF NOT KEYWORD_SET(quiet) THEN BEGIN
    PRINT, ' Pixels used:', n_elements(ind)
    PRINT, ' Peak (x,y):', xpeak, ypeak
    PRINT, ' Mean (x,y):', xmed, ymed, FORMAT='(a,2f10.2)'
    PRINT, ' Theta (deg):', theta, FORMAT='(a,f10.2)'
    PRINT, ' Eps:', eps, FORMAT='(a,f10.3)'
    PRINT, ' Sigma along major axis (pix):', majoraxis, FORMAT='(a,f10.2)'
ENDIF

IF KEYWORD_SET(plot) THEN BEGIN
    LOADCT, 2, /SILENT
    a = BYTARR(s[1],s[2])
    a[ind] = 1
    display, a, ratio
    mjr = 3.5*majoraxis/ratio
    xc = xmed/ratio
    yc = ymed/ratio
    coef = tan(theta/!radeg)
    tvellipse, mjr, mjr*(1-eps), xc, yc, theta-90, COLOR=73, THICK=3
    xx = FINDGEN(s[1])/ratio
    PLOTS, xx, yc - (xx-xc)/coef, COLOR=28, /DEVICE, THICK=3
    PLOTS, xx, yc + (xx-xc)*coef, COLOR=28, /DEVICE, THICK=3
ENDIF

END
;-------------------------------------------------------------------------
