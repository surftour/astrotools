;######################################################################
;
; Copyright (C) 1999-2002, Michele Cappellari
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
;     SECTORS_PHOTOMETRY_TWIST
;
; AUTHOR:
;       Michele Cappellari, Astrophysics Sub-department, University of Oxford, UK
;
; PURPOSE:
;       Perform photometry of a galaxy image along sectors equally spaced in
;       angle. This routine assumes point-symmetry, so measurements symmetric
;       with respect to the center are averaged together. This routine is
;       useful to generate the input photometry required by the MGE fitting
;       routine MGE_FIT_SECTORS_TWIST.
;
; EXPLANATION:
;       Further information on SECTORS_PHOTOMETRY_TWIST is available in
;       Cappellari M., 2002, MNRAS, 333, 400
;
; CALLING SEQUENCE:
;       SECTORS_PHOTOMETRY_TWIST, Img, Ang, Xc, Yc, Radius, Angle, Counts,
;           N_SECTORS=n_sectors, SECTOR_WIDTH=sector_width,
;           BADPIXELS=badpixels, MINLEVEL=minlevel
;
; INPUTS:
;       Img = The galaxy images as a 2D array.
;       Ang = Position angle measured from the image Y axis,
;           counterclockwise to the galaxy major axis.
;           All angles are measured with respect to this direction.
;       Xc = X coordinate of the galaxy center in pixels.
;       Yc = Y coordinate of the galaxy center in pixels.
;           The three quantities (Ang,Xc,Yc) can be
;           measured with the routine FIND_GALAXY.
;
; OUTPUTS:
;       Radius = Vector containing the radius of the surface brightness
;               measurements, taken from the galaxy center. This is given
;               in units of PIXELS (!).
;       Angle = Vector containing the polar angle of the surface brightness
;               measurements, taken from the galaxy major axis.
;       Counts = Vector containing the actual surface brightness measurements
;               in COUNTS (!) at the polar coordinates specified by the vectors
;               Radius and Angle. These three vectors need to have the same
;               number of elements.
;
; OPTIONAL INPUT KEYWORDS:
;       N_SECTORS - Number of the sectors, equally spaced in angle, from
;               the galaxy major axis to the minor axis (one quadrant).
;               (default: 36 to cover the whole image with 5 degrees width).
;       SECTOR_WIDTH - Scalar giving the width of the sectors
;               in degrees (default: 5 degrees)
;       BADPIXELS - Vector containing the indices of the pixels in the
;               array Img that have to be excluded from the photometry
;               calculation. Pixels indices can be selected using the
;               WHERE or POLYFILLV commands or interactively with
;               the DEFROI routine.
;       MINLEVEL - Scalar giving the minimum COUNTS level to include
;               in the photometry. The measurement along one profile
;               will stop when the counts go below this level.
;
; EXAMPLE:
;       The sequence of commands below was used to generate the complete
;       MGE model of the Figures 11-12 in Cappellari (2002).
;       1) The FITS file is read sky, the bad regions are masked and
;           the image is sky subtracted;
;       2) the photometry along sectors is performed with SECTORS_PHOTOMETRY_TWIST;
;       3) the resulting measurements are fitted with MGE_FIT_SECTORS_TWIST;
;       4) the contour are printed on a PS file with MGE_PRINT_CONTOURS_TWIST.
;       The geometric parameters of the galaxy (ang,xc,yc) were
;       previously determined using FIND_GALAXY.
;
;           fits_read, 'ngc5831_f702w_mosaic.fits', img, h
;
;           badpixels = where(img le 0) ; mask the unexposed part of the mosaic
;           ;
;           ; Mask a galaxy and two bright stars
;           ;
;           s = SIZE(img)
;           dist_circle, r, s[1:2], 1408.09, 357.749
;           badpixels = [badpixels,WHERE(r LT 200)]
;           dist_circle, r, s[1:2], 539.243, 161.264
;           badpixels = [badpixels,WHERE(r LT 50)]
;           dist_circle, r, s[1:2], 393.028, 713.681
;           badpixels = [badpixels,WHERE(r LT 50)]
;
;           skylevel = 13.0 ; In counts
;           img = img - skyLevel
;           scale = 0.100   ; WFPC2/WF CCD pixel scale arcsec/pixels
;
;           sigmaPSF = 0.4  ; WF pixels: here I use one single Gaussian PSF
;           eps = 0.28      ; These values were measured with FIND_GALAXY
;           ang = 141.0     ; major axis PA in the inner regions
;           xc = 969
;           yc = 974
;
;           sectors_photometry_twist, img, ang, xc, yc, radius, counts, angle, $
;               BADPIXELS=badpixelsMosaic, MINLEVEL=2.0
;
;           MGE_fit_sectors_twist, radius, angle, counts, eps, $
;               NGAUSS=11, SIGMAPSF=sigmaPSF, SOL=sol, /PRINT, $
;               SCALE=scale, QBOUNDS=[0.2,1]
;
;           MGE_print_contours_twist, img, ang, xc, yc, sol, $
;               FILE='NGC5831_Mosaic.ps', BINNING=5, MAGRANGE=9.5, $
;               SIGMAPSF=sigmaPSF, SCALE=scale
;
; PROCEDURES USED:
;       The following procedures are contained in the main SECTORS_PHOTOMETRY program.
;           ANGLES_IMAGE_TWIST               -- construct an image of position angles
;           ELLIPSES_IMAGE_TWIST             -- construct an image of elliptical radii
;           SECTOR_PHOTOMETRIC_PROFILE_TWIST -- perform photometry along ONE position angle
;
;       Astronomy User's Library (http://idlastro.gsfc.nasa.gov/) routines needed:
;           BIWEIGHT_MEAN  -- Computes biweight mean of vectors.
;
; MODIFICATION HISTORY:
;       V1.0 First implementation for the NGC2681 photometric modeling.
;           Michele Cappellari, ESO Garching, 27 september 1999
;       V2.0 Major revisions, to use it with MGE_FIT_SECTORS.
;           Leiden, January 2000, MC
;       V2.1 Further updates, Padova, August 2000, MC
;       V2.11 Added compilation options, MC, Leiden 20 May 2002
;       V2.12 Changed input parameters order according to the
;           documentation to be consistent with SECTORS_PHOTOMETRY
;           MC, Leiden, 6 June 2002
;       V2.13 Use biweight mean insted of sigma-clipped mean.
;           MC, Leiden, 30 April 2004
;       V2.14 Removed the un-necessary parameter EPS from this routine.
;           MC, Leiden, 1 September 2004
;       V2.15 Reduced amount of verbose output. MC, Leiden, 24 October 2004
;       V2.16 Replaced LOGRANGE keyword in example with the new MAGRANGE.
;           MC, Leiden, 1 May 2005
;-
;----------------------------------------------------------------------------
FUNCTION angles_image_twist, pos_ang, xc, yc, n
;
; By Michele Cappellari, Padova, August 2000
;   Modified for triaxial case, August 2001
;
COMPILE_OPT IDL2, HIDDEN

rot = pos_ang/!RADEG
x = FINDGEN(n[0]) - xc
y = FINDGEN(n[1]) - yc

im = FLTARR(n[0], n[1], /NOZERO)
FOR i=0,n[1]-1 DO im[0,i] = ATAN(x,y[i])
;
; Put all angles in the range [-90,+90] deg.
; It has point-symmetry as a triaxial Galaxy
;
ang = ATAN(TAN(im + rot))*!RADEG

RETURN, ang
END
;----------------------------------------------------------------------------
FUNCTION ellipses_image_twist, ellip, pos_ang, xc, yc, n
;
; Returns a 2D image with the elliptical radius corresponding to every pixel.
; The image has size N[0]xN[1], center (XC,YC), ellipticity
; ELLIP (= 1-b/a) and position angle POS_ANG, measured from the positive
; Y axis to the ellipse major axis (positive counter-clockwise).
; By Michele Cappellari, Leiden, January 2000
;
COMPILE_OPT IDL2, HIDDEN

ang = pos_ang/!RADEG

x = FINDGEN(n[0]) - xc
y = FINDGEN(n[1]) - yc

xcosang = COS(ang) * x
xsinang = SIN(ang) * (1.0-ellip) * x
ycosang = COS(ang) * (1.0-ellip) * y
ysinang = SIN(ang) * y

im = FLTARR(n[0], n[1], /NOZERO)

FOR i=0,n[1]-1 DO BEGIN
   xtemp =  xcosang + ysinang[i]
   ytemp = -xsinang + ycosang[i]
   im[0,i] = xtemp^2 + ytemp^2
ENDFOR

RETURN, SQRT(im)
END
;----------------------------------------------------------------------------
PRO sector_photometric_profile_twist, data, ang, xc, yc, $
    radius, cnt, img2, img3, sigma_cnt, BADPIXELS = badpixels, GAIN = gain, $
    SECTOR = sector, PRINT = print, MINLEVEL = minlevel

COMPILE_OPT IDL2, HIDDEN

COMMON profilo, tmp, m1

IF NOT KEYWORD_SET(gain) THEN gain = 1.0
IF NOT KEYWORD_SET(sector) THEN sector = [-90.0, 90.0] ; All image by default
IF NOT KEYWORD_SET(minlevel) THEN minlevel = 0.0

s = SIZE(data)
xc = ROUND(xc)
yc = ROUND(yc)

img3[xc,yc] = sector[0] ; all sectors include the central pixel
;
; This is a simple trick to make sure that the innermost
; 21x21 pixels are properly sampled in the profile
;
j = INDGEN(21)-10
m = TAN((ang-sector[0])/!RADEG)
IF ABS(m) LT 1 THEN $
    img3[ROUND(xc-j*m),yc+j] = sector[0] $
ELSE $
    img3[xc+j,ROUND(yc-j/m)] = sector[0]

IF N_ELEMENTS(badpixels) GT 1 THEN data[badpixels] = -123 ; is an arbitrary flag value
IF ABS(sector[0]) EQ 90 THEN BEGIN ; avoids counting twice -90 and +90
    good = WHERE(ABS(img3) GE ABS(sector[0]) - sector[1]/2.0 AND data NE -123, m)
ENDIF ELSE BEGIN
    good = WHERE(img3 GE sector[0] - sector[1]/2.0 AND $
                 img3 LE sector[0] + sector[1]/2.0 AND $
                 data NE -123, m)
ENDELSE
IF m EQ 0 THEN MESSAGE, 'Input error: No good pixels in the image!'

img = ROUND(24.2*ALOG10(img2[good])) ; 24 isophotes per decade: factor 1.1 spacing
h = HISTOGRAM(img, REVERSE=ind) ; Get population and members of each "ellipse"
w = WHERE(h NE 0, m)

cnt = FLTARR(m)
sigma_cnt = cnt
radius = cnt

IF KEYWORD_SET(print) THEN $                ; debug
    IF sector[0] EQ -90 THEN BEGIN          ; debug
        s = SIZE(data)                      ; debug
        WINDOW, XSIZE=s[1], YSIZE=s[2]      ; debug
        tmp = data                          ; debug
        m1 = MAX(MEDIAN(tmp[good],3))/10    ; debug
    ENDIF                                   ; debug

FOR j=0,N_ELEMENTS(w)-1 DO BEGIN
    sub = ind[ind[w[j]]:ind[w[j]+1]-1]  ; Find subscripts of members of region "r[j]"
    sub = good[sub]                     ; Use the indices of the original image
    n = N_ELEMENTS(sub)
    IF  KEYWORD_SET(print) THEN BEGIN       ; debug
        tmp[sub] = m1*(w[j] MOD 2)          ; debug
        TVSCL, ALOG(tmp>1e-3)               ; debug
    ENDIF                                   ; debug
    IF n GE 10 THEN BEGIN ; Evaluate a biweight mean
        cnt[j] = BIWEIGHT_MEAN(data[sub],sigma)
        sigma_cnt[j] = sigma/sqrt(n)
    ENDIF ELSE BEGIN
        cnt[j] = TOTAL(data[sub])/n ; Usual mean
        sigma_cnt[j] = SQRT(TOTAL(data[sub])*gain)/n ; Poissonian error
    ENDELSE
    IF cnt[j] LT minlevel THEN BEGIN
        cnt = cnt[0:j-1]
        radius = radius[0:j-1]
        GOTO, fine ; In IDL 5.4 this was a BREAK
    ENDIF
    radius[j] = TOTAL(img2[sub]*data[sub])/total(data[sub])  ; Weighted average radius in pixels
ENDFOR
fine:

j = SORT(radius)
cnt = cnt[j]
radius = radius[j]

END
;----------------------------------------------------------------------------
PRO sectors_photometry_twist, img, ang, xc, yc, radius, angle, cnt, $
    BADPIXELS=badpixels, SECTOR_WIDTH=sector_width, $
    N_SECTORS=n_sectors, MINLEVEL=minlevel
;
; This routine performs photometry along sectors linearly spaced
; in angle between the major and minor axis of a galaxy.
; In output it returns the three vectors RADIUS, ANGLE, CNT,
; containing the photometric measurements in polar coordinates.
;
COMPILE_OPT IDL2
ON_ERROR, 2

IF N_ELEMENTS(n_sectors) EQ 0 THEN n_sectors = 36 ELSE n_sectors = n_sectors<72
IF N_ELEMENTS(sector_width) EQ 0 THEN sector_width = 5 ; degrees (5==>[-2.5,+2.5])
IF N_ELEMENTS(minlevel) EQ 0 THEN minlevel = 0.0

s = SIZE(img)
img2 = ellipses_image_twist(0.0, 0.0, xc, yc, s[1:2]) ; eps=0 here: this is the true radius
img2[xc,yc] = 0.38  ; This is the average radius within the central pixel
img3 = angles_image_twist(ang, xc, yc, s[1:2])

PRINT, 'Measuring ', STRTRIM(n_sectors,2), ' profiles of width ', $
    STRING(sector_width,FORMAT='(G0.2)'), ', equally spaced by ', $
    STRING(180./n_sectors,FORMAT='(G0.2)'), ' from -90 to +90 degrees...'

; Samples profiles from -90 to 90 degrees but 90 coincides
; with -90 and is not computed twice.
;
FOR j=0,n_sectors-1 DO BEGIN
    sector_ang = 180.0/n_sectors*j - 90.0
    sector_photometric_profile_twist, img, ang, xc, yc, radius1, cnt1, img2, img3, $
        SECTOR=[sector_ang, sector_width], BADPIXELS=badpixels, MINLEVEL=minlevel
    IF j EQ 0 THEN BEGIN
        radius = radius1
        cnt = cnt1
        angle = REPLICATE(sector_ang,N_ELEMENTS(radius1))
    ENDIF ELSE BEGIN
        radius = [radius, radius1]
        cnt = [cnt, cnt1]
        angle = [angle, REPLICATE(sector_ang,N_ELEMENTS(radius1))]
    ENDELSE
ENDFOR

END
;----------------------------------------------------------------------------
