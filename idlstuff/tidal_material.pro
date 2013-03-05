;============================================================
;
;
;
;
;
;
;============================================================

pro process_hstimg, junk

skipmapping= 0


;--------------------------------------
;
;  load HST i-band image of NGC 2623
;

datadir= "/n/home/tcox/Papers/tidal_tails/"
galdir= "ngc2623/"
HSTibandfile= "ngc2623-f814w.fits"
skipmasktrim= 0
;HSTibandmaskfile= "ngc2623-f814w-mask.fits"
HSTibandmaskfile= "ngc2623-f814w-tjmask.fits"  & skipmasktrim= 1
;HSTibandmaskfile= "ngc2623-f814w-tjgal.fits"   & skipmasktrim= 1
file= datadir+galdir+HSTibandfile
image = mrdfits (file, 0)
;help, image
;print, size(image)


; now load the mask
file= datadir+galdir+HSTibandmaskfile
imagemask = mrdfits (file, 0)


;---------------------
;
;   crop image
;
;
;  note, the hst image have 0.049'' per pixel
;
;  for ngc 2623, cz= 5550 km/sec, or z=0.0185
;  
;  this yields, 1'' = 373.7 pc, or 18.3118 pc/pixel,
;  so we're cropping a section that is 
;      51.8 kpc x 30.1 kpc
;

;
; make it 2830 x 1650
;
croppedimage= image[2000:4829, 2200:3849]
image= croppedimage
writefits, datadir+galdir+"ngc2623-f814w-cutout.fits", image
if skipmasktrim eq 0 then begin
	croppedimagemask= imagemask[2000:4829, 2200:3849]
	imagemask= croppedimagemask
endif



;--------------------
;
;  replace mask stuff with smooth
;  background from upper left
;

bckgrnd_mean= mean(image[0:100,0:100])
bckgrnd_std= sqrt(variance(image[0:100,0:100]))
;bckgrnd_mean= mean(image[2730:2829,1550:1649])
;bckgrnd_std= sqrt(variance(image[2730:2829,1550:1649]))
print, "top 100x100 mean/std= ", bckgrnd_mean, bckgrnd_std
;
;  background (i.e., sky) is usually ~0.1 +/- 0.01
;
;
sz= size(image)
bckgrnd_array= randomn(seed, sz[2], sz[3]) * bckgrnd_std / 2. + bckgrnd_mean

idx=where(imagemask eq 1)
image(idx)= bckgrnd_array(idx)



;--------------------
;
;  subtract sky
;

;image= image - 0.1
;idx=where(image le 0.0)
;if idx(0) ne -1 then image(idx)= 2.0e-2



;--------------------
;
;  smooth??
;

; use median filter
;
;image= median(image, 5)
;image= median(image, 15)

; use average filter
;
image= smooth(image, 10, /edge_truncate)
;image= smooth(image, 15, /edge_truncate)


; sobel?
;image= sobel(image)


; unsharp mask
;image= image - smoothedimage


;---------------------
;
; generate a mask
;

; designed for the smoothedimage
skipmapping= 1
cutoff= 0.12
cutoff= 0.115
;cutoff= 0.015   ; if we subtract out the background
image(where(image ge cutoff))= 2    ; make "galaxy" black
image(where(image lt cutoff))= 1    ; make non-galaxy white



; not too bad with smoothed by 10, and cutoff at 0.12


; add Barry's mask
;image(where(imagemask eq 1))= 150


;---------------------
;
;  map to 256 colors
;

cols= 256


if not skipmapping then begin
	idx= where(image le 0.0)
	if idx(0) ne -1 then image(idx)= min(image(where(image gt 0.0)))
	print, "N(<=0.0)= ", n_elements(idx)

	logpic= alog10(image)
	;image= cols * (logpic - min(logpic)) / (max(logpic) - min(logpic))
	image= cols * (logpic - min(logpic)) / (max(logpic) - min(logpic)) / 1.05   ; scale down a bit
endif



;-------------------------
;
;  map to 3-color image
;

LOADCT, 4
;LOADCT, 13
v1=[0,255]
v2=[0,255]
v3=[0,255]
tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)
TVLCT, c0, c1, c2, /GET


sz= size(image)
img3d= bytarr(3,sz[1],sz[2])
img3d(0,*,*)= c0(image)
img3d(1,*,*)= c1(image)
img3d(2,*,*)= c2(image)



;------------
;
; write jpg
;

jpgfile= "test.jpg"
print, "Saving image " + jpgfile
write_jpeg, jpgfile, img3d, true=1, quality=99

end





;============================================================================




pro process_hstimg2, junk



;--------------------------------------
;
;  load HST i-band image of NGC 2623
;

datadir= "/n/home/tcox/Papers/tidal_tails/"
galdir= "ngc2623/"
HSTibandfile= "ngc2623-f814w-cutout.fits"
file= datadir+galdir+HSTibandfile
image = mrdfits (file, 0)
;help, image
;print, size(image)

idx= where(image le 0.0)
if idx(0) ne -1 then image(idx)= min(image(where(image gt 0.0)))


;--------------------------------------
;
; find radius and mass of each pixel?
;
;

;for now, we'll do it manually
;
; HST image center= 1330, 845

zeropix_x= 1330
zeropix_y= 845

radius= image*0.0
for i= 0, (size(image))[1]-1 do begin
  for j= 0, (size(image))[2]-1 do begin
        dx= 1.*i - 1.*zeropix_x
        dy= 1.*j - 1.*zeropix_y
        radius[i,j]= sqrt(dx*dx + dy*dy) * 18.3118 * 1.0e-3    ; convert to kpc
  endfor
endfor


;
; do this later
;
;idx=where((radius gt 9.95) and (radius lt 10.05))
;image(idx)= 1




;---------------------
;
;  map to 256 colors
;

cols= 256

logpic= alog10(image)
;image= cols * (logpic - min(logpic)) / (max(logpic) - min(logpic))
image= cols * (logpic - min(logpic)) / (max(logpic) - min(logpic)) / 1.05   ; scale down a bit




;---------------------
;
;  draw a white circle
;
idx=where((radius gt 2.98) and (radius lt 3.02))
image(idx)= 1

idx=where((radius gt 7.48) and (radius lt 7.52))
image(idx)= 1


;-------------------------
;
;  map to 3-color image
;

LOADCT, 4
;LOADCT, 13
v1=[0,255]
v2=[0,255]
v3=[0,255]
tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)
TVLCT, c0, c1, c2, /GET


sz= size(image)
img3d= bytarr(3,sz[1],sz[2])
img3d(0,*,*)= c0(image)
img3d(1,*,*)= c1(image)
img3d(2,*,*)= c2(image)






;------------
;
; write jpg
;

jpgfile= "test.jpg"
print, "Saving image " + jpgfile
write_jpeg, jpgfile, img3d, true=1, quality=99

end







;============================================================================




pro process_VATT_img, junk


;--------------------------------------
;
;  load HST i-band image of NGC 2623
;

datadir= "/n/home/tcox/Papers/tidal_tails/"
galdir= "ngc2623/"
;VATTbandfile= "NGC_2623_I_B_tjw2005.fits"
;VATTbandfile= "NGC_2623_I_R_tjw2005.fits"
;VATTbandfile= "NGC_2623_I_U_tjw2005.fits"
;VATTbandfile= "NGC_2623_I_V_tjw2005.fits"
VATTbandfile= "NGC_2623_I_B_tjw2005-cutout.fits"
;VATTbandfile= "NGC_2623_I_R_tjw2005.fits"
;VATTbandfile= "NGC_2623_I_U_tjw2005.fits"
;VATTbandfile= "NGC_2623_I_V_tjw2005.fits"

file= datadir+galdir+VATTbandfile
image = mrdfits (file, 0)
;help, image
;print, size(image)

;HSTibandmaskfile= "ngc2623-f814w-mask.fits"
;HSTibandmaskfile= "ngc2623-f814w-tjmask.fits"  & skipmasktrim= 1
;HSTibandmaskfile= "ngc2623-f814w-tjgal.fits"   & skipmasktrim= 1
;help, image
;print, size(image)


; now load the mask


;---------------------
;
;   crop image
;
;
;
;croppedimage= image[440:809, 560:779]
;image= croppedimage
;writefits, datadir+galdir+"NGC_2623_I_V_tjw2005-cutout.fits", image

end



;============================================================================



pro convert_jpeg_to_fits, junk

datadir= "/n/home/tcox/Papers/tidal_tails/"
galdir= "ngc2623/"

;
; this is the actual mask, i.e., all the stars and
; background or foreground sources
;
;maskjpegfile= datadir+galdir+"maskmaking_Jan20/test_sm10_mask0.12_gimp2.jpg"
;maskfitsfile= datadir+galdir+"ngc2623-f814w-tjmask.fits"


;
; this is the galaxy itself
;
;maskjpegfile= datadir+galdir+"maskmaking_Jan20/test_sm10_mask0.12_gimp0.jpg"
;maskfitsfile= datadir+galdir+"ngc2623-f814w-tjgal.fits"



;
; this is an outline of the galaxy's pixels
;
;maskjpegfile= datadir+galdir+"maskmaking_Jan20/test_sm10_mask0.12_gimp4.jpg"
;maskfitsfile= datadir+galdir+"ngc2623-f814w-tjgaloutline.fits"


read_jpeg, maskjpegfile, image

help, image

sz= size(image)
newimage= bytarr(sz[2],sz[3])
newimage[*,*]= image[0,*,*]
newimage(where(newimage lt 10))= 1
newimage(where(newimage gt 10))= 0

writefits, maskfitsfile, newimage

end





;============================================================================
;============================================================================
;============================================================================
;============================================================================
;============================================================================







pro plot_m_v_r, junk


;--------------------------------------
;
;  load HST i-band image of NGC 2623
;   and the appropriate mask
;

datadir= "/n/home/tcox/Papers/tidal_tails/"
galdir= "ngc2623/"
DDist= 77.0833d+6    ; in pc

; actually, don't need the following two things.
;PHOTFLAM= 6.92557d-20
;PHOTFLAM= 7.07236d-20     ; from header, inverse sensitivity, ergs/s/cm2/Ang/electron
VEGAmag= 25.53561
Msun_band= 4.10          ; I-band luminosity of sun
HSTibandfile= "ngc2623-f814w-cutout.fits"
file= datadir+galdir+HSTibandfile
image = mrdfits (file, 0)

idx= where(image le 0.0)
if idx(0) ne -1 then image(idx)= min(image(where(image gt 0.0)))

HSTibandmaskfile= "ngc2623-f814w-tjgal.fits"
file= datadir+galdir+HSTibandmaskfile
galpiximage = mrdfits (file, 0)


HSTibandmaskfile= "ngc2623-f814w-tjmask.fits"
file= datadir+galdir+HSTibandmaskfile
tjmaskpiximage = mrdfits (file, 0)

;--------------------------------------
;
; find radius and mass of each pixel?
;
;

; hard to get this automatically
;zeropix_idx= where(image eq max(image))
;zeropix_flux= image(zeropix_idx)

;for now, we'll do it manually
;
; HST image center= 1330, 845

zeropix_x= 1330
zeropix_y= 845
print, "zeropix x,y= ", zeropix_x, zeropix_y
zeropix_flux= image[zeropix_x, zeropix_y]
print, zeropix_flux


radius= image*0.0
appmags= image*0.0
absmags= image*0.0
mass= image*0.0
flux= image*0.0
for i= 0, (size(image))[1]-1 do begin
  for j= 0, (size(image))[2]-1 do begin
	dx= 1.*i - 1.*zeropix_x
	dy= 1.*j - 1.*zeropix_y
	radius[i,j]= sqrt(dx*dx + dy*dy) * 18.3118 * 1.0e-3    ; convert to kpc
	;flux[i,j]= image[i,j] * PHOTFLAM                      ; in erg/s/cm^2/Ang
	flux[i,j]= image[i,j]                                  ; in electrons (which the zeropoint is calibrated for)
	appmags[i,j]= -2.5 * alog10(flux[i,j]) + VEGAmag       ; Vega magnitudes
	absmags[i,j]= appmags[i,j] - 5.0 * alog10(DDist) + 5.0
	mass[i,j]= 10.0^((absmags[i,j] - Msun_band) / (-1.0 * 2.5))
  endfor
endfor


;
;  within 3.0 kpc
; ------------------
galidx= where(radius lt 3.0)
print, "NGC 2623 (inner 3 kpc)"
print, "========"
print, "percentage of pixels= ", 1.0*n_elements(galidx)/n_elements(image)
print, "total counts (electrons)= ", total(image(galidx))
print, "max/min counts (electron)= ", max(image(galidx)), min(image(galidx))
print, "total flux= ", total(flux(galidx))


totmag= -2.5 * alog10(total(flux(galidx))) + VEGAmag
print, " "
print, "Apparent Mag = ", -2.5 * alog10(total(flux(galidx))) + VEGAmag
print, "Absolute Mag = ", -2.5 * alog10(total(flux(galidx))) + VEGAmag - 5.0 * alog10(DDist) + 5.0
print, " "

totmass= total(mass(galidx))
print, " "
print, "Mass= ", total(mass(galidx))
print, " "


;
; within 7.5 kpc
; ----------------
galidx= where(radius lt 7.5)
print, "NGC 2623 (inner 7.5 kpc)"
print, "========"
print, "percentage of pixels= ", 1.0*n_elements(galidx)/n_elements(image)
print, "total counts (electrons)= ", total(image(galidx))
print, "max/min counts (electron)= ", max(image(galidx)), min(image(galidx))
print, "total flux= ", total(flux(galidx))


totmag= -2.5 * alog10(total(flux(galidx))) + VEGAmag
print, " "
print, "Apparent Mag = ", -2.5 * alog10(total(flux(galidx))) + VEGAmag
print, "Absolute Mag = ", -2.5 * alog10(total(flux(galidx))) + VEGAmag - 5.0 * alog10(DDist) + 5.0
print, " "

totmass= total(mass(galidx))
print, " "
print, "Mass= ", total(mass(galidx))
print, " "




;
; total magnitude of (entire) system
; --------
galidx= where((galpiximage ge 1) and (tjmaskpiximage lt 1))     ; select the "galaxy's" pixels, but not mask
;galidx= where(galpiximage ge 1)     ; select the "galaxy's" pixels
;galidx= where(galpiximage ge 0)     ; grab them all
totalflux= total(flux(galidx))
print, "================"
print, "NGC 2623 (total)"
print, "================"
print, "percentage of pixels= ", 1.0*n_elements(galidx)/n_elements(flux)
print, "total counts (electrons)= ", totalflux
print, "max/min counts (electron)= ", max(flux(galidx)), min(flux(galidx))
radius= radius(galidx)
flux= flux(galidx)
mass= mass(galidx)
print, "total flux= ", total(flux)


totmag= -2.5 * alog10(total(flux)) + VEGAmag
print, " "
print, "Apparent Mag = ", -2.5 * alog10(total(flux)) + VEGAmag
print, "Absolute Mag = ", -2.5 * alog10(total(flux)) + VEGAmag - 5.0 * alog10(DDist) + 5.0
print, " "

print, " "
print, "Mass= ", total(mass)
print, " "



;
; just outer parts
; -----------------
;galidx= where(radius gt 8.0)
;galidx= where((radius gt 8.0) and (galpiximage ge 1))
;tailflux= total(flux(galidx))
;print, "percentage of pixels= ", 1.0*n_elements(galidx)/n_elements(flux)
;print, "total counts (electrons)= ", tailflux
;print, "tail percentage = ", tailflux/totalflux
;print, "max/min counts (electron)= ", max(flux(galidx)), min(flux(galidx))


;totmag= -2.5 * alog10(total(flux(galidx))) + VEGAmag
;print, " "
;print, "NGC 2623 (outer 8 kpc)"
;print, "========"
;print, "Apparent Mag = ", -2.5 * alog10(total(flux(galidx))) + VEGAmag
;print, "Absolute Mag = ", -2.5 * alog10(total(flux(galidx))) + VEGAmag - 5.0 * alog10(DDist) + 5.0
;print, " "

;radius= radius(galidx)
;flux= flux(galidx)



;stop

flux= alog10(flux)
mass= alog10(mass)

;--------------------------------------
;
; now plot it up
;
;

filename="m_v_r.eps"

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize= 25.0, newysize= 25.0


xmax= 25.0
xmin= 0.0
;
; flux  (ymax/ymin not in log - will be "logged" below)
;ymax= 18.0
ymax= 1.2
ymin= 7.0e-2
;
; mass
ymax= 2.0e+5  & ymax= alog10(ymax)
ymin= 9.0e+3  & ymin= alog10(ymin)
x0=0.12
y0=0.10
x1=0.98
y1=0.98
!p.position=[x0,y0,x1,y1]

;
;
; flux/mass should be in log here 
contour_makegeneralpic, radius, mass, xmax, xmin, ymax, ymin, $
;contour_makegeneralpic, radius, flux, xmax, xmin, ymax, ymin, $
;contour_makegeneralpic, radius, mass, xmax, xmin, alog10(ymax), alog10(ymin), $
;contour_makegeneralpic, radius, flux, xmax, xmin, alog10(ymax), alog10(ymin), $
                                pixels= 156, $
                                NxNImage=NxNImage

tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal



plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
        xcharsize=2.5, ycharsize=2.5, xstyle=1, ystyle=1, /normal, /nodata, $
	charthick=7.0, $
	;charthick=7.0, /ylog, $
	;ytickformat='exp_label', $
        ;xthick=3.0, ythick=3.0, xtitle='R (kpc)', ytitle='Mass (M!D!9n!6!N)'
        xthick=7.0, ythick=7.0, xtitle='!6Radius (kpc)', ytitle='!6Log Mass (M!D!9n!6!N)'
        ;xthick=7.0, ythick=7.0, xtitle='!6Radius (kpc)', ytitle='!6Mass (M!D!9n!6!N)'
        ;xthick=7.0, ythick=7.0, xtitle='!6Radius (kpc)', ytitle='!6Flux'


;oplot, radius, flux, psym=3, thick=5.0, color= 220, symsize= 2.0


;
; flux cutoff
;oplot, [xmin,xmax], [0.115, 0.115], psym=-3, linestyle= 2, color= 0, thick= 3.0
masscut= -2.5 * alog10(0.115) + VEGAmag - 5.0 * alog10(DDist) + 5.0
masscut= 10.0^((masscut - Msun_band) / (-1.0 * 2.5))
;oplot, [xmin,xmax], [masscut, masscut], psym=-3, linestyle= 2, color= 0, thick= 3.0
oplot, [xmin,xmax], [alog10(masscut), alog10(masscut)], psym=-3, linestyle= 2, color= 0, thick= 3.0


;
; inner disk 
oplot, [3.0,3.0], [ymin, ymax], psym=-3, linestyle= 1, color= 0, thick= 3.0

;
; flux cutoff
oplot, [7.5,7.5], [ymin, ymax], psym=-3, linestyle= 1, color= 0, thick= 3.0




; done, close this up
; --------------------
device, /close





end





