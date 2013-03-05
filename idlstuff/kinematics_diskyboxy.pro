;==============================================================
;
;
;  Look at the isophotal shape from many angles
;   and write to file.
;
;
;
;
;
;
;
;==============================================================
pro diskyboxy, frun, snapnum, xlen

;frun= '/raid4/tcox/vc3vc3b'
;snapnum= 30
;xlen= 10.0

if not keyword_set(frun) then begin
	print, " "
	print, " diskyboxy, frun, snapnum, xlen"
	print, " "
	print, " "
	return
endif



; open the file up
;   and load variables
; ---------------------
ok= fload_snapshot_bh(frun, snapnum)
center=fload_center_alreadycomp(1)

ndisk= fload_npart(2)                               ; disk particles
nbulge= fload_npart(3)                              ; bulge particles
nstars= fload_npart(4)
npart= long(ndisk) + long(nbulge) + long(nstars)
print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars

x= fload_allstars_xyz('x',center=[0,0,0])
y= fload_allstars_xyz('y',center=[0,0,0])
z= fload_allstars_xyz('z',center=[0,0,0])
m= fload_allstars_mass(1)

; need the stellar density (for smoothing the images)
;hsml= fload_allstars_hsml(1)
hsml= m*0.0 + 0.2
;hsml= m*0.0 + 0.1

hsml= hsml*2.0
print, "WARNING: multiplying the stellar hsml by 2.0"


fitsdir= frun+'/snap_'+fload_finfo_exts(1)+'_info/'
print, "fitsdir= ", fitsdir
print, "making directory ", fitsdir
spawn, "mkdir "+fitsdir


; --------------------------------------------------------

;  read the list of theta's and phi's

anglefile= '/n/home/tcox/unitsphereangles.txt'

spawn, "wc "+anglefile,result
lines=long(result)
Nangles=lines(0)-5
if Nangles GT 0 then angle_data= fltarr(2,Nangles)

openr, 1, anglefile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, angle_data
close, 1


theta= angle_data[0,*]
phi= angle_data[1,*]


; --------------------------------------------------------
;
;   Loop through the angles, project and calculate
;  the disky/boxy - ness
;

fit_a= fltarr(Nangles)
fit_b= fltarr(Nangles)
ellip= fltarr(Nangles)
a4diva= fltarr(Nangles)



for i= 0, Nangles-1 do begin

	print, " ------------------------------------------------------------------ "

	rotate_phi= phi[i]
	rotate_theta= theta[i]


	; ------------------------------------------------------------------

	set_maxden=10.0
	set_dynrng= 1.0e+5

	istring= strcompress(string(i), /remove_all)
	;fitsfilename= frun+'/fits/massdensity'+istring+'.eps'    ; kind of funky becuase
	fitsfilename= fitsdir+'/massdensity'+istring+'.eps'    ; kind of funky becuase
							; contour_makepic assumes
							; the ending is eps

        contour_makepic, x, y, z, m, xlen, $
                        hsml=hsml, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        center=center, filename=fitsfilename, /fitstoo, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage, HalfMassSB=HalfMassSB

        Pic= NxNImage
        PIXELS= n_elements(Pic[0,*])



		; do some slight processing
		; --------------------------
		idx=where(Pic eq 1)
                Pic= 256-Pic
                if idx(0) ne -1 then Pic(idx)= 1
                print, "before: min/max= ",min(Pic), max(Pic)
                print, "n > 240= ", n_elements(where(Pic gt 240))
                print, "n < 20= ", n_elements(where(Pic lt 20))
                ; smooth?
                width= 10      ; try other things?
                newMap= smooth(Pic, width, /edge_truncate)
                Pic= newMap
                print, "after: min/max= ",min(Pic), max(Pic)



	; ------------------------------------------------------------------



        ; --------------------------------
        ;
        ;  Fit isophote & calc deviations
        ;
        ; --------------------------------
        calc_isophote_boxydisky= 1
        if calc_isophote_boxydisky eq 1 then begin

            Thresh= 256-HalfMassSB
            SurfaceBrightness= Pic
            Scale= 2.0*xlen

            Xfit=fltarr(360)
            YFit=fltarr(360)

            Xcon=fltarr(360)
            Ycon=fltarr(360)

            AFit=0.0
            BFit=0.0
            X0Fit=0.0
            Y0Fit=0.0
            PhiFit=0.0
            a4Fit=0.0


            spawn, 'echo $HOME', result
            homedir= strcompress(result,/remove_all)
            libfile= homedir+'/Tools/C-Routines_for_IDL/ComputeIsoPhotFitIDL/isophotfit.so'
            S = CALL_EXTERNAL(libfile, $
                      'fit_isophot', $
                      PIXELS, $
                      SurfaceBrightness,$
                      Thresh,$
                      Scale, $
                      Xfit, Yfit,$
                      PhiFit, AFit, BFit, X0Fit, Y0Fit, a4Fit, $
                      Xcon, Ycon)

           ; disky/boxy
           print, "a4/a= ", a4Fit/AFit
	   fit_a[i]= AFit
	   fit_b[i]= BFit
	   if AFit lt BFit then ellip[i]= 1-(AFit/BFit) else ellip[i]= 1-(BFit/AFit)
	   a4diva[i]= a4Fit/AFit


	endif


	; ------------------------------------------------------------------

endfor


; ------------------------------------------------------------------

openw, 1, frun+'/diskyboxy.txt', ERROR=err
        
printf, 1, "#   diskyboxy.txt, "+frun+" snap= "+string(snapnum)
printf, 1, "#   "
printf, 1, "#   "
printf, 1, "#  theta      phi        a       b   "
printf, 1, "#  (deg)    (deg)   (kpc/h) (kpc/h)     ellip     a4/a  "

for i=0,Nangles-1 do begin
	printf, 1, FORMAT= '(2(F8.3," "),4(F8.5,"  "))', $
                theta[i], phi[i], fit_a[i], fit_b[i], ellip[i], a4diva[i]
endfor
close, 1




; ------------------------------------------------------------------



end









; =================================================================


