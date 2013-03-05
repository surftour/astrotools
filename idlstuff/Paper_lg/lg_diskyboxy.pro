pro diskyboxy, frun, snapnum, xlen


;frun= '/raid4/tcox/vc3vc3b'
;snapnum= 30
;xlen= 10.0



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
hsml= fload_allstars_hsml(1)

hsml= hsml*2.0
print, "WARNING: multiplying the stellar hsml by 2.0"


; --------------------------------------------------------

;  read the list of theta's and phi's

anglefile= '/home/tcox/unitsphereangles.txt'

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

	rotate_phi= phi[i]
	rotate_theta= theta[i]


	; ------------------------------------------------------------------

	set_maxden=10.0
	set_dynrng= 1.0e+5

	istring= strcompress(string(i), /remove_all)
	fitsfilename= frun+'/fits/massdensity'+istring+'.eps'    ; kind of funky becuase
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


            S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/ComputeIsoPhotFitIDL/isophotfit.so', $
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
        
printf, 1, "#   diskyboxy.txt"
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


pro read_diskyboxy, frun, aa, bb, ellip, a4diva

dbfile= frun+'/diskyboxy.txt'

spawn, "wc "+dbfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then db_data= fltarr(6,lines)

openr, 1, dbfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, db_data
close, 1


aa= db_data[2,*]
bb= db_data[3,*]
ellip= db_data[4,*]
a4diva= db_data[5,*]

end



; =================================================================





;------------------------------------------------------------------------
;
;     Density Map of (v/sig)* and vmin/vtot versus a4/a
;     --------------------------------------------------
;
;
;------------------------------------------------------------------------

pro dbmap, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "db_map1, junk"
   print, "  "
   return 
endif

frun="/raid4/tcox/localgroup/v5"

filename='lgdbmap.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3


;--------------------------------------
;--------------------------------------




read_kin_file, frun, theta, phi, ellipi, pai, vmax, sig
thisvsigst= (vmax/sig)/(sqrt(ellipi/(1-ellipi)))
;thisvsigst= [vsigst,transpose(thisvsigst)]
allvsigst_log= alog10(thisvsigst)
print, frun, "   <v/sig*>= ", mean(thisvsigst), "   +/-", sqrt(variance(thisvsigst))


;-------------------------------------

read_diskyboxy, frun, aa, bb, ellip, a4diva
a4s= transpose(a4diva)
print, frun, "    <a4/a>= ", mean(100*a4s), "   +/- ", sqrt(variance(100*a4s))
alla4diva=   100.0 * a4s





;----------------------------------------
; Plot:  (v/sig)* versus a4/a
;----------------------------------------

xaxistitle= '!6100 a!D4!N/a'
xmax= 3.5
xmin= -2.0

yaxistitle= '!6Log (V/!7r!6)!E*!N'
ymax = 0.5
ymin = -1.6

bins= 30

contour_makegeneralpic, alla4diva, allvsigst_log, xmax, xmin, ymax, ymin, $
                                pixels= bins, $
                                NxNImage=NxNImage

; --------------------


x0= 0.20 & y0= 0.15
x1= 0.98 & y1= 0.98

tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

!p.position=[x0, y0, x1, y1]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, color= 0, /noerase, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle=xaxistitle, ytitle=yaxistitle, $
        ycharsize=1.5, xcharsize=1.5, xthick=4.0, ythick=4.0, charthick=3.0


; --------------------

y= [ymin,ymax]
x= [0.0,0.0]
oplot, x, y, psym=-3,color=0,thick=4.0,linestyle=2

y= [0.0,0.0]
x= [xmin,xmax]
oplot, x, y, psym=-3,color=0,thick=4.0,linestyle=2




; ------------------------------------------------------------------



        ; -------------------------
        ; read bender current data
	read_bender_current, name, ellip, a4a, vmaj, vmin, sig, pK, dPA, M_t, SB_e, L_R, L_X

        ; make sure there is a major axis velocity
        idx= where(vmaj lt 800.0)
        if idx(0) ne -1 then begin
                vmaj= vmaj(idx)
                vmin= vmin(idx)
                a4a= a4a(idx)
                ellip= ellip(idx)
		sig= sig(idx)
        endif

        ; make lower limits for vmaj
        idx= where(vmaj lt 0.0)
        if idx(0) ne -1 then vmaj(idx)= -1.0 * vmaj(idx)


	vsigst= vmaj/sig/sqrt(ellip/(1-ellip))
        vsigst= alog10(vsigst)

        print, "a4/a   (max/min) = ", max(a4a), min(a4a)
        idx= where(a4a ge 0.0)
        print, "a4/a percent >= 0.0: ", 100.0*n_elements(idx)/n_elements(a4a)
        idx= where(a4a gt 0.0) 
        print, "a4/a percent >  0.0: ", 100.0*n_elements(idx)/n_elements(a4a)

        f = (2.0*!pi/16.0)*findgen(17)
        usersym,0.7*cos(f),0.7*sin(f),/fill
        oplot, a4a, vsigst, psym=8, color=0, thick=3.0

        ;oplot, [1.05,1.18], [0.86,0.86], psym=-3, color= 0, thick=12.0, linestyle= 1
        ;xyouts, 0.67, 0.73, 'BBF (1992)', /normal, charthick=4.0, size=1.5, color=0



; --------------------

xyouts, 0.8, 0.9, '(d)', /normal, size=1.8, color= 0, charthick=3.0


; --------------------

device, /close


end









