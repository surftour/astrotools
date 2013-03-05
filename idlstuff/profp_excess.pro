;------------------------------------------------------------------------
;
;    Random Procedures related to Profiles Paper
;
;
;
;
;------------------------------------------------------------------------


;   average sigmas
; --------------------------
pro fload_fruns_sigma_e, frun, snapnum, xvar2

   ns= n_elements(frun)
   allsig= fltarr(ns)

   for i=0,ns-1 do begin

        read_sigma_file, frun[i], time, Asigxy, Asigxz, Asigyz, Asigavg, Asigerr, $
                                Bsigxy, Bsigxz, Bsigyz, Bsigavg, Bsigerr, $
                                gsigavg, gsigerr
        slstidx= n_elements(time)-1
        allsig[i]= Asigavg[slstidx]
   endfor

   xvar= allsig

end





pro fload_fruns_sersicn, frun, snapnum, sersicn, sersic_1sig, $
			excesslight, excesslight_1sig, $
			re, re_1sig

; make sure these are same as
; calling routine
;
bins = 50

xmax = 2.4
;xmax = 4.2
;xmin = 0.3
xmin = 0.5
;ymax = 6.0
;ymin = -1.0


   ns= n_elements(frun)
   sersicn= fltarr(ns)
   sersic_1sig= fltarr(ns)
   excesslight= fltarr(ns)
   excesslight_1sig= fltarr(ns)
   re= fltarr(ns)
   re_1sig= fltarr(ns)

   for i=0,ns-1 do begin
        ;
        ; total (all star) profile
        ; --------------------------

        ok= fload_snapshot_bh(frun[i],snapnum[i])

        x=fload_allstars_xyz('x')
        y=fload_allstars_xyz('y')
        z=fload_allstars_xyz('z')
        m=fload_allstars_mass(1)


        r_devac = fltarr(bins)
        mass_sd_avg = fltarr(bins)
        mass_sd_1sig = fltarr(bins)

        process_averages_formanyprojections, x, y, z, m, $
		xmin, xmax, bins, /x_is_devac, $
		mass_sd_avg, mass_sd_1sig, r_devac, $
		sersic_n_avg=sersic_n_avg, /fitsersic, $
		sersic_n_1sig=sersic_n_1sig, /dontoverplot, $
		exslt_avg=exslt_avg, exslt_1sig= exslt_1sig,$
		r_e_avg=r_e_avg, r_e_1sig=r_e_1sig

	sersicn[i]= sersic_n_avg
	sersic_1sig[i]= sersic_n_1sig
	excesslight[i]= exslt_avg
	excesslight_1sig[i]= exslt_1sig
	re[i]= r_e_avg
	re_1sig[i]= r_e_1sig

   endfor

end





;============================================================================





pro n_excess, junk



if not keyword_set(junk) then begin
        print, " "
        print, " n_excess, junk"
        print, " "
        print, " "
        return
endif



;------------------------------------------------------------------------------------


; gas fraction
; -------------
;do_gf= 0
do_gf= 1
if do_gf eq 1 then begin
	xvar1= [0,0.05,0.20,0.40,0.60,0.80,1.00]
	xvar2= xvar1
	xvar3= xvar1
	xmax= 1.1 & xmin= -0.1
	xaxistitle='Progenitor Gas Fraction'


	; e
	;
	frun= ['/raid4/tcox/cvc3vc3e', $
		'/raid4/tcox/gfs/vc3vc3z_e', $
		'/raid4/tcox/gfs/vc3vc3y_e', $
		'/raid4/tcox/gfs/vc3vc3x2_e', $
		'/raid4/tcox/gfs/vc3vc3w_e', $
		'/raid4/tcox/gfs/vc3vc3v_e', $
		'/raid4/tcox/gfs/vc3vc3u_e']
	snapnum= [30,30,30,30,30,30,30]
	fload_fruns_sersicn, frun, snapnum, sersicn1, sersicn_1sig1, $
		excesslight1, exslgt_1sig1, $
		res1, res_1sig1

	;sersicn1= [2.71,3.24, 3.46, 4.79, 4.60, 4.49, 4.25]
	;excesslight1= [2.7,-0.5,-0.1,-4.4,4.5,0.5,10.5]
	;res1= [4.4,3.8, 3.4, 2.1, 1.9, 1.3, 1.1]
	; not minimum radius for fits
	;sersicn1= [3.05,3.52, 3.88, 7.09, 6.14, 7.80, 6.17]
	;excesslight1= [2.7,-0.5,-0.1,-4.3,4.5,0.5,10.6]
	;res1= [4.4,3.7, 3.2, 1.7, 1.5, 0.7, 0.6]
	; min r= 1.0 kpc/h
	;sersicn1= [3.83,3.80, 4.02, 4.63, 4.55, 6.02, 3.67]
        ;excesslight1= [-16.2,-8.0,-3.4,26.7,28.2,25.9,50.5]
        ;res1= [3.8,3.6, 3.1, 2.6, 2.1, 1.2, 1.4]


	; h
	;
	frun= ['/raid4/tcox/cvc3vc3h', $
		'/raid4/tcox/gfs/vc3vc3z_h', $
		'/raid4/tcox/gfs/vc3vc3y_h', $
		'/raid4/tcox/gfs/vc3vc3x2_h', $
		'/raid4/tcox/gfs/vc3vc3w_h', $
		'/raid4/tcox/gfs/vc3vc3v_h', $
		'/raid4/tcox/gfs/vc3vc3u_h']
	snapnum= [30,30,30,30,30,30,30]
	fload_fruns_sersicn, frun, snapnum, sersicn2, sersicn_1sig2, $
		excesslight2, exslgt_1sig2, $
		res2, res_1sig2

	; no minimum radius for fit
	;sersicn2= [ 3.12, 4.07, 5.59, 6.68, 7.6, 5.6, 4.49]
	;excesslight2= [ 13.8, 4.8, 2.1, 0.1, 0.02, 4.1, 3.35]
	;res2= [ 4.46, 4.06, 3.55, 2.28, 1.29, 1.25, 0.95]
	; 1 kpc/h is the minimum radius for fit
	;sersicn2= [ 4.84, 5.15, 4.99, 6.15, 5.58, 4.64, 2.76]
	;excesslight2= [ -26.4, -15.0, -4.1, 7.4, 24.0, 21.7, 40.5]
	;res2= [ 3.86, 3.74, 3.43, 2.48, 1.94, 1.63, 1.79]


	; k
	;
	frun= ['/raid4/tcox/cvc3vc3k', $
		'/raid4/tcox/gfs/vc3vc3z_k', $
		'/raid4/tcox/gfs/vc3vc3y_k', $
		'/raid4/tcox/gfs/vc3vc3x2_k', $
		'/raid4/tcox/gfs/vc3vc3w_k', $
		'/raid4/tcox/gfs/vc3vc3v_k', $
		'/raid4/tcox/gfs/vc3vc3u_k']
	snapnum= [30,30,30,30,30,30,30]
	fload_fruns_sersicn, frun, snapnum, sersicn3, sersicn_1sig3, $
		excesslight3, exslgt_1sig3, $
		res3, res_1sig3

	; no minimum radius for fit
	;sersicn3= [ 2.44, 3.04, 3.90, 5.73, 7.31, 3.77, 4.75]
	;excesslight3= [ 11.2, 1.0, -0.1, -1.5, 0.5, 27.7, 10.1]
	;res3= [ 4.22, 3.76, 2.9, 2.4, 1.5, 2.1, 1.0]
	; 1 kpc/h is the minimum radius for fit
	;sersicn3= [ 2.57, 2.59, 3.28, 3.19, 4.28, 2.83, 2.78]
	;excesslight3= [ 6.4, 15.8, 15.3, 38.3, 32.9, 45.0, 47.4]
	;res3= [ 4.13, 4.07, 3.28, 3.48, 2.40, 2.71, 2.0]

endif



; do masses
; -------------
;do_mass= 1
do_mass= 0
if do_mass eq 1 then begin
	xvar= [62.0,105.0,146.0,178.0,241.0,367.0,552.0]
	xvar1= xvar
	xvar2= xvar
	xvar3= xvar
        xmax= 600.0 & xmin= 30.0
        xaxistitle='!7r!6 (km s!E-1!N)'


        ; e
        ;
	frun= ['/raid4/tcox/ds/d0e', $
		'/raid4/tcox/ds/d1e', $
		'/raid4/tcox/ds/d2e', $
		'/raid4/tcox/ds/d3e7', $
		'/raid4/tcox/ds/d4e', $
		'/raid4/tcox/ds/d5e', $
		'/raid4/tcox/ds/d6e']
	snapnum= [50,30,30,30,30,30,30]
	;fload_fruns_sigma_e, frun, snapnum, xvar1
	fload_fruns_sersicn, frun, snapnum, sersicn1, sersicn_1sig1, $
		excesslight1, exslgt_1sig1, $
		res1, res_1sig1

	; no minimum radius for fit
        ; min r= 1.0 kpc/h
        ;sersicn1= [3.12, 3.48, 2.43, 5.93, 6.7, 8.78, 10.0]
        ;excesslight1= [20.6, 35.1, 38.6, -6.1, 7.1, -4.6, -56.9]
        ;res1= [ 1.5, 1.9, 2.7, 2.4, 4.1, 6.6, 12.4]


        ; h
        ;
	frun= ['/raid4/tcox/ds/d0h', $
		'/raid4/tcox/ds/d1h', $
		'/raid4/tcox/ds/d2h', $
		'/raid4/tcox/ds/d3h7', $
		'/raid4/tcox/ds/d4h', $
		'/raid4/tcox/ds/d5h', $
		'/raid4/tcox/ds/d6h']
	snapnum= [50,30,30,30,30,30,30]
	;fload_fruns_sigma_e, frun, snapnum, xvar2
	fload_fruns_sersicn, frun, snapnum, sersicn2, sersicn_1sig2, $
		excesslight2, exslgt_1sig2, $
		res2, res_1sig2

        ; 1 kpc/h is the minimum radius for fit
        ;sersicn2= [ 4.10, 4.42, 5.52, 7.46, 5.47, 10.0, 10.0]
        ;excesslight2= [ 2.98, 29.9, 29.8, -4.8, 33.3, -58.6, -48.4]
        ;res2= [ 1.2, 1.9, 2.4, 2.96, 5.98, 5.1, 13.1]


        ; k
        ;
	frun= ['/raid4/tcox/ds/d0k', $
		'/raid4/tcox/ds/d1k', $
		'/raid4/tcox/ds/d2k', $
		'/raid4/tcox/ds/d3k7', $
		'/raid4/tcox/ds/d4k', $
		'/raid4/tcox/ds/d5k', $
		'/raid4/tcox/ds/d6k']
	snapnum= [50,30,30,30,30,30,30]
	;fload_fruns_sigma_e, frun, snapnum, xvar3
	fload_fruns_sersicn, frun, snapnum, sersicn3, sersicn_1sig3, $
		excesslight3, exslgt_1sig3, $
		res3, res_1sig3

        ; 1 kpc/h is the minimum radius for fit
        ;sersicn3= [ 2.83, 2.43, 2.73, 3.84, 5.30, 5.50, 8.31]
        ;excesslight3= [ 22.8, 47.1, 41.9, 18.5, -2.05, -6.45, -39.8]
        ;res3= [ 1.2, 2.2, 2.9, 3.34, 3.7, 5.5, 11.2]

endif


;------------------------------------------------------------------------------------



; Print thie mess up
; -------------------
filename='profinfo.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12, newysize=20
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12, newysize=18

x0= 0.18 & x1= 0.99
;y0= 0.09 & ysize=0.3      ; 1 x 3  panel
y0= 0.09 & ysize=0.45       ; 1 x 2 panel
y1= y0+ysize
y2= y0+ysize+ysize
y3= y0+ysize+ysize+ysize
y4= y0+ysize+ysize+ysize+ysize
y5= y0+ysize+ysize+ysize+ysize+ysize


f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill



; zeroth: Re
; -------------
;yaxistitle='!6R!D!17e!6!N'
;ymax = 14.6
;ymin= 0.6
;!p.position= [x0, y2, x1, y3]
;plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
;        xrange=[xmin,xmax], yrange=[ymin,ymax], $
;        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
;        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
;
;oplot, xvar1, res1, psym=-2, color=150, thick= 3.0, symsize= 1.5
;oplot, xvar2, res2, psym=-8, color=50, thick= 3.0, symsize= 1.5
;oplot, xvar3, res3, psym=-5, color=100, thick= 3.0, symsize= 1.5
;
;oploterror, xvar1, res1, res_1sig1, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0
;oploterror, xvar2, res2, res_1sig2, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0
;oploterror, xvar3, res3, res_1sig3, psym=-3, errcolor=100, color=100, thick= 3.0, errthick= 3.0

; -------------

;oplot, [0.66,0.74], [4.15,4.15], psym=-2, color=150, thick= 3.0, symsize= 1.5
;oplot, [0.66,0.74], [3.75,3.75], psym=-8, color=50, thick= 3.0, symsize= 1.5
;oplot, [0.66,0.74], [3.35,3.35], psym=-5, color=100, thick= 3.0, symsize= 1.5

xyouts, 0.78, 0.95, 'orbit e', /normal, charthick=3.0, size=1.5, color=150
xyouts, 0.78, 0.92, 'orbit h', /normal, charthick=3.0, size=1.5, color=50
xyouts, 0.78, 0.89, 'orbit k', /normal, charthick=3.0, size=1.5, color=100

; --------------------------------------------------------------------------------------


; zeroth: Sersic n
; -------------
yaxistitle='!6Sersic !17n!6'
ymax = 10.8
ymin= 2.1
!p.position= [x0, y1, x1, y2]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, xvar1, sersicn1, psym=-2, color=150, thick= 3.0, symsize= 1.5
oplot, xvar2, sersicn2, psym=-8, color=50, thick= 3.0, symsize= 1.5
oplot, xvar3, sersicn3, psym=-5, color=100, thick= 3.0, symsize= 1.5

oploterror, xvar1, sersicn1, sersicn_1sig1, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0
oploterror, xvar2, sersicn2, sersicn_1sig2, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0
oploterror, xvar3, sersicn3, sersicn_1sig3, psym=-3, errcolor=100, color=100, thick= 3.0, errthick= 3.0

; --------------------------------------------------------------------------------------


; first: Excess Light
; -------------------
yaxistitle='!6Percent Extra Light'
ymax = 22.0
ymin= -22.0
!p.position= [x0, y0, x1, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase

oplot, xvar1, excesslight1, psym=-2, color=150, thick= 3.0, symsize= 1.5
oplot, xvar2, excesslight2, psym=-8, color=50, thick= 3.0, symsize= 1.5
oplot, xvar3, excesslight3, psym=-5, color=100, thick= 3.0, symsize= 1.5

oploterror, xvar1, excesslight1, exslgt_1sig1, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0
oploterror, xvar2, excesslight2, exslgt_1sig2, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0
oploterror, xvar3, excesslight3, exslgt_1sig3, psym=-3, errcolor=100, color=100, thick= 3.0, errthick= 3.0


oplot, [xmin,xmax], [0.0,0.0], psym=-3, color=0, thick=2.0, linestyle=1




; -------------
;  Done
; -------------

device, /close



end







;================================================================================





