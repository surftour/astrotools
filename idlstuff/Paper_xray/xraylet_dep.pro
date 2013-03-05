
; =================================================
;
;  various relations regarding hot gas
;
;
;
; ==================================================




; ----------------------------
;  Read wind.txt file
; ----------------------------
pro read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z

filename= frun+'/wind.txt'

spawn, "wc "+filename,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then filedata= fltarr(6,lines)

openr, 1, filename

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, filedata
close, 1


time= filedata[0,*]
mass_gtx= filedata[1,*]
mass_egt= filedata[2,*]
wind_ke= filedata[3,*]
wind_rke= filedata[4,*]
wind_z= filedata[5,*]


end



; -----------------------------------------------------------------




;   xray info
; ---------------
pro fload_xrayinfo, frun, xraygasmass, xraygasz, lx_max, lx_rem


   ns= n_elements(frun)
   allxmass= fltarr(ns)
   allxz= fltarr(ns)
   allxmax= fltarr(ns)
   allxrem= fltarr(ns)

   for i=0,ns-1 do begin

	read_xrays_file, frun[i], time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

        allxmass[i]= mass_xraygas(n_elements(time)-1)
        allxz[i]= z_X(n_elements(time)-1)
        allxmax[i]= max(xray_rs_s)
        allxrem[i]= xray_rs_s(n_elements(time)-1)
   endfor

   xraygasmass= allxmass
   xraygasz= allxz
   lx_max= allxmax
   lx_rem= allxrem


end




;   hotgas info
; --------------------------
pro fload_hotgasinfo, frun, init_gasmass, rem_gasmass

   ns= n_elements(frun)
   allinitm= fltarr(ns)
   allremm= fltarr(ns)

   for i=0,ns-1 do begin


   read_hotgas_file, frun[i], time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf


        allinitm[i]= gas_tot(0)
        allremm[i]= gas_tot(n_elements(time)-1)
   endfor

   init_gasmass= allinitm
   rem_gasmass= allremm

end





;   wind info
; --------------------------
pro fload_windinfo, frun, windmass, windz

   ns= n_elements(frun)
   allwmass= fltarr(ns)
   allwz= fltarr(ns)

   for i=0,ns-1 do begin

	read_wind_file, frun[i], time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
	allwmass[i]= [mass_egt(n_elements(time)-1)]
	allwz[i]= [wind_z(n_elements(time)-1)]
   endfor

   windmass= allwmass
   windz= allwz

end






;   average sigmas
; --------------------------
pro fload_fruns_sigma_e, frun, xvar

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






;============================================================================








;--------------------------------------------------------------------------------


;===================================
;
;  Relations  1
;
;===================================


;
pro r1, junk


if not keyword_set(junk) then begin
	print, " "
	print, " r1, junk"
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
        xvar1= [0.05,0.20,0.40,0.60,0.80,1.00]
        xvar2= xvar1
        xvar3= xvar1
        xmax= 1.1 & xmin= -0.1
        xaxistitle='Progenitor Gas Fraction'


        ; e
        ;
        frun= ['/raid4/tcox/gfs/vc3vc3z_e', $
                '/raid4/tcox/gfs/vc3vc3y_e', $
                '/raid4/tcox/gfs/vc3vc3x2_e', $
                '/raid4/tcox/gfs/vc3vc3w_e', $
                '/raid4/tcox/gfs/vc3vc3v_e', $
                '/raid4/tcox/gfs/vc3vc3u_e']

	fload_windinfo, frun, windmass1, windz1
	fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1
	fload_hotgasinfo, frun, init_gasmass1, rem_gasmass1

	gasconsumption1= rem_gasmass1/init_gasmass1

	;xraygm1= xraygasmass1/rem_gasmass1
	;windgm1= windmass1/rem_gasmass1
	xraygm1= xraygasmass1/init_gasmass1
	windgm1= windmass1/init_gasmass1


        ; h
        ;
        frun= ['/raid4/tcox/gfs/vc3vc3z_h', $
                '/raid4/tcox/gfs/vc3vc3y_h', $
                '/raid4/tcox/gfs/vc3vc3x2_h', $
                '/raid4/tcox/gfs/vc3vc3w_h', $
                '/raid4/tcox/gfs/vc3vc3v_h', $
                '/raid4/tcox/gfs/vc3vc3u_h']

	fload_windinfo, frun, windmass2, windz2
	fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2
	fload_hotgasinfo, frun, init_gasmass2, rem_gasmass2

	gasconsumption2= rem_gasmass2/init_gasmass2

	;xraygm2= xraygasmass2/rem_gasmass2
	;windgm2= windmass2/rem_gasmass2
	xraygm2= xraygasmass2/init_gasmass2
	windgm2= windmass2/init_gasmass2



        ; k
        ;
        frun= ['/raid4/tcox/gfs/vc3vc3z_k', $
                '/raid4/tcox/gfs/vc3vc3y_k', $
                '/raid4/tcox/gfs/vc3vc3x2_k', $
                '/raid4/tcox/gfs/vc3vc3w_k', $
                '/raid4/tcox/gfs/vc3vc3v_k', $
                '/raid4/tcox/gfs/vc3vc3u_k']

	fload_windinfo, frun, windmass3, windz3
	fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3
	fload_hotgasinfo, frun, init_gasmass3, rem_gasmass3

	gasconsumption3= rem_gasmass3/init_gasmass3

	;xraygm3= xraygasmass3/rem_gasmass3
	;windgm3= windmass3/rem_gasmass3
	xraygm3= xraygasmass3/init_gasmass3
	windgm3= windmass3/init_gasmass3


endif





;  size
; -------------
do_size= 0
;do_size= 1
if do_size eq 1 then begin
        xvar1= [50,80.0,115.0,160.0,225.0,320.0,500.0]
        xvar2= xvar1
        xvar3= xvar1
        ;xmax= 580.0 & xmin= 20.0            ; linear scale
        xmax= 800.0 & xmin= 40.0             ; log scale
        xaxistitle='!7r!6 (km sec!E-1!N)'


        ; e
        ;
        frun= ['/raid4/tcox/ds/d0e', $
                '/raid4/tcox/ds/d1e', $
                '/raid4/tcox/ds/d2e', $
         ;       '/raid4/tcox/ds/d3e7', $
         ;       '/raid4/tcox/vc3vc3e_2', $
                '/raid4/tcox/gfs/vc3vc3x2_e', $
                '/raid4/tcox/ds/d4e', $
                '/raid4/tcox/ds/d5e', $
                '/raid4/tcox/ds/d6e']

	fload_windinfo, frun, windmass1, windz1
	fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1
	fload_hotgasinfo, frun, init_gasmass1, rem_gasmass1
	fload_fruns_sigma_e, frun, xvar1

	gasconsumption1= rem_gasmass1/init_gasmass1

	;xraygm1= xraygasmass1/rem_gasmass1
	;windgm1= windmass1/rem_gasmass1
	xraygm1= xraygasmass1/init_gasmass1
	windgm1= windmass1/init_gasmass1



        ; h
        ;
        frun= ['/raid4/tcox/ds/d0h', $
                '/raid4/tcox/ds/d1h', $
                '/raid4/tcox/ds/d2h', $
        ;        '/raid4/tcox/ds/d3h7', $
                '/raid4/tcox/gfs/vc3vc3x2_h', $
                '/raid4/tcox/ds/d4h', $
                '/raid4/tcox/ds/d5h', $
                '/raid4/tcox/ds/d6h']

	fload_windinfo, frun, windmass2, windz2
	fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2
	fload_hotgasinfo, frun, init_gasmass2, rem_gasmass2
	fload_fruns_sigma_e, frun, xvar2

	gasconsumption2= rem_gasmass2/init_gasmass2

	;xraygm2= xraygasmass2/rem_gasmass2
	;windgm2= windmass2/rem_gasmass2
	xraygm2= xraygasmass2/init_gasmass2
	windgm2= windmass2/init_gasmass2




        ; k
        ;
        frun= ['/raid4/tcox/ds/d0k', $
                '/raid4/tcox/ds/d1k', $
                '/raid4/tcox/ds/d2k', $
        ;        '/raid4/tcox/ds/d3k7', $
                '/raid4/tcox/gfs/vc3vc3x2_k', $
                '/raid4/tcox/ds/d4k', $
                '/raid4/tcox/ds/d5k', $
                '/raid4/tcox/ds/d6k']

	fload_windinfo, frun, windmass3, windz3
	fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3
	fload_hotgasinfo, frun, init_gasmass3, rem_gasmass3
	fload_fruns_sigma_e, frun, xvar3

	gasconsumption3= rem_gasmass3/init_gasmass3

	;xraygm3= xraygasmass3/rem_gasmass3
	;windgm3= windmass3/rem_gasmass3
	xraygm3= xraygasmass3/init_gasmass3
	windgm3= windmass3/init_gasmass3


endif




;  size II
; -------------
do_size= 0
;do_size= 1
if do_size eq 1 then begin
        xvar1= [50,80.0,115.0,160.0,225.0,320.0,500.0]
        xvar2= xvar1
        xvar3= xvar1
        xmax= 580.0 & xmin= 20.0
        xaxistitle='!7r!6 (km sec!E-1!N)'


        ; e
        ;
        frun= ['/raid4/tcox/As/A0', $
                '/raid4/tcox/As/A1', $
                '/raid4/tcox/As/A2', $
                '/raid4/tcox/As/A3', $
                '/raid4/tcox/As/A4', $
                '/raid4/tcox/As/A5', $
                '/raid4/tcox/As/A6']

	fload_windinfo, frun, windmass1, windz1
	fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1
	fload_hotgasinfo, frun, init_gasmass1, rem_gasmass1
	fload_fruns_sigma_e, frun, xvar1

	gasconsumption1= rem_gasmass1/init_gasmass1

	;xraygm1= xraygasmass1/rem_gasmass1
	;windgm1= windmass1/rem_gasmass1
	xraygm1= xraygasmass1/init_gasmass1
	windgm1= windmass1/init_gasmass1


	xraygm= xraygm1
	windgm= windgm1
	xraygasz= xraygasz1
	windz= windz1


endif



; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(3)
   n_f= n_elements(xvar1)
   xraygm= fltarr(n_f)
   xraygm_err= fltarr(n_f)
   windgm= fltarr(n_f)
   windgm_err= fltarr(n_f)
   xraygasz= fltarr(n_f)
   xraygasz_err= fltarr(n_f)
   windz= fltarr(n_f)
   windz_err= fltarr(n_f)

   for i=0,n_f-1 do begin
	vars[0]= xraygm1[i]
	vars[1]= xraygm2[i]
	vars[2]= xraygm3[i]
	xraygm[i]= mean(vars)
	xraygm_err[i]= sqrt(variance(vars))

	vars[0]= windgm1[i]
	vars[1]= windgm2[i]
	vars[2]= windgm3[i]
	windgm[i]= mean(vars)
	windgm_err[i]= sqrt(variance(vars))

	vars[0]= xraygasz1[i]
	vars[1]= xraygasz2[i]
	vars[2]= xraygasz3[i]
	xraygasz[i]= mean(vars)
	xraygasz_err[i]= sqrt(variance(vars))

	vars[0]= windz1[i]
	vars[1]= windz2[i]
	vars[2]= windz3[i]
	windz[i]= mean(vars)
	windz_err[i]= sqrt(variance(vars))

   endfor
endif



;------------------------------------------------------------------------------------



; Print thie mess up
; -------------------
filename='hotgasmasses.eps'

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




; --------------------------------------------------------------------------------------


; zeroth: Mass fractions
; ---------------------------
yaxistitle='!6f'
ymax = 0.44
;ymax = 0.25
ymin= 0.0
!p.position= [x0, y1, x1, y2]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $  ;/xlog, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

;oplot, xvar1, xraygm1, psym=-2, color=150, thick= 3.0, symsize= 1.5
;oplot, xvar1, windgm1, psym=-2, color=150, thick= 3.0, symsize= 1.5

;oplot, xvar2, xraygm2, psym=-8, color=50, thick= 3.0, symsize= 1.5
;oplot, xvar2, windgm2, psym=-8, color=50, thick= 3.0, symsize= 1.5

;oplot, xvar3, xraygm3, psym=-5, color=100, thick= 3.0, symsize= 1.5
;oplot, xvar3, windgm3, psym=-5, color=100, thick= 3.0, symsize= 1.5

oplot, xvar1, xraygm, psym=-2, color=150, thick= 3.0, symsize= 1.5
oploterror, xvar1, xraygm, xraygm_err, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0

oplot, xvar1, windgm, psym=-8, color=50, thick= 3.0, symsize= 1.5
oploterror, xvar1, windgm, windgm_err, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0


; plot the sum as a dashed line
oplot, xvar1, xraygm+windgm, psym=-3, color=0, thick= 2.0, linestyle= 1






; --------------------------------------------------------------------------------------



; first: Excess Light
; -------------------
yaxistitle='!6<Z> (Z!D!9n!6!N)'
ymax = 3.0
ymin= 0.05
!p.position= [x0, y0, x1, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
	/ylog, $ ;/xlog, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase

;oplot, xvar1, xraygasz1, psym=-2, color=150, thick= 3.0, symsize= 1.5
;oplot, xvar1, windz1, psym=-2, color=150, thick= 3.0, symsize= 1.5

;oplot, xvar2, xraygasz2, psym=-8, color=50, thick= 3.0, symsize= 1.5
;oplot, xvar2, windz2, psym=-8, color=50, thick= 3.0, symsize= 1.5

;oplot, xvar3, xraygasz3, psym=-5, color=100, thick= 3.0, symsize= 1.5
;oplot, xvar3, windz3, psym=-5, color=100, thick= 3.0, symsize= 1.5

oplot, xvar1, xraygasz, psym=-2, color=150, thick= 3.0, symsize= 1.5
oploterror, xvar1, xraygasz, xraygasz_err, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0

oplot, xvar1, windz, psym=-8, color=50, thick= 3.0, symsize= 1.5
oploterror, xvar1, windz, windz_err, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0



;-----------
xyouts, 0.68, 0.25, 'X-ray', /normal, charthick=3.0, size=1.5, color=150
xyouts, 0.68, 0.18, 'Unbound', /normal, charthick=3.0, size=1.5, color=50

; sigma - log
;oplot, [55.0,70.0], [0.24,0.24], psym=-2, color=150, thick= 3.0, symsize= 1.5
;oplot, [55.0,70.0], [0.125,0.125], psym=-8, color=50, thick= 3.0, symsize= 1.5

; gas frac - linear
oplot, [0.4,0.48], [0.24,0.24], psym=-2, color=150, thick= 3.0, symsize= 1.5
oplot, [0.4,0.48], [0.125,0.125], psym=-8, color=50, thick= 3.0, symsize= 1.5




;---------
device, /close

end











;================================================================================================










;===================================
;
;  Relations  3
;
;===================================


;
pro r3, junk


if not keyword_set(junk) then begin
	print, " "
	print, " r3, junk"
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
        xvar1= [0.05,0.20,0.40,0.60,0.80,1.00]
        xvar2= xvar1
        xvar3= xvar1
        xmax= 1.1 & xmin= -0.1
        xaxistitle='Progenitor Gas Fraction'


        ; e
        ;
        frun= ['/raid4/tcox/gfs/vc3vc3z_e', $
                '/raid4/tcox/gfs/vc3vc3y_e', $
                '/raid4/tcox/gfs/vc3vc3x2_e', $
                '/raid4/tcox/gfs/vc3vc3w_e', $
                '/raid4/tcox/gfs/vc3vc3v_e', $
                '/raid4/tcox/gfs/vc3vc3u_e']

	fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1


        ; h
        ;
        frun= ['/raid4/tcox/gfs/vc3vc3z_h', $
                '/raid4/tcox/gfs/vc3vc3y_h', $
                '/raid4/tcox/gfs/vc3vc3x2_h', $
                '/raid4/tcox/gfs/vc3vc3w_h', $
                '/raid4/tcox/gfs/vc3vc3v_h', $
                '/raid4/tcox/gfs/vc3vc3u_h']

	fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2



        ; k
        ;
        frun= ['/raid4/tcox/gfs/vc3vc3z_k', $
                '/raid4/tcox/gfs/vc3vc3y_k', $
                '/raid4/tcox/gfs/vc3vc3x2_k', $
                '/raid4/tcox/gfs/vc3vc3w_k', $
                '/raid4/tcox/gfs/vc3vc3v_k', $
                '/raid4/tcox/gfs/vc3vc3u_k']

	fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3


endif





;  size
; -------------
do_size= 0
;do_size= 1
if do_size eq 1 then begin
        xvar1= [50,80.0,115.0,160.0,225.0,320.0,500.0]
        xvar2= xvar1
        xvar3= xvar1
        xmax= 580.0 & xmin= 20.0
        ;xmax= 800.0 & xmin= 20.0
        xaxistitle='!7r!6 (km sec!E-1!N)'


        ; e
        ;
        frun= ['/raid4/tcox/ds/d0e', $
                '/raid4/tcox/ds/d1e', $
                '/raid4/tcox/ds/d2e', $
         ;       '/raid4/tcox/ds/d3e7', $
         ;       '/raid4/tcox/vc3vc3e_2', $
                '/raid4/tcox/gfs/vc3vc3x2_e', $
                '/raid4/tcox/ds/d4e', $
                '/raid4/tcox/ds/d5e', $
                '/raid4/tcox/ds/d6e']

	fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1
	fload_fruns_sigma_e, frun, xvar1



        ; h
        ;
        frun= ['/raid4/tcox/ds/d0h', $
                '/raid4/tcox/ds/d1h', $
                '/raid4/tcox/ds/d2h', $
        ;        '/raid4/tcox/ds/d3h7', $
                '/raid4/tcox/gfs/vc3vc3x2_h', $
                '/raid4/tcox/ds/d4h', $
                '/raid4/tcox/ds/d5h', $
                '/raid4/tcox/ds/d6h']

	fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2
	fload_fruns_sigma_e, frun, xvar2



        ; k
        ;
        frun= ['/raid4/tcox/ds/d0k', $
                '/raid4/tcox/ds/d1k', $
                '/raid4/tcox/ds/d2k', $
        ;        '/raid4/tcox/ds/d3k7', $
                '/raid4/tcox/gfs/vc3vc3x2_k', $
                '/raid4/tcox/ds/d4k', $
                '/raid4/tcox/ds/d5k', $
                '/raid4/tcox/ds/d6k']

	fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3
	fload_fruns_sigma_e, frun, xvar3


endif




; average the runs
; ------------------
average_three_runs= 0
;average_three_runs= 1
if average_three_runs eq 1 then begin
   vars= fltarr(3)
   n_f= n_elements(xvar1)
   lx_max= fltarr(n_f)
   lx_max_err= fltarr(n_f)
   lx_rem= fltarr(n_f)
   lx_rem_err= fltarr(n_f)
   sigma= fltarr(n_f)
   sigma_err= fltarr(n_f)

   for i=0,n_f-1 do begin
	vars[0]= lx_max1[i]
	vars[1]= lx_max2[i]
	vars[2]= lx_max3[i]
	lx_max[i]= mean(vars)
	lx_max_err[i]= sqrt(variance(vars))

	vars[0]= lx_rem1[i]
	vars[1]= lx_rem2[i]
	vars[2]= lx_rem3[i]
	lx_rem[i]= mean(vars)
	lx_rem_err[i]= sqrt(variance(vars))

	vars[0]= xvar1[i]
	vars[1]= xvar2[i]
	vars[2]= xvar3[i]
	sigma[i]= mean(vars)
	sigma_err[i]= sqrt(variance(vars))

   endfor

   xvar= sigma
endif



;------------------------------------------------------------------------------------



; Print thie mess up
; -------------------
filename='xrays.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12, newysize=20
;setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12, newysize=18
setup_plot_stuff, 'ps', filename=filename, colortable= 4

x0= 0.18 & x1= 0.98
y0= 0.15 & ysize=0.83      ; 1 x 1  std figure
;y0= 0.09 & ysize=0.45       ; 1 x 2 panel
y1= y0+ysize
y2= y0+ysize+ysize
y3= y0+ysize+ysize+ysize
y4= y0+ysize+ysize+ysize+ysize
y5= y0+ysize+ysize+ysize+ysize+ysize


f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill



xyouts, 0.68, 0.35, 'Remnant', /normal, charthick=3.0, size=1.5, color=150
xyouts, 0.68, 0.29, 'Maximum', /normal, charthick=3.0, size=1.5, color=50

; --------------------------------------------------------------------------------------


yaxistitle='!6L!DX!N (ergs s!E-1!N)'
ymax = 44.0
ymin= 36.0
!p.position= [x0, y0, x1, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $   ; /xlog, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase

;oplot, xvar1, lx_max1, psym=-2, color=150, thick= 3.0, symsize= 1.5
;oplot, xvar1, lx_rem1, psym=-2, color=150, thick= 3.0, symsize= 1.5

;oplot, xvar2, lx_max2, psym=-8, color=50, thick= 3.0, symsize= 1.5
;oplot, xvar2, lx_rem2, psym=-8, color=50, thick= 3.0, symsize= 1.5

;oplot, xvar3, lx_max3, psym=-5, color=100, thick= 3.0, symsize= 1.5
;oplot, xvar3, lx_rem3, psym=-5, color=100, thick= 3.0, symsize= 1.5


oplot, xvar, lx_rem, psym=-2, color=150, thick= 3.0, symsize= 1.5
oploterror, xvar, lx_rem, lx_rem_err, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0

oplot, xvar, lx_max, psym=-8, color=50, thick= 3.0, symsize= 1.5
oploterror, xvar, lx_max, lx_max_err, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0








;---------
device, /close

end











;================================================================================================




;===================================
;
;  Relations  All-in-one
;
;===================================


;
pro r5, junk


if not keyword_set(junk) then begin
	print, " "
	print, " r5, junk"
	print, " "
	print, " "
	return
endif



;------------------------------------------------------------------------------------


xvar1= [50,80.0,115.0,160.0,225.0,320.0,500.0]
xvar2= xvar1
xvar3= xvar1
;xmax= 580.0 & xmin= 20.0            ; linear scale
xmax= 800.0 & xmin= 40.0             ; log scale
xaxistitle='!7r!6 (km sec!E-1!N)'



;------------------------------------------------------------------------------------


print, '5 percent gas '

;
;  5 % Gas
;----------------

; e
;
frun= ['/raid4/tcox/zs/z0e', '/raid4/tcox/zs/z1e', '/raid4/tcox/zs/z2e', '/raid4/tcox/zs/z3e', $
		'/raid4/tcox/zs/z4e', '/raid4/tcox/zs/z5e', '/raid4/tcox/zs/z6e']

fload_windinfo, frun, windmass1, windz1
fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1
fload_hotgasinfo, frun, init_gasmass1, rem_gasmass1
;fload_fruns_sigma_e, frun, xvar1

gasconsumption1= rem_gasmass1/init_gasmass1
xraygm1= xraygasmass1/init_gasmass1
windgm1= windmass1/init_gasmass1



; h
;
frun= ['/raid4/tcox/zs/z0h', '/raid4/tcox/zs/z1h', '/raid4/tcox/zs/z2h', '/raid4/tcox/zs/z3h', $
		'/raid4/tcox/zs/z4h', '/raid4/tcox/zs/z5h', '/raid4/tcox/zs/z6h']

fload_windinfo, frun, windmass2, windz2
fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2
fload_hotgasinfo, frun, init_gasmass2, rem_gasmass2
;fload_fruns_sigma_e, frun, xvar2

gasconsumption2= rem_gasmass2/init_gasmass2

xraygm2= xraygasmass2/init_gasmass2
windgm2= windmass2/init_gasmass2




; k
;
frun= ['/raid4/tcox/zs/z0k', '/raid4/tcox/zs/z1k', '/raid4/tcox/zs/z2k', '/raid4/tcox/zs/z3k', $
		'/raid4/tcox/zs/z4k', '/raid4/tcox/zs/z5k', '/raid4/tcox/zs/z6k']

fload_windinfo, frun, windmass3, windz3
fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3
fload_hotgasinfo, frun, init_gasmass3, rem_gasmass3
;fload_fruns_sigma_e, frun, xvar3

gasconsumption3= rem_gasmass3/init_gasmass3

xraygm3= xraygasmass3/init_gasmass3
windgm3= windmass3/init_gasmass3



; f
;
frun= ['/raid4/tcox/zs/z0f', '/raid4/tcox/zs/z1f', '/raid4/tcox/zs/z2f', '/raid4/tcox/zs/z3f', $
                '/raid4/tcox/zs/z4f', '/raid4/tcox/zs/z5f', '/raid4/tcox/zs/z6f']

fload_windinfo, frun, windmass4, windz4
fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4
fload_hotgasinfo, frun, init_gasmass4, rem_gasmass4
;fload_fruns_sigma_e, frun, xvar4

gasconsumption4= rem_gasmass4/init_gasmass4

xraygm4= xraygasmass4/init_gasmass4
windgm4= windmass4/init_gasmass4



; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(xvar1)
   xraygm_5= fltarr(n_f)
   xraygm_err_5= fltarr(n_f)
   windgm_5= fltarr(n_f)
   windgm_err_5= fltarr(n_f)
   xraygasz_5= fltarr(n_f)
   xraygasz_err_5= fltarr(n_f)
   windz_5= fltarr(n_f)
   windz_err_5= fltarr(n_f)
   lx_max_5= fltarr(n_f)
   lx_max_err_5= fltarr(n_f)
   lx_rem_5= fltarr(n_f)
   lx_rem_err_5= fltarr(n_f)

   for i=0,n_f-1 do begin
	vars[0]= xraygm1[i]
	vars[1]= xraygm2[i]
	vars[2]= xraygm3[i]
	vars[3]= xraygm4[i]
	xraygm_5[i]= mean(vars)
	xraygm_err_5[i]= sqrt(variance(vars))

	vars[0]= windgm1[i]
	vars[1]= windgm2[i]
	vars[2]= windgm3[i]
	vars[3]= windgm4[i]
	windgm_5[i]= mean(vars)
	windgm_err_5[i]= sqrt(variance(vars))

	vars[0]= xraygasz1[i]
	vars[1]= xraygasz2[i]
	vars[2]= xraygasz3[i]
	vars[3]= xraygasz4[i]
	xraygasz_5[i]= mean(vars)
	xraygasz_err_5[i]= sqrt(variance(vars))

	vars[0]= windz1[i]
	vars[1]= windz2[i]
	vars[2]= windz3[i]
	vars[3]= windz4[i]
	windz_5[i]= mean(vars)
	windz_err_5[i]= sqrt(variance(vars))

	vars[0]= lx_max1[i]
	vars[1]= lx_max2[i]
	vars[2]= lx_max3[i]
	vars[3]= lx_max4[i]
	lx_max_5[i]= mean(vars)
	lx_max_err_5[i]= sqrt(variance(vars))

	vars[0]= lx_rem1[i]
	vars[1]= lx_rem2[i]
	vars[2]= lx_rem3[i]
	vars[3]= lx_rem4[i]
	lx_rem_5[i]= mean(vars)
	lx_rem_err_5[i]= sqrt(variance(vars))

	;vars[0]= xvar1[i]
	;vars[1]= xvar2[i]
	;vars[2]= xvar3[i]
	;vars[3]= xvar4[i]
	;sigma_5[i]= mean(vars)
	;sigma_err_5[i]= sqrt(variance(vars))

   endfor
endif





; -----------------------------------------------------------------------


print, '20 percent gas '

;
;  20 % Gas
;----------------

; e
;
;frun= ['/raid4/tcox/vcs/vc0vc0e', '/raid4/tcox/vcs/vc1vc1e', '/raid4/tcox/vcs/vc2vc2e', '/raid4/tcox/vcs/vc3vc3_e', $
;		'/raid4/tcox/vcs/vc4vc4e', '/raid4/tcox/vcs/vc5vc5e', '/raid4/tcox/vcs/vc6vc6e']
frun= ['/raid4/tcox/es/e0e', '/raid4/tcox/es/e1e', '/raid4/tcox/es/e2e', '/raid4/tcox/es/e3e', $
		'/raid4/tcox/es/e4e', '/raid4/tcox/es/e5e', '/raid4/tcox/es/e6e']

fload_windinfo, frun, windmass1, windz1
fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1
fload_hotgasinfo, frun, init_gasmass1, rem_gasmass1
;fload_fruns_sigma_e, frun, xvar1

gasconsumption1= rem_gasmass1/init_gasmass1
xraygm1= xraygasmass1/init_gasmass1
windgm1= windmass1/init_gasmass1



; h
;
;frun= ['/raid4/tcox/vcs/vc0vc0h', '/raid4/tcox/vcs/vc1vc1', '/raid4/tcox/vcs/vc2vc2', '/raid4/tcox/vcs/vc3vc3', $
;		'/raid4/tcox/vcs/vc4vc4a', '/raid4/tcox/vcs/vc5vc5a', '/raid4/tcox/vcs/vc6vc6a']
frun= ['/raid4/tcox/es/e0h', '/raid4/tcox/es/e1h', '/raid4/tcox/es/e2h', '/raid4/tcox/es/e3h', $
		'/raid4/tcox/es/e4h', '/raid4/tcox/es/e5h', '/raid4/tcox/es/e6h']

fload_windinfo, frun, windmass2, windz2
fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2
fload_hotgasinfo, frun, init_gasmass2, rem_gasmass2
;fload_fruns_sigma_e, frun, xvar2

gasconsumption2= rem_gasmass2/init_gasmass2

xraygm2= xraygasmass2/init_gasmass2
windgm2= windmass2/init_gasmass2



; f
;
frun= ['/raid4/tcox/es/e0f', '/raid4/tcox/es/e1f', '/raid4/tcox/es/e2f', '/raid4/tcox/es/e3f', $
                '/raid4/tcox/es/e4f', '/raid4/tcox/es/e5f', '/raid4/tcox/es/e6f']

fload_windinfo, frun, windmass3, windz3
fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3
fload_hotgasinfo, frun, init_gasmass3, rem_gasmass3
;fload_fruns_sigma_e, frun, xvar3

gasconsumption3= rem_gasmass3/init_gasmass3

xraygm3= xraygasmass3/init_gasmass3
windgm3= windmass3/init_gasmass2



; k
;
frun= ['/raid4/tcox/es/e0k', '/raid4/tcox/es/e1k', '/raid4/tcox/es/e2k', '/raid4/tcox/es/e3k', $
                '/raid4/tcox/es/e4k', '/raid4/tcox/es/e5k', '/raid4/tcox/es/e6k']

fload_windinfo, frun, windmass4, windz4
fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4
fload_hotgasinfo, frun, init_gasmass4, rem_gasmass4
;fload_fruns_sigma_e, frun, xvar4

gasconsumption4= rem_gasmass4/init_gasmass4

xraygm4= xraygasmass4/init_gasmass4
windgm4= windmass4/init_gasmass4





; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(xvar1)
   xraygm_20= fltarr(n_f)
   xraygm_err_20= fltarr(n_f)
   windgm_20= fltarr(n_f)
   windgm_err_20= fltarr(n_f)
   xraygasz_20= fltarr(n_f)
   xraygasz_err_20= fltarr(n_f)
   windz_20= fltarr(n_f)
   windz_err_20= fltarr(n_f)
   lx_max_20= fltarr(n_f)
   lx_max_err_20= fltarr(n_f)
   lx_rem_20= fltarr(n_f)
   lx_rem_err_20= fltarr(n_f)

   for i=0,n_f-1 do begin
	vars[0]= xraygm1[i]
	vars[1]= xraygm2[i]
	vars[2]= xraygm3[i]
	vars[3]= xraygm4[i]
	xraygm_20[i]= mean(vars)
	xraygm_err_20[i]= sqrt(variance(vars))

	vars[0]= windgm1[i]
	vars[1]= windgm2[i]
	vars[2]= windgm3[i]
	vars[3]= windgm4[i]
	windgm_20[i]= mean(vars)
	windgm_err_20[i]= sqrt(variance(vars))

	vars[0]= xraygasz1[i]
	vars[1]= xraygasz2[i]
	vars[2]= xraygasz3[i]
	vars[3]= xraygasz4[i]
	xraygasz_20[i]= mean(vars)
	xraygasz_err_20[i]= sqrt(variance(vars))

	vars[0]= windz1[i]
	vars[1]= windz2[i]
	vars[2]= windz3[i]
	vars[3]= windz4[i]
	windz_20[i]= mean(vars)
	windz_err_20[i]= sqrt(variance(vars))

	vars[0]= lx_max1[i]
	vars[1]= lx_max2[i]
	vars[2]= lx_max3[i]
	vars[3]= lx_max4[i]
	lx_max_20[i]= mean(vars)
	lx_max_err_20[i]= sqrt(variance(vars))

	vars[0]= lx_rem1[i]
	vars[1]= lx_rem2[i]
	vars[2]= lx_rem3[i]
	vars[3]= lx_rem4[i]
	lx_rem_20[i]= mean(vars)
	lx_rem_err_20[i]= sqrt(variance(vars))

	;vars[0]= xvar1[i]
	;vars[1]= xvar2[i]
	;sigma_20[i]= mean(vars)
	;sigma_err_20[i]= sqrt(variance(vars))

   endfor
endif




; -----------------------------------------------------------------------


print, '40 percent gas '

;
;  40 % Gas
;----------------

; e
;
;frun= ['/raid4/tcox/ds/d0e', '/raid4/tcox/ds/d1e', '/raid4/tcox/ds/d2e', '/raid4/tcox/gfs/vc3vc3x2_e', $
;		'/raid4/tcox/ds/d4e', '/raid4/tcox/ds/d5e', '/raid4/tcox/ds/d6e']
frun= ['/raid4/tcox/ds/d0e2_q', '/raid4/tcox/ds/d1e2_q', '/raid4/tcox/ds/d2e2_q', '/raid4/tcox/ds/d3e7', $
		'/raid4/tcox/ds/d4e2_q', '/raid4/tcox/ds/d5e2_q', '/raid4/tcox/ds/d6e2_q']

fload_windinfo, frun, windmass1, windz1
fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1
fload_hotgasinfo, frun, init_gasmass1, rem_gasmass1
fload_fruns_sigma_e, frun, xvar1

gasconsumption1= rem_gasmass1/init_gasmass1
xraygm1= xraygasmass1/init_gasmass1
windgm1= windmass1/init_gasmass1



; h
;
;frun= ['/raid4/tcox/ds/d0h', '/raid4/tcox/ds/d1h', '/raid4/tcox/ds/d2h', '/raid4/tcox/gfs/vc3vc3x2_h', $
;		'/raid4/tcox/ds/d4h', '/raid4/tcox/ds/d5h', '/raid4/tcox/ds/d6h']
frun= ['/raid4/tcox/ds/d0h2_q', '/raid4/tcox/ds/d1h2_q', '/raid4/tcox/ds/d2h2_q', '/raid4/tcox/ds/d3h7', $
		'/raid4/tcox/ds/d4h2_q', '/raid4/tcox/ds/d5h2_q', '/raid4/tcox/ds/d6h2_q']

fload_windinfo, frun, windmass2, windz2
fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2
fload_hotgasinfo, frun, init_gasmass2, rem_gasmass2
fload_fruns_sigma_e, frun, xvar2

gasconsumption2= rem_gasmass2/init_gasmass2

xraygm2= xraygasmass2/init_gasmass2
windgm2= windmass2/init_gasmass2




; k
;
;frun= ['/raid4/tcox/ds/d0k', '/raid4/tcox/ds/d1k', '/raid4/tcox/ds/d2k', '/raid4/tcox/gfs/vc3vc3x2_k', $
;		'/raid4/tcox/ds/d4k', '/raid4/tcox/ds/d5k', '/raid4/tcox/ds/d6k']
frun= ['/raid4/tcox/ds/d0k2_q', '/raid4/tcox/ds/d1k2_q', '/raid4/tcox/ds/d2k2_q', '/raid4/tcox/ds/d3k7', $
		'/raid4/tcox/ds/d4k2_q', '/raid4/tcox/ds/d5k2_q', '/raid4/tcox/ds/d6k2_q']

fload_windinfo, frun, windmass3, windz3
fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3
fload_hotgasinfo, frun, init_gasmass3, rem_gasmass3
fload_fruns_sigma_e, frun, xvar3

gasconsumption3= rem_gasmass3/init_gasmass3

xraygm3= xraygasmass3/init_gasmass3
windgm3= windmass3/init_gasmass3



; f
;
frun= ['/raid4/tcox/ds/d0f2_q', '/raid4/tcox/ds/d1f2_q', '/raid4/tcox/ds/d2f2_q', '/raid4/tcox/ds/d3f7', $
                '/raid4/tcox/ds/d4f2_q', '/raid4/tcox/ds/d5f2_q', '/raid4/tcox/ds/d6f2_q']

fload_windinfo, frun, windmass4, windz4
fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4
fload_hotgasinfo, frun, init_gasmass4, rem_gasmass4
fload_fruns_sigma_e, frun, xvar4

gasconsumption4= rem_gasmass4/init_gasmass4

xraygm4= xraygasmass4/init_gasmass4
windgm4= windmass4/init_gasmass4




; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(xvar1)
   xraygm_40= fltarr(n_f)
   xraygm_err_40= fltarr(n_f)
   windgm_40= fltarr(n_f)
   windgm_err_40= fltarr(n_f)
   xraygasz_40= fltarr(n_f)
   xraygasz_err_40= fltarr(n_f)
   windz_40= fltarr(n_f)
   windz_err_40= fltarr(n_f)
   lx_max_40= fltarr(n_f)
   lx_max_err_40= fltarr(n_f)
   lx_rem_40= fltarr(n_f)
   lx_rem_err_40= fltarr(n_f)

   for i=0,n_f-1 do begin
	vars[0]= xraygm1[i]
	vars[1]= xraygm2[i]
	vars[2]= xraygm3[i]
	vars[3]= xraygm4[i]
	xraygm_40[i]= mean(vars)
	xraygm_err_40[i]= sqrt(variance(vars))

	vars[0]= windgm1[i]
	vars[1]= windgm2[i]
	vars[2]= windgm3[i]
	vars[3]= windgm4[i]
	windgm_40[i]= mean(vars)
	windgm_err_40[i]= sqrt(variance(vars))

	vars[0]= xraygasz1[i]
	vars[1]= xraygasz2[i]
	vars[2]= xraygasz3[i]
	vars[3]= xraygasz4[i]
	xraygasz_40[i]= mean(vars)
	xraygasz_err_40[i]= sqrt(variance(vars))

	vars[0]= windz1[i]
	vars[1]= windz2[i]
	vars[2]= windz3[i]
	vars[3]= windz4[i]
	windz_40[i]= mean(vars)
	windz_err_40[i]= sqrt(variance(vars))

	vars[0]= lx_max1[i]
	vars[1]= lx_max2[i]
	vars[2]= lx_max3[i]
	vars[3]= lx_max4[i]
	lx_max_40[i]= mean(vars)
	lx_max_err_40[i]= sqrt(variance(vars))

	vars[0]= lx_rem1[i]
	vars[1]= lx_rem2[i]
	vars[2]= lx_rem3[i]
	vars[3]= lx_rem4[i]
	lx_rem_40[i]= mean(vars)
	lx_rem_err_40[i]= sqrt(variance(vars))

	;vars[0]= xvar1[i]
	;vars[1]= xvar2[i]
	;vars[2]= xvar3[i]
	;sigma_40[i]= mean(vars)
	;sigma_err_40[i]= sqrt(variance(vars))

   endfor
endif





; -----------------------------------------------------------------------


;print, '100 percent gas '
print, '80 percent gas '

;
;  100 % Gas
; -------------
xvar1= [50,80.0,115.0,160.0,225.0,320.0,500.0]
xvar2= xvar1
xvar3= xvar1


; e
;
;frun= ['/raid4/tcox/As/A0e', '/raid4/tcox/As/A1e', '/raid4/tcox/As/A2e', '/raid4/tcox/As/A3e', $
;                '/raid4/tcox/As/A4e', '/raid4/tcox/As/A5e', '/raid4/tcox/As/A6e']
frun= ['/raid4/tcox/bs/b0e', '/raid4/tcox/bs/b1e', '/raid4/tcox/bs/b2e', '/raid4/tcox/bs/b3e', $
                '/raid4/tcox/bs/b4e', '/raid4/tcox/bs/b5e', '/raid4/tcox/bs/b6e']

fload_windinfo, frun, windmass1, windz1
fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1
fload_hotgasinfo, frun, init_gasmass1, rem_gasmass1
;fload_fruns_sigma_e, frun, xvar1

gasconsumption1= rem_gasmass1/init_gasmass1

xraygm1= xraygasmass1/init_gasmass1
windgm1= windmass1/init_gasmass1



; h
;
;frun= ['/raid4/tcox/As/A0', '/raid4/tcox/As/A1', '/raid4/tcox/As/A2', '/raid4/tcox/As/A3', $
;                '/raid4/tcox/As/A4', '/raid4/tcox/As/A5', '/raid4/tcox/As/A6']
frun= ['/raid4/tcox/bs/b0h', '/raid4/tcox/bs/b1h', '/raid4/tcox/bs/b2h', '/raid4/tcox/bs/b3h', $
                '/raid4/tcox/bs/b4h', '/raid4/tcox/bs/b5h', '/raid4/tcox/bs/b6h']

fload_windinfo, frun, windmass2, windz2
fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2
fload_hotgasinfo, frun, init_gasmass2, rem_gasmass2
;fload_fruns_sigma_e, frun, xvar2

gasconsumption2= rem_gasmass2/init_gasmass2

xraygm2= xraygasmass2/init_gasmass2
windgm2= windmass2/init_gasmass2



; f
;
frun= ['/raid4/tcox/bs/b0f', '/raid4/tcox/bs/b1f', '/raid4/tcox/bs/b2f', '/raid4/tcox/bs/b3f', $
                '/raid4/tcox/bs/b4f', '/raid4/tcox/bs/b5f', '/raid4/tcox/bs/b6f']

fload_windinfo, frun, windmass3, windz3
fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3
fload_hotgasinfo, frun, init_gasmass3, rem_gasmass3
;fload_fruns_sigma_e, frun, xvar3

gasconsumption3= rem_gasmass3/init_gasmass3

xraygm3= xraygasmass3/init_gasmass3
windgm3= windmass3/init_gasmass3



; k
;
frun= ['/raid4/tcox/bs/b0k', '/raid4/tcox/bs/b1k', '/raid4/tcox/bs/b2k', '/raid4/tcox/bs/b3k', $
                '/raid4/tcox/bs/b4k', '/raid4/tcox/bs/b5k', '/raid4/tcox/bs/b6k']

fload_windinfo, frun, windmass4, windz4
fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4
fload_hotgasinfo, frun, init_gasmass4, rem_gasmass4
;fload_fruns_sigma_e, frun, xvar4

gasconsumption4= rem_gasmass4/init_gasmass4

xraygm4= xraygasmass4/init_gasmass4
windgm4= windmass4/init_gasmass4





; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(xvar1)
   xraygm_100= fltarr(n_f)
   xraygm_err_100= fltarr(n_f)
   windgm_100= fltarr(n_f)
   windgm_err_100= fltarr(n_f)
   xraygasz_100= fltarr(n_f)
   xraygasz_err_100= fltarr(n_f)
   windz_100= fltarr(n_f)
   windz_err_100= fltarr(n_f)
   lx_max_100= fltarr(n_f)
   lx_max_err_100= fltarr(n_f)
   lx_rem_100= fltarr(n_f)
   lx_rem_err_100= fltarr(n_f)

   for i=0,n_f-1 do begin
	vars[0]= xraygm1[i]
	vars[1]= xraygm2[i]
	vars[2]= xraygm3[i]
	vars[3]= xraygm4[i]
	xraygm_100[i]= mean(vars)
	xraygm_err_100[i]= sqrt(variance(vars))

	vars[0]= windgm1[i]
	vars[1]= windgm2[i]
	vars[2]= windgm3[i]
	vars[3]= windgm4[i]
	windgm_100[i]= mean(vars)
	windgm_err_100[i]= sqrt(variance(vars))

	vars[0]= xraygasz1[i]
	vars[1]= xraygasz2[i]
	vars[2]= xraygasz3[i]
	vars[3]= xraygasz4[i]
	xraygasz_100[i]= mean(vars)
	xraygasz_err_100[i]= sqrt(variance(vars))

	vars[0]= windz1[i]
	vars[1]= windz2[i]
	vars[2]= windz3[i]
	vars[3]= windz4[i]
	windz_100[i]= mean(vars)
	windz_err_100[i]= sqrt(variance(vars))

	vars[0]= lx_max1[i]
	vars[1]= lx_max2[i]
	vars[2]= lx_max3[i]
	vars[3]= lx_max4[i]
	lx_max_100[i]= mean(vars)
	lx_max_err_100[i]= sqrt(variance(vars))

	vars[0]= lx_rem1[i]
	vars[1]= lx_rem2[i]
	vars[2]= lx_rem3[i]
	vars[3]= lx_rem4[i]
	lx_rem_100[i]= mean(vars)
	lx_rem_err_100[i]= sqrt(variance(vars))

	;vars[0]= xvar1[i]
	;vars[1]= xvar2[i]
	;sigma_100[i]= mean(vars)
	;sigma_err_100[i]= sqrt(variance(vars))

   endfor
endif






;------------------------------------------------------------------------------------

xvar1= [50,80.0,115.0,160.0,225.0,320.0,500.0]


; Print thie mess up
; -------------------
filename='hotgasmasses.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=22, newysize=18

x0= 0.10 & xsize= 0.225
x1= x0+xsize
x2= x0+xsize+xsize
x3= x0+xsize+xsize+xsize
x4= x0+xsize+xsize+xsize+xsize


y0= 0.09 & ysize=0.3
y1= y0+ysize
y2= y0+ysize+ysize
y3= y0+ysize+ysize+ysize




; --------------------------------------------------------------------------------------



;  X-ray Luminosities
; ---------------------------
yaxistitle='!6L!DX!N (ers s!E-1!N)'
ymax = 45.0
ymin= 36.0


overplot_panel, xvar1, lx_max_5, lx_max_err_5, lx_rem_5, lx_rem_err_5, $
                                x0, y2, x1, y3, $
                                xmax, xmin, ymax, ymin, $
				 yaxistitle=yaxistitle, /plot_xray_lum, /plabel

overplot_panel, xvar1, lx_max_20, lx_max_err_20, lx_rem_20, lx_rem_err_20, $
                                x1, y2, x2, y3, $ 
                                xmax, xmin, ymax, ymin, /plot_xray_lum

overplot_panel, xvar1, lx_max_40, lx_max_err_40, lx_rem_40, lx_rem_err_40, $
                                x2, y2, x3, y3, $ 
                                xmax, xmin, ymax, ymin, /plot_xray_lum

overplot_panel, xvar1, lx_max_100, lx_max_err_100, lx_rem_100, lx_rem_err_100, $
                                x3, y2, x4, y3, $ 
                                xmax, xmin, ymax, ymin, /plot_xray_lum




; --------------------------------------------------------------------------------------



;  Mass fractions
; ---------------------------
yaxistitle='!6mass fraction'
ymax = 0.44
;ymax = 0.25
ymin= 0.0


overplot_panel, xvar1, xraygm_5, xraygm_err_5, windgm_5, windgm_err_5, $
				x0, y1, x1, y2, $
				xmax, xmin, ymax, ymin, $
				yaxistitle=yaxistitle, /showtotal

overplot_panel, xvar1, xraygm_20, xraygm_err_20, windgm_20, windgm_err_20, $
				x1, y1, x2, y2, $
				xmax, xmin, ymax, ymin, /showtotal

overplot_panel, xvar1, xraygm_40, xraygm_err_40, windgm_40, windgm_err_40, $
				x2, y1, x3, y2, $
				xmax, xmin, ymax, ymin, /showtotal, /plabel

overplot_panel, xvar1, xraygm_100, xraygm_err_100, windgm_100, windgm_err_100, $
				x3, y1, x4, y2, $
				xmax, xmin, ymax, ymin, /showtotal






; --------------------------------------------------------------------------------------



;  Metallicity
; -------------------
yaxistitle='!6Metallicity (Z!D!9n!6!N)'
ymax = 3.0
ymin= 0.05


overplot_panel, xvar1, xraygasz_5, xraygasz_err_5, windz_5, windz_err_5, $
				x0, y0, x1, y1, $
				xmax, xmin, ymax, ymin, /ylog, $
				yaxistitle=yaxistitle, $
				xaxistitle=xaxistitle

overplot_panel, xvar1, xraygasz_20, xraygasz_err_20, windz_20, windz_err_20, $
				x1, y0, x2, y1, $
				xmax, xmin, ymax, ymin, /ylog, $
				xaxistitle=xaxistitle

overplot_panel, xvar1, xraygasz_40, xraygasz_err_40, windz_40, windz_err_40, $
				x2, y0, x3, y1, $
				xmax, xmin, ymax, ymin, /ylog, $
				xaxistitle=xaxistitle

overplot_panel, xvar1, xraygasz_100, xraygasz_err_100, windz_100, windz_err_100, $
				x3, y0, x4, y1, $
				xmax, xmin, ymax, ymin, /ylog, $
				xaxistitle=xaxistitle






; --------------------------------------------------------------------------------------

xyouts, 0.12, 0.94, '5 %', /normal, charthick=8.0, size=1.6, color=0
xyouts, 0.345, 0.94, '20 %', /normal, charthick=8.0, size=1.6, color=0
xyouts, 0.57, 0.94, '40 %', /normal, charthick=8.0, size=1.6, color=0
;xyouts, 0.795, 0.94, '100 %', /normal, charthick=8.0, size=1.6, color=0
xyouts, 0.795, 0.94, '80 %', /normal, charthick=8.0, size=1.6, color=0

xyouts, 0.62, 0.64, 'X-ray', /normal, charthick=4.0, size=1.3, color=150
xyouts, 0.62, 0.61, 'Unbound', /normal, charthick=4.0, size=1.3, color=50

xyouts, 0.225, 0.770, '!6L!DX,max!N', /normal, charthick=4.0, size=1.5, color=0
xyouts, 0.225, 0.735, '!6L!DX,remnant!N', /normal, charthick=4.0, size=1.5, color=0

; sigma - log
;oplot, [55.0,70.0], [0.24,0.24], psym=-2, color=150, thick= 3.0, symsize= 1.5
;oplot, [55.0,70.0], [0.125,0.125], psym=-8, color=50, thick= 3.0, symsize= 1.5

; gas frac - linear
;oplot, [0.60,0.64], [0.53,0.53], psym=-2, color=150, thick= 3.0, symsize= 1.5
;oplot, [0.60,0.64], [0.49,0.49], psym=-8, color=50, thick= 3.0, symsize= 1.5




;---------
device, /close

end




; --------------------------------------------------------------------------------------



pro overplot_panel, xvar1, yvar1, yvar1_err, yvar2, yvar2_err, $
				x0, y0, x1, y1, $
				xmax, xmin, ymax, ymin, $
				yaxistitle=yaxistitle, $
				xaxistitle=xaxistitle, ylog=ylog, $
				showtotal=showtotal, $
				plot_xray_lum=plot_xray_lum, $
				plabel=plabel


f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill

!p.position= [x0, y0, x1, y1]

if keyword_set(xaxistitle) and keyword_set(yaxistitle) then begin
   if keyword_set(ylog) then begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
                xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, /ylog, $
                xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
                xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase
   endif else begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
	        xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, $
	        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase
   endelse
endif

if keyword_set(xaxistitle) and (not keyword_set(yaxistitle)) then begin
   if keyword_set(ylog) then begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
                xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, /ylog, $
                xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
                xtitle=xaxistitle, ytickformat='(a1)', /nodata, /noerase
   endif else begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
	        xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, $
	        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	        xtitle=xaxistitle, ytickformat='(a1)', /nodata, /noerase
   endelse
endif

if not keyword_set(xaxistitle) and keyword_set(yaxistitle) then begin
   if keyword_set(ylog) then begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
                xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, /ylog, $
                xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
                xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
   endif else begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
                xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, $
                xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
                xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
   endelse
endif

if not keyword_set(xaxistitle) and keyword_set(yaxistitle) then begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
	        xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, $
	        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
endif

if not keyword_set(xaxistitle) and (not keyword_set(yaxistitle)) then begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
                xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, $
                xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
                xtickformat='(a1)', ytickformat='(a1)', /nodata, /noerase
endif


if keyword_set(plot_xray_lum) then begin
	oplot, xvar1, yvar1, psym=-3, color=0, thick= 7.0, symsize= 1.5
	oploterror, xvar1, yvar1, yvar1_err, psym=-3, errcolor=0, color=0, thick= 3.0, errthick= 3.0

	oplot, xvar1, yvar2, psym=-3, color=0, thick= 7.0, symsize= 1.5, linestyle= 2
	oploterror, xvar1, yvar2, yvar2_err, psym=-3, errcolor=0, color=0, thick= 3.0, errthick= 3.0, linestyle= 2

	; for labels
	if keyword_set(plabel) then begin
		oplot, [115.0,200.0], [38.5,38.5], psym=-3, color=0, thick= 7.0, symsize= 1.5
		oplot, [115.0,200.0], [37.5,37.5], psym=-3, color=0, thick= 7.0, symsize= 1.5, linestyle= 2
	endif
endif else begin
	oplot, xvar1, yvar1, psym=-2, color=150, thick= 3.0, symsize= 1.5
	oploterror, xvar1, yvar1, yvar1_err, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0

	f = (2.0*!pi/16.0)*findgen(17)
	usersym,cos(f),sin(f),/fill

	oplot, xvar1, yvar2, psym=-8, color=50, thick= 3.0, symsize= 1.5
	oploterror, xvar1, yvar2, yvar2_err, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0

	; for labels
	if keyword_set(plabel) then begin
		oplot, [55.0,80.0], [0.38,0.38], psym=-2, color=150, thick= 3.0, symsize= 1.5
		oplot, [55.0,80.0], [0.33,0.33], psym=-8, color=50, thick= 3.0, symsize= 1.5
	endif
endelse


; plot the sum as a dashed line
if keyword_set(showtotal) then begin
	oplot, xvar1, yvar1+yvar2, psym=-3, color=0, thick= 2.0, linestyle= 1
endif



end



;================================================================================================
