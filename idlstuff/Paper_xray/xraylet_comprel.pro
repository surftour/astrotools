


; ----------------------------
;  Read hotgas.txt file
; ----------------------------
pro read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, $
				gas_tot, gas_hot, gas_cold, gas_sf

hgasfile= frun+'/hotgas.txt'

spawn, "wc "+hgasfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then hgas_data= fltarr(9,lines)

openr, 1, hgasfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, hgas_data
close, 1


time= hgas_data[0,*]
temp_keV_X= hgas_data[1,*]
temp_keV= hgas_data[2,*]
temp_K= hgas_data[3,*]
entropy= hgas_data[4,*]
gas_tot= hgas_data[5,*]
gas_hot= hgas_data[6,*]
gas_cold= hgas_data[7,*]
gas_sf= hgas_data[8,*]


end




; ----------------------------
;  Read xrays.txt file
; ----------------------------
pro read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, $
				xray, xray_sf, $
				xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

hgasfile= frun+'/xrays.txt'

spawn, "wc "+hgasfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then hgas_data= fltarr(11,lines)

openr, 1, hgasfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, hgas_data
close, 1


time= hgas_data[0,*]
temp_keV_X= hgas_data[1,*]

; with new xrays.txt file (i.e. it has X-ray emission weighted Z)
z_X= hgas_data[2,*]
mass_xraygas= hgas_data[3,*]
entropy= hgas_data[4,*]
xray= hgas_data[5,*]
xray_sf= hgas_data[6,*]
xray_rs_s= hgas_data[7,*]
xray_rs_h= hgas_data[8,*]
xray_rs0_s= hgas_data[9,*]
xray_rs0_h= hgas_data[10,*]

end





; -------------------------------------------------------------------------------------------





;   xray info
; ---------------
pro fload_xrayinfo, frun, xraygasmass, xraygasz, lx_max, lx_rem, tx


   ns= n_elements(frun)
   allxmass= fltarr(ns)
   allxz= fltarr(ns)
   allxmax= fltarr(ns)
   allxrem= fltarr(ns)
   alltx= fltarr(ns)

   for i=0,ns-1 do begin

        read_xrays_file, frun[i], time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

        allxmass[i]= mass_xraygas(n_elements(time)-1)
        allxz[i]= z_X(n_elements(time)-1)
        allxmax[i]= max(xray_rs_s)
        allxrem[i]= xray_rs_s(n_elements(time)-1)
	alltx[i]= temp_keV_X(n_elements(time)-1)
   endfor

   xraygasmass= allxmass
   xraygasz= allxz
   lx_max= allxmax
   lx_rem= allxrem
   tx= alltx


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




;   blue luminosity
; --------------------------
pro fload_lbinfo, frun, lb

   ns= n_elements(frun)
   alllb= fltarr(ns)

   for i=0,ns-1 do begin

    ; L_b
    ; -----
    read_colors_file, frun[i], time,  mags
    lstidx= n_elements(time)-1
    Mb= mags[2,*]
    alllb[i]= -0.4*(Mb[lstidx]-5.51)
   endfor

   lb= alllb

end














; -------------------------------------------------------------------------------------------






;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------
;
;       L_x   -   T_x  
;
;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------

pro plot_Lx_Tx, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_Lx_Tx, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='LxTx.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
yaxistitle = "!6Log L!DX!N (ergs s!E-1!N)"
xaxistitle = "!6Temperature (keV)"

;xmax = 4.0
xmax = 3.0
;xmin = 0.10
xmin = 0.02
;xmin = 0.04

ymax = 43.5
;ymin = 38.6
ymin = 37.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
	/xlog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata





; -----------------------------------------------------------------------


print, '5 percent gas '

;
;  5 % Gas
;----------------

; e
;
frun= ['/raid4/tcox/zs/z0e', '/raid4/tcox/zs/z1e', '/raid4/tcox/zs/z2e', '/raid4/tcox/zs/z3e', $
                '/raid4/tcox/zs/z4e', '/raid4/tcox/zs/z5e', '/raid4/tcox/zs/z6e']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1



; h
;
frun= ['/raid4/tcox/zs/z0h', '/raid4/tcox/zs/z1h', '/raid4/tcox/zs/z2h', '/raid4/tcox/zs/z3h', $
                '/raid4/tcox/zs/z4h', '/raid4/tcox/zs/z5h', '/raid4/tcox/zs/z6h']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2


; f
;
frun= ['/raid4/tcox/zs/z0f', '/raid4/tcox/zs/z1f', '/raid4/tcox/zs/z2f', '/raid4/tcox/zs/z3f', $
                '/raid4/tcox/zs/z4f', '/raid4/tcox/zs/z5f', '/raid4/tcox/zs/z6f']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3


; k
;
frun= ['/raid4/tcox/zs/z0k', '/raid4/tcox/zs/z1k', '/raid4/tcox/zs/z2k', '/raid4/tcox/zs/z3k', $
                '/raid4/tcox/zs/z4k', '/raid4/tcox/zs/z5k', '/raid4/tcox/zs/z6k']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4




; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(frun)
   ;lx_max_5= fltarr(n_f)
   ;lx_max_err_5= fltarr(n_f)
   lx_rem_5= fltarr(n_f)
   lx_rem_err_5= fltarr(n_f)
   tx_5= fltarr(n_f)
   tx_err_5= fltarr(n_f)
   sigma_5= fltarr(n_f)
   sigma_err_5= fltarr(n_f)
   lb_5= fltarr(n_f)
   lb_err_5= fltarr(n_f)

   for i=0,n_f-1 do begin
        ;vars[0]= lx_max1[i]
        ;vars[1]= lx_max2[i]
        ;vars[2]= lx_max3[i]
        ;lx_max_5[i]= mean(vars)
        ;lx_max_err_5[i]= sqrt(variance(vars))

        vars[0]= lx_rem1[i]
        vars[1]= lx_rem2[i]
        vars[2]= lx_rem3[i]
        vars[3]= lx_rem4[i]
        lx_rem_5[i]= mean(vars)
        lx_rem_err_5[i]= sqrt(variance(vars))

        vars[0]= tx1[i]
        vars[1]= tx2[i]
        vars[2]= tx3[i]
        vars[3]= tx4[i]
        tx_5[i]= mean(vars)
        tx_err_5[i]= sqrt(variance(vars))

   endfor
endif


	usersym,[-1,-1,1,1,-1],[-1,1,1,-1,-1],thick=4.0,/fill

        oplot, tx_5, lx_rem_5, psym=-8, color=100, thick= 3.0, symsize= 1.5
        oploterror, tx_5, lx_rem_5, lx_rem_err_5, psym=-3, errcolor=100, color=100, thick= 3.0, errthick= 3.0
   



; -----------------------------------------------------------------------


print, '20 percent gas '

;
;  20 % Gas
;----------------
xvar1= [50.0, 80.0, 115.0, 160.0, 225.0, 320.0, 500.0]
xvar2= xvar1
xvar3= xvar1
xvar4= xvar1

; e
;
;frun= ['/raid4/tcox/vcs/vc0vc0e', '/raid4/tcox/vcs/vc1vc1e', '/raid4/tcox/vcs/vc2vc2e', '/raid4/tcox/vcs/vc3vc3_e', $
;                '/raid4/tcox/vcs/vc4vc4e', '/raid4/tcox/vcs/vc5vc5e', '/raid4/tcox/vcs/vc6vc6e']
frun= ['/raid4/tcox/es/e0e', '/raid4/tcox/es/e1e', '/raid4/tcox/es/e2e', '/raid4/tcox/es/e3e', $
                '/raid4/tcox/es/e4e', '/raid4/tcox/es/e5e', '/raid4/tcox/es/e6e']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
;fload_fruns_sigma_e, frun, xvar1
;fload_lbinfo, frun, lb1



; h
;
;frun= ['/raid4/tcox/vcs/vc0vc0h', '/raid4/tcox/vcs/vc1vc1', '/raid4/tcox/vcs/vc2vc2', '/raid4/tcox/vcs/vc3vc3', $
;                '/raid4/tcox/vcs/vc4vc4a', '/raid4/tcox/vcs/vc5vc5a', '/raid4/tcox/vcs/vc6vc6a']
frun= ['/raid4/tcox/es/e0h', '/raid4/tcox/es/e1h', '/raid4/tcox/es/e2h', '/raid4/tcox/es/e3h', $
                '/raid4/tcox/es/e4h', '/raid4/tcox/es/e5h', '/raid4/tcox/es/e6h']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
;fload_fruns_sigma_e, frun, xvar2
;fload_lbinfo, frun, lb2



; f
;
frun= ['/raid4/tcox/es/e0f', '/raid4/tcox/es/e1f', '/raid4/tcox/es/e2f', '/raid4/tcox/es/e3f', $
                '/raid4/tcox/es/e4f', '/raid4/tcox/es/e5f', '/raid4/tcox/es/e6f']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
;fload_lbinfo, frun, lb3



; k
;
frun= ['/raid4/tcox/es/e0k', '/raid4/tcox/es/e1k', '/raid4/tcox/es/e2k', '/raid4/tcox/es/e3k', $
                '/raid4/tcox/es/e4k', '/raid4/tcox/es/e5k', '/raid4/tcox/es/e6k']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
;fload_lbinfo, frun, lb4






; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(frun)
   ;lx_max_20= fltarr(n_f)
   ;lx_max_err_20= fltarr(n_f)
   lx_rem_20= fltarr(n_f)
   lx_rem_err_20= fltarr(n_f)
   tx_20= fltarr(n_f)
   tx_err_20= fltarr(n_f)
   sigma_20= fltarr(n_f)
   sigma_err_20= fltarr(n_f)
   lb_20= fltarr(n_f)
   lb_err_20= fltarr(n_f)

   for i=0,n_f-1 do begin
        ;vars[0]= lx_max1[i]
        ;vars[1]= lx_max2[i]
        ;vars[2]= lx_max3[i]
        ;vars[3]= lx_max4[i]
        ;lx_max_20[i]= mean(vars)
        ;lx_max_err_20[i]= sqrt(variance(vars))

        vars[0]= lx_rem1[i]
        vars[1]= lx_rem2[i]
        vars[2]= lx_rem3[i]
        vars[3]= lx_rem4[i]
        lx_rem_20[i]= mean(vars)
        lx_rem_err_20[i]= sqrt(variance(vars))

        vars[0]= tx1[i]
        vars[1]= tx2[i]
        vars[2]= tx3[i]
        vars[3]= tx4[i]
        tx_20[i]= mean(vars)
        tx_err_20[i]= sqrt(variance(vars))

        vars[0]= xvar1[i]
        vars[1]= xvar2[i]
        vars[2]= xvar3[i]
        vars[3]= xvar4[i]
        sigma_20[i]= mean(vars)
        sigma_err_20[i]= sqrt(variance(vars))

        ;vars[0]= lb1[i]
        ;vars[1]= lb2[i]
        ;vars[2]= lb3[i]
        ;lb_20[i]= mean(vars)
        ;lb_err_20[i]= sqrt(variance(vars))

   endfor
endif


        oplot, tx_20, lx_rem_20, psym=-2, color=50, thick= 3.0, symsize= 1.5
        oploterror, tx_20, lx_rem_20, lx_rem_err_20, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0
   


; -----------------------------------------------------------------------



print, '40 percent gas '

;
;  40 % Gas
;----------------

; e
;
frun= ['/raid4/tcox/ds/d0e2_q', '/raid4/tcox/ds/d1e2_q', '/raid4/tcox/ds/d2e2_q', '/raid4/tcox/ds/d3e7', $
                '/raid4/tcox/ds/d4e2_q', '/raid4/tcox/ds/d5e2_q', '/raid4/tcox/ds/d6e2_q']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
;fload_fruns_sigma_e, frun, xvar1
;fload_lbinfo, frun, lb1



; h
;
frun= ['/raid4/tcox/ds/d0h2_q', '/raid4/tcox/ds/d1h2_q', '/raid4/tcox/ds/d2h2_q', '/raid4/tcox/ds/d3h7', $
                '/raid4/tcox/ds/d4h2_q', '/raid4/tcox/ds/d5h2_q', '/raid4/tcox/ds/d6h2_q']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
;fload_fruns_sigma_e, frun, xvar2
;fload_lbinfo, frun, lb2



; k
;
frun= ['/raid4/tcox/ds/d0k2_q', '/raid4/tcox/ds/d1k2_q', '/raid4/tcox/ds/d2k2_q', '/raid4/tcox/ds/d3k7', $
                '/raid4/tcox/ds/d4k2_q', '/raid4/tcox/ds/d5k2_q', '/raid4/tcox/ds/d6k2_q']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
;fload_fruns_sigma_e, frun, xvar3
;fload_lbinfo, frun, lb3



; f
;
frun= ['/raid4/tcox/ds/d0f2_q', '/raid4/tcox/ds/d1f2_q', '/raid4/tcox/ds/d2f2_q', '/raid4/tcox/ds/d3k7', $
                '/raid4/tcox/ds/d4f2_q', '/raid4/tcox/ds/d5f2_q', '/raid4/tcox/ds/d6f2_q']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
;fload_fruns_sigma_e, frun, xvar3
;fload_lbinfo, frun, lb3





; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(frun)
   ;lx_max_40= fltarr(n_f)
   ;lx_max_err_40= fltarr(n_f)
   lx_rem_40= fltarr(n_f)
   lx_rem_err_40= fltarr(n_f)
   tx_40= fltarr(n_f)
   tx_err_40= fltarr(n_f)
   sigma_40= fltarr(n_f)
   sigma_err_40= fltarr(n_f)
   lb_40= fltarr(n_f)
   lb_err_40= fltarr(n_f)

   for i=0,n_f-1 do begin
        ;vars[0]= lx_max1[i]
        ;vars[1]= lx_max2[i]
        ;vars[2]= lx_max3[i]
        ;lx_max_40[i]= mean(vars)
        ;lx_max_err_40[i]= sqrt(variance(vars))

        vars[0]= lx_rem1[i]
        vars[1]= lx_rem2[i]
        vars[2]= lx_rem3[i]
        vars[3]= lx_rem4[i]
        lx_rem_40[i]= mean(vars)
        lx_rem_err_40[i]= sqrt(variance(vars))

        vars[0]= tx1[i]
        vars[1]= tx2[i]
        vars[2]= tx3[i]
        vars[3]= tx4[i]
        tx_40[i]= mean(vars)
        tx_err_40[i]= sqrt(variance(vars))

        ;vars[0]= xvar1[i]
        ;vars[1]= xvar2[i]
        ;vars[2]= xvar3[i]
        ;sigma_40[i]= mean(vars)
        ;sigma_err_40[i]= sqrt(variance(vars))

        ;vars[0]= lb1[i]
        ;vars[1]= lb2[i]
        ;vars[2]= lb3[i]
        ;lb_40[i]= mean(vars)
        ;lb_err_40[i]= sqrt(variance(vars))

   endfor
endif


        f = (2.0*!pi/16.0)*findgen(17)
        usersym,cos(f),sin(f),/fill

        oplot, tx_40, lx_rem_40, psym=-8, color=150, thick= 3.0, symsize= 1.5
        oploterror, tx_40, lx_rem_40, lx_rem_err_40, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0
   


; -----------------------------------------------------------------------


;print, '100 percent gas '
print, '80 percent gas '

;
;  100 % Gas
;----------------

; e
;
;frun= ['/raid4/tcox/As/A0e', '/raid4/tcox/As/A1e', '/raid4/tcox/As/A2e', '/raid4/tcox/As/A3e', $
;                '/raid4/tcox/As/A4e', '/raid4/tcox/As/A5e', '/raid4/tcox/As/A6e']
frun= ['/raid4/tcox/bs/b0e', '/raid4/tcox/bs/b1e', '/raid4/tcox/bs/b2e', '/raid4/tcox/bs/b3e', $
                '/raid4/tcox/bs/b4e', '/raid4/tcox/bs/b5e', '/raid4/tcox/bs/b6e']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
;tx1[2]= 0.15
;tx1[3]= 0.25



; h
;
;frun= ['/raid4/tcox/As/A0', '/raid4/tcox/As/A1', '/raid4/tcox/As/A2', '/raid4/tcox/As/A3', $
;                '/raid4/tcox/As/A4', '/raid4/tcox/As/A5', '/raid4/tcox/As/A6']
frun= ['/raid4/tcox/bs/b0h', '/raid4/tcox/bs/b1h', '/raid4/tcox/bs/b2h', '/raid4/tcox/bs/b3h', $
                '/raid4/tcox/bs/b4h', '/raid4/tcox/bs/b5h', '/raid4/tcox/bs/b6h']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
;tx2[2]= 0.25
;tx2[3]= 0.35



; f
;
frun= ['/raid4/tcox/bs/b0f', '/raid4/tcox/bs/b1f', '/raid4/tcox/bs/b2f', '/raid4/tcox/bs/b3f', $
                '/raid4/tcox/bs/b4f', '/raid4/tcox/bs/b5f', '/raid4/tcox/bs/b6f']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3



; k
;
frun= ['/raid4/tcox/bs/b0k', '/raid4/tcox/bs/b1k', '/raid4/tcox/bs/b2k', '/raid4/tcox/bs/b3k', $
                '/raid4/tcox/bs/b4k', '/raid4/tcox/bs/b5k', '/raid4/tcox/bs/b6k']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4





; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(frun)
   ;lx_max_100= fltarr(n_f)
   ;lx_max_err_100= fltarr(n_f)
   lx_rem_100= fltarr(n_f)
   lx_rem_err_100= fltarr(n_f)
   tx_100= fltarr(n_f)
   tx_err_100= fltarr(n_f)
   sigma_100= fltarr(n_f)
   sigma_err_100= fltarr(n_f)
   lb_100= fltarr(n_f)
   lb_err_100= fltarr(n_f)

   for i=0,n_f-1 do begin
        ;vars[0]= lx_max1[i]
        ;vars[1]= lx_max2[i]
        ;vars[2]= lx_max3[i]
        ;lx_max_100[i]= mean(vars)
        ;lx_max_err_100[i]= sqrt(variance(vars))

        vars[0]= lx_rem1[i]
        vars[1]= lx_rem2[i]
        vars[2]= lx_rem3[i]
        vars[3]= lx_rem4[i]
        lx_rem_100[i]= mean(vars)
        lx_rem_err_100[i]= sqrt(variance(vars))

        vars[0]= tx1[i]
        vars[1]= tx2[i]
        vars[2]= tx3[i]
        vars[3]= tx4[i]
        tx_100[i]= mean(vars)
        tx_err_100[i]= sqrt(variance(vars))

        ;vars[0]= lb1[i]
        ;vars[1]= lb2[i]
        ;vars[2]= lb3[i]
        ;lb_100[i]= mean(vars)
        ;lb_err_100[i]= sqrt(variance(vars))

   endfor
endif


        oplot, tx_100, lx_rem_100, psym=-6, color=200, thick= 3.0, symsize= 1.5
        oploterror, tx_100, lx_rem_100, lx_rem_err_100, psym=-3, errcolor=200, color=200, thick= 3.0, errthick= 3.0
   



; -----------------------------------------------------------------------





; extras
; ---------
;xyouts, 0.2, 0.92, 'bulge', /normal, charthick=3.0, size=1.33, color= 0


; Put actual data on there
; -------------------------
read_osullivan_03, loglx, loglb, tempx
;symsize= 0.2
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;oplot, tempx, loglx, psym=8, color= 0
oplot, tempx, loglx, psym=7, color= 0   ;, symsize=0.5


; a little key
;x0=0.06  &  x1= 0.28      ; first two for xmin/xmax= 0.04/4.0
;y0=42.6  &  y1= 43.0
x0=0.03415  &  x1= 0.16       ;  for xmin/xmax = 0.08/3.0
y0=42.4  &  y1= 42.8
;xyouts, 0.7, 0.28, 'OFP (2001)', color= 0, charthick=2.0, /normal
;oplot, [x0+0.13], [y0+0.75], psym=8, color= 0
xyouts, 0.31, 0.83, 'OPC (2003)', color= 0, charthick=2.0, /normal
oplot, [x0+0.007], [y0+0.19], psym=7, color= 0
oplot, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], psym=-3, color=0



; Lx_T fit
; ----------
x= [0.01,100.0]
y= 42.45 + 4.8*(alog10(x))

oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0

; draw arrow for merger time
;timemerge=1.05/0.7
;arrow, timemerge, 36.8, timemerge, 37.5, COLOR=0, THICK=3.0, hthick=3.0, /data




; Data key
; ----------
x0= 0.9  &  y0= 38.8   & dy= 0.40
xn= 0.8  &  yn= 0.36
usersym,[-1,-1,1,1,-1],[-1,1,1,-1,-1],thick=4.0,/fill
oplot, [x0], [y0], psym=-8, color=100, thick= 3.0, symsize= 1.5
xyouts, xn, yn, '5%', color=100, charthick=2.0, /normal, size= 1.2
;----
oplot, [x0], [y0-dy], psym=-2, color=50, thick= 3.0, symsize= 1.5
xyouts, xn, yn-0.05, '20%', color=50, charthick=2.0, /normal, size= 1.2
;----
f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill
oplot, [x0], [y0-2*dy], psym=-8, color=150, thick= 3.0, symsize= 1.5
xyouts, xn, yn-0.10, '40%', color=150, charthick=2.0, /normal, size= 1.2
;----
oplot, [x0], [y0-3*dy], psym=-6, color=200, thick= 3.0, symsize= 1.5
xyouts, xn, yn-0.15, '80%', color=200, charthick=2.0, /normal, size= 1.2

; enclosing box
oplot, [x0-0.15,x0-0.15,x0+1.05,x0+1.05,x0-0.15], $
	[y0-1.4,y0+0.2,y0+0.2,y0-1.4,y0-1.4], psym=-3, color=0




; helpful info
; --------------

; bh seed mass
;xyouts, 0.13, 42.6, 'BH seed mass', color= 0, charthick=1.2, /data
;arrow, 0.13, 42.5, 0.13, 42.0, COLOR=0, THICK=3.0, hthick=3.0, /data


;xyouts, 0.22, 0.22, 'Without discrete source contribution', /normal, color= 0, size=1.4


device, /close




end





;=====================================================================================









;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------
;
;       L_x   -   L_b  
;
;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------

pro plot_Lx_Lb, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_Lx_Lb, junk"
	print, " "
	print, " "
	return
endif

filename='LxLb.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
yaxistitle = "!6Log L!DX!N (ergs s!E-1!N)"
xaxistitle = "!6Log L!DB!N (L!D!9n!6!N)"
xmax = 12.5
xmin = 8.5
ymax = 43.5
ymin = 37.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
	;/xlog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; -------------------------------------------------------------------------


print, '5 percent gas '

;
;  5 % Gas
;----------------

; e
;
frun= ['/raid4/tcox/zs/z0e', '/raid4/tcox/zs/z1e', '/raid4/tcox/zs/z2e', '/raid4/tcox/zs/z3e', $
                '/raid4/tcox/zs/z4e', '/raid4/tcox/zs/z5e', '/raid4/tcox/zs/z6e']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
fload_lbinfo, frun, lb1


; h
;
frun= ['/raid4/tcox/zs/z0h', '/raid4/tcox/zs/z1h', '/raid4/tcox/zs/z2h', '/raid4/tcox/zs/z3h', $
                '/raid4/tcox/zs/z4h', '/raid4/tcox/zs/z5h', '/raid4/tcox/zs/z6h']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
fload_lbinfo, frun, lb2


; f
;
frun= ['/raid4/tcox/zs/z0f', '/raid4/tcox/zs/z1f', '/raid4/tcox/zs/z2f', '/raid4/tcox/zs/z3f', $
                '/raid4/tcox/zs/z4f', '/raid4/tcox/zs/z5f', '/raid4/tcox/zs/z6f']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
fload_lbinfo, frun, lb3


; k
;
frun= ['/raid4/tcox/zs/z0k', '/raid4/tcox/zs/z1k', '/raid4/tcox/zs/z2k', '/raid4/tcox/zs/z3k', $
                '/raid4/tcox/zs/z4k', '/raid4/tcox/zs/z5k', '/raid4/tcox/zs/z6k']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
fload_lbinfo, frun, lb4




; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(frun)
   ;lx_max_5= fltarr(n_f)
   ;lx_max_err_5= fltarr(n_f)
   lx_rem_5= fltarr(n_f)
   lx_rem_err_5= fltarr(n_f)
   tx_5= fltarr(n_f)
   tx_err_5= fltarr(n_f)
   sigma_5= fltarr(n_f)
   sigma_err_5= fltarr(n_f)
   lb_5= fltarr(n_f)
   lb_err_5= fltarr(n_f)

   for i=0,n_f-1 do begin
        ;vars[0]= lx_max1[i]
        ;vars[1]= lx_max2[i]
        ;vars[2]= lx_max3[i]
        ;lx_max_5[i]= mean(vars)
        ;lx_max_err_5[i]= sqrt(variance(vars))

        vars[0]= lx_rem1[i]
        vars[1]= lx_rem2[i]
        vars[2]= lx_rem3[i]
        vars[3]= lx_rem4[i]
        lx_rem_5[i]= mean(vars)
        lx_rem_err_5[i]= sqrt(variance(vars))

        vars[0]= tx1[i]
        vars[1]= tx2[i]
        vars[2]= tx3[i]
        vars[3]= tx4[i]
        tx_5[i]= mean(vars)
        tx_err_5[i]= sqrt(variance(vars))

        vars[0]= lb1[i]
        vars[1]= lb2[i]
        vars[2]= lb3[i]
        vars[3]= lb4[i]
        lb_5[i]= mean(vars)
        lb_err_5[i]= sqrt(variance(vars))

   endfor
endif


	usersym,[-1,-1,1,1,-1],[-1,1,1,-1,-1],thick=4.0,/fill

        oplot, lb_5, lx_rem_5, psym=-8, color=100, thick= 3.0, symsize= 1.5
        oploterror, lb_5, lx_rem_5, lx_rem_err_5, psym=-3, errcolor=100, color=100, thick= 3.0, errthick= 3.0
   



; -----------------------------------------------------------------------


print, '20 percent gas '

;
;  20 % Gas
;----------------

; e
;
;frun= ['/raid4/tcox/vcs/vc0vc0e', '/raid4/tcox/vcs/vc1vc1e', '/raid4/tcox/vcs/vc2vc2e', '/raid4/tcox/vcs/vc3vc3_e', $
;                '/raid4/tcox/vcs/vc4vc4e', '/raid4/tcox/vcs/vc5vc5e', '/raid4/tcox/vcs/vc6vc6e']
frun= ['/raid4/tcox/es/e0e', '/raid4/tcox/es/e1e', '/raid4/tcox/es/e2e', '/raid4/tcox/es/e3e', $
                '/raid4/tcox/es/e4e', '/raid4/tcox/es/e5e', '/raid4/tcox/es/e6e']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
fload_lbinfo, frun, lb1



; h
;
;frun= ['/raid4/tcox/vcs/vc0vc0h', '/raid4/tcox/vcs/vc1vc1', '/raid4/tcox/vcs/vc2vc2', '/raid4/tcox/vcs/vc3vc3', $
;                '/raid4/tcox/vcs/vc4vc4a', '/raid4/tcox/vcs/vc5vc5a', '/raid4/tcox/vcs/vc6vc6a']
frun= ['/raid4/tcox/es/e0h', '/raid4/tcox/es/e1h', '/raid4/tcox/es/e2h', '/raid4/tcox/es/e3h', $
                '/raid4/tcox/es/e4h', '/raid4/tcox/es/e5h', '/raid4/tcox/es/e6h']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
fload_lbinfo, frun, lb2



; f
;
frun= ['/raid4/tcox/es/e0f', '/raid4/tcox/es/e1f', '/raid4/tcox/es/e2f', '/raid4/tcox/es/e3f', $
                '/raid4/tcox/es/e4f', '/raid4/tcox/es/e5f', '/raid4/tcox/es/e6f']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
fload_lbinfo, frun, lb3



; k
;
frun= ['/raid4/tcox/es/e0k', '/raid4/tcox/es/e1k', '/raid4/tcox/es/e2k', '/raid4/tcox/es/e3k', $
                '/raid4/tcox/es/e4k', '/raid4/tcox/es/e5k', '/raid4/tcox/es/e6k']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
fload_lbinfo, frun, lb4



; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(frun)
   ;lx_max_20= fltarr(n_f)
   ;lx_max_err_20= fltarr(n_f)
   lx_rem_20= fltarr(n_f)
   lx_rem_err_20= fltarr(n_f)
   tx_20= fltarr(n_f)
   tx_err_20= fltarr(n_f)
   sigma_20= fltarr(n_f)
   sigma_err_20= fltarr(n_f)
   lb_20= fltarr(n_f)
   lb_err_20= fltarr(n_f)

   for i=0,n_f-1 do begin
        ;vars[0]= lx_max1[i]
        ;vars[1]= lx_max2[i]
        ;vars[2]= lx_max3[i]
        ;lx_max_20[i]= mean(vars)
        ;lx_max_err_20[i]= sqrt(variance(vars))

        vars[0]= lx_rem1[i]
        vars[1]= lx_rem2[i]
        vars[2]= lx_rem3[i]
        vars[3]= lx_rem4[i]
        lx_rem_20[i]= mean(vars)
        lx_rem_err_20[i]= sqrt(variance(vars))

        vars[0]= tx1[i]
        vars[1]= tx2[i]
        vars[2]= tx3[i]
        vars[3]= tx4[i]
        tx_20[i]= mean(vars)
        tx_err_20[i]= sqrt(variance(vars))

        ;vars[0]= xvar1[i]
        ;vars[1]= xvar2[i]
        ;vars[2]= xvar3[i]
        ;sigma_20[i]= mean(vars)
        ;sigma_err_20[i]= sqrt(variance(vars))

        vars[0]= lb1[i]
        vars[1]= lb2[i]
        vars[2]= lb3[i]
        vars[3]= lb4[i]
        lb_20[i]= mean(vars)
        lb_err_20[i]= sqrt(variance(vars))

   endfor
endif


        oplot, lb_20, lx_rem_20, psym=-2, color=50, thick= 3.0, symsize= 1.5
        oploterror, lb_20, lx_rem_20, lx_rem_err_20, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0
   


; -----------------------------------------------------------------------



print, '40 percent gas '

;
;  40 % Gas
;----------------

; e
;
frun= ['/raid4/tcox/ds/d0e2_q', '/raid4/tcox/ds/d1e2_q', '/raid4/tcox/ds/d2e2_q', '/raid4/tcox/ds/d3e7', $
                '/raid4/tcox/ds/d4e2_q', '/raid4/tcox/ds/d5e2_q', '/raid4/tcox/ds/d6e2_q']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
fload_lbinfo, frun, lb1



; h
;
frun= ['/raid4/tcox/ds/d0h2_q', '/raid4/tcox/ds/d1h2_q', '/raid4/tcox/ds/d2h2_q', '/raid4/tcox/ds/d3h7', $
                '/raid4/tcox/ds/d4h2_q', '/raid4/tcox/ds/d5h2_q', '/raid4/tcox/ds/d6h2_q']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
fload_lbinfo, frun, lb2



; k
;
frun= ['/raid4/tcox/ds/d0k2_q', '/raid4/tcox/ds/d1k2_q', '/raid4/tcox/ds/d2k2_q', '/raid4/tcox/ds/d3k7', $
                '/raid4/tcox/ds/d4k2_q', '/raid4/tcox/ds/d5k2_q', '/raid4/tcox/ds/d6k2_q']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
fload_lbinfo, frun, lb3



; f
;
frun= ['/raid4/tcox/ds/d0f2_q', '/raid4/tcox/ds/d1f2_q', '/raid4/tcox/ds/d2f2_q', '/raid4/tcox/ds/d3f7', $
                '/raid4/tcox/ds/d4f2_q', '/raid4/tcox/ds/d5f2_q', '/raid4/tcox/ds/d6f2_q']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
fload_lbinfo, frun, lb4





; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(frun)
   ;lx_max_40= fltarr(n_f)
   ;lx_max_err_40= fltarr(n_f)
   lx_rem_40= fltarr(n_f)
   lx_rem_err_40= fltarr(n_f)
   tx_40= fltarr(n_f)
   tx_err_40= fltarr(n_f)
   sigma_40= fltarr(n_f)
   sigma_err_40= fltarr(n_f)
   lb_40= fltarr(n_f)
   lb_err_40= fltarr(n_f)

   for i=0,n_f-1 do begin
        ;vars[0]= lx_max1[i]
        ;vars[1]= lx_max2[i]
        ;vars[2]= lx_max3[i]
        ;lx_max_40[i]= mean(vars)
        ;lx_max_err_40[i]= sqrt(variance(vars))

        vars[0]= lx_rem1[i]
        vars[1]= lx_rem2[i]
        vars[2]= lx_rem3[i]
        vars[3]= lx_rem4[i]
        lx_rem_40[i]= mean(vars)
        lx_rem_err_40[i]= sqrt(variance(vars))

        vars[0]= tx1[i]
        vars[1]= tx2[i]
        vars[2]= tx3[i]
        vars[3]= tx4[i]
        tx_40[i]= mean(vars)
        tx_err_40[i]= sqrt(variance(vars))

        vars[0]= lb1[i]
        vars[1]= lb2[i]
        vars[2]= lb3[i]
        vars[3]= lb4[i]
        lb_40[i]= mean(vars)
        lb_err_40[i]= sqrt(variance(vars))

   endfor
endif


        f = (2.0*!pi/16.0)*findgen(17)
        usersym,cos(f),sin(f),/fill

        oplot, lb_40, lx_rem_40, psym=-8, color=150, thick= 3.0, symsize= 1.5
        oploterror, lb_40, lx_rem_40, lx_rem_err_40, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0



; -----------------------------------------------------------------------


;print, '100 percent gas '
print, '80 percent gas '

;
;  100 % Gas
;----------------

; e
;
;frun= ['/raid4/tcox/As/A0e', '/raid4/tcox/As/A1e', '/raid4/tcox/As/A2e', '/raid4/tcox/As/A3e', $
;                '/raid4/tcox/As/A4e', '/raid4/tcox/As/A5e', '/raid4/tcox/As/A6e']
frun= ['/raid4/tcox/bs/b0e', '/raid4/tcox/bs/b1e', '/raid4/tcox/bs/b2e', '/raid4/tcox/bs/b3e', $
                '/raid4/tcox/bs/b4e', '/raid4/tcox/bs/b5e', '/raid4/tcox/bs/b6e']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
fload_lbinfo, frun, lb1



; h
;
;frun= ['/raid4/tcox/As/A0', '/raid4/tcox/As/A1', '/raid4/tcox/As/A2', '/raid4/tcox/As/A3', $
;                '/raid4/tcox/As/A4', '/raid4/tcox/As/A5', '/raid4/tcox/As/A6']
frun= ['/raid4/tcox/bs/b0h', '/raid4/tcox/bs/b1h', '/raid4/tcox/bs/b2h', '/raid4/tcox/bs/b3h', $
                '/raid4/tcox/bs/b4h', '/raid4/tcox/bs/b5h', '/raid4/tcox/bs/b6h']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
fload_lbinfo, frun, lb2



; f
;
frun= ['/raid4/tcox/bs/b0f', '/raid4/tcox/bs/b1f', '/raid4/tcox/bs/b2f', '/raid4/tcox/bs/b3f', $
                '/raid4/tcox/bs/b4f', '/raid4/tcox/bs/b5f', '/raid4/tcox/bs/b6f']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
fload_lbinfo, frun, lb3



; k
;
frun= ['/raid4/tcox/bs/b0k', '/raid4/tcox/bs/b1k', '/raid4/tcox/bs/b2k', '/raid4/tcox/bs/b3k', $
                '/raid4/tcox/bs/b4k', '/raid4/tcox/bs/b5k', '/raid4/tcox/bs/b6k']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
fload_lbinfo, frun, lb4







; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(frun)
   ;lx_max_100= fltarr(n_f)
   ;lx_max_err_100= fltarr(n_f)
   lx_rem_100= fltarr(n_f)
   lx_rem_err_100= fltarr(n_f)
   tx_100= fltarr(n_f)
   tx_err_100= fltarr(n_f)
   sigma_100= fltarr(n_f)
   sigma_err_100= fltarr(n_f)
   lb_100= fltarr(n_f)
   lb_err_100= fltarr(n_f)

   for i=0,n_f-1 do begin
        ;vars[0]= lx_max1[i]
        ;vars[1]= lx_max2[i]
        ;vars[2]= lx_max3[i]
        ;lx_max_100[i]= mean(vars)
        ;lx_max_err_100[i]= sqrt(variance(vars))

        vars[0]= lx_rem1[i]
        vars[1]= lx_rem2[i]
        vars[2]= lx_rem3[i]
        vars[3]= lx_rem4[i]
        lx_rem_100[i]= mean(vars)
        lx_rem_err_100[i]= sqrt(variance(vars))

        vars[0]= tx1[i]
        vars[1]= tx2[i]
        vars[2]= tx3[i]
        vars[3]= tx4[i]
        tx_100[i]= mean(vars)
        tx_err_100[i]= sqrt(variance(vars))

        vars[0]= lb1[i]
        vars[1]= lb2[i]
        vars[2]= lb3[i]
        vars[3]= lb4[i]
        lb_100[i]= mean(vars)
        lb_err_100[i]= sqrt(variance(vars))

   endfor
endif


        oplot, lb_100, lx_rem_100, psym=-6, color=200, thick= 3.0, symsize= 1.5
        oploterror, lb_100, lx_rem_100, lx_rem_err_100, psym=-3, errcolor=200, color=200, thick= 3.0, errthick= 3.0
   




; -------------------------------------------------------------------------




; extras
; ---------
;xyouts, 0.2, 0.92, 'bulge', /normal, charthick=3.0, size=1.33, color= 0


; Put actual data on there
; -------------------------
read_osullivan_01, loglb, loglx, ttype
symsize= 0.2
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
oplot, loglb, loglx, psym=8, color= 0
;oplot, loglb, loglx, psym=7, color= 0, symsize= 0.5

read_osullivan_03, loglx, loglb, tempx
oplot, loglb, loglx, psym=7, color= 0


; a little key
x0=10.9  &  x1= 12.15
y0=37.4  &  y1= 38.4
xyouts, 0.7, 0.28, 'OFP (2001)', color= 0, charthick=2.0, /normal
oplot, [x0+0.13], [y0+0.75], psym=8, color= 0
xyouts, 0.7, 0.22, 'OPC (2003)', color= 0, charthick=2.0, /normal
oplot, [x0+0.13], [y0+0.25], psym=7, color= 0
oplot, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], psym=-3, color=0





; Lx_Lb fit
; ----------
; fit for all E's (OFP '01)
;x= [1.0,10.0,100.0]
;y= 17.98 + 2.17*x
;oplot, x, y, psym=-3, linestyle= 3, thick=2.0, color= 0

; bright E galaxy slope (OPC '03)
x= [1.0,10.0,100.0]
;y= 11.9 + 2.7*x     ; this is the OPC fit, but B lums are for h=0.5
y= 12.5 + 2.7*x     ; corrected for this plot because we use B lums from OFP, h=0.75
oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0



; discrete source Lx  (Ciotti et al. 1991)
x= [1.0,10.0,100.0]
y= 29.45 + x
oplot, x, y, psym=-3, linestyle= 1, thick=2.0, color= 0


; helpful info
oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0




; discrete source Lx  (Ciotti et al. 1991)
; -----------------------
x= [1.0,10.0,100.0]
y= 29.45 + x
oplot, x, y, psym=-3, linestyle= 1, thick=2.0, color= 0




; Data key
; ----------
x0= 8.9  &  y0= 43.0   & dy= 0.4
xn= 0.30  &  yn= 0.88
usersym,[-1,-1,1,1,-1],[-1,1,1,-1,-1],thick=4.0,/fill
oplot, [x0], [y0], psym=-8, color=100, thick= 3.0, symsize= 1.5
xyouts, xn, yn, '5%', color=100, charthick=2.0, /normal, size= 1.2
;----
oplot, [x0], [y0-dy], psym=-2, color=50, thick= 3.0, symsize= 1.5
xyouts, xn, yn-0.05, '20%', color=50, charthick=2.0, /normal, size= 1.2
;----
f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill
oplot, [x0], [y0-2*dy], psym=-8, color=150, thick= 3.0, symsize= 1.5
xyouts, xn, yn-0.10, '40%', color=150, charthick=2.0, /normal, size= 1.2
;----
oplot, [x0], [y0-3*dy], psym=-6, color=200, thick= 3.0, symsize= 1.5
xyouts, xn, yn-0.15, '80%', color=200, charthick=2.0, /normal, size= 1.2

; enclosing box
oplot, [x0-0.15,x0-0.15,x0+0.65,x0+0.65,x0-0.15], $
	[y0-1.45,y0+0.2,y0+0.2,y0-1.45,y0-1.45], psym=-3, color=0




; helpful info
; --------------

; bh seed mass
;xyouts, 8.8, 42.6, 'BH seed mass', color= 0, charthick=1.2, /data
;arrow, 8.8, 42.5, 8.8, 42.0, COLOR=0, THICK=3.0, hthick=3.0, /data


;xyouts, 0.22, 0.22, 'Without discrete source contribution', /normal, color= 0, size=1.4

device, /close


end








;=====================================================================================







;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------
;
;       Sigma   -   T_x  
;
;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------

pro plot_sigma_Tx, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_sigma_Tx, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='SigmaTx.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
yaxistitle = "!7r!6 (km s!E-1!N)"
xaxistitle = "!6Temperature (keV)"
xmax = 3.0
xmin = 0.02
ymax = 450.0
;ymin = 38.6
ymin = 40.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
	/xlog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; ----------------------------------------------------------------------------------------



print, '5 percent gas '

;
;  5 % Gas
;----------------

; e
;
frun= ['/raid4/tcox/zs/z0e', '/raid4/tcox/zs/z1e', '/raid4/tcox/zs/z2e', '/raid4/tcox/zs/z3e', $
                '/raid4/tcox/zs/z4e', '/raid4/tcox/zs/z5e', '/raid4/tcox/zs/z6e']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
fload_fruns_sigma_e, frun, xvar1


; h
;
frun= ['/raid4/tcox/zs/z0h', '/raid4/tcox/zs/z1h', '/raid4/tcox/zs/z2h', '/raid4/tcox/zs/z3h', $
                '/raid4/tcox/zs/z4h', '/raid4/tcox/zs/z5h', '/raid4/tcox/zs/z6h']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
fload_fruns_sigma_e, frun, xvar2


; f
;
frun= ['/raid4/tcox/zs/z0f', '/raid4/tcox/zs/z1f', '/raid4/tcox/zs/z2f', '/raid4/tcox/zs/z3f', $
                '/raid4/tcox/zs/z4f', '/raid4/tcox/zs/z5f', '/raid4/tcox/zs/z6f']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
fload_fruns_sigma_e, frun, xvar3


; k
;
frun= ['/raid4/tcox/zs/z0k', '/raid4/tcox/zs/z1k', '/raid4/tcox/zs/z2k', '/raid4/tcox/zs/z3k', $
                '/raid4/tcox/zs/z4k', '/raid4/tcox/zs/z5k', '/raid4/tcox/zs/z6k']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
fload_fruns_sigma_e, frun, xvar4




; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(frun)
   ;lx_max_5= fltarr(n_f)
   ;lx_max_err_5= fltarr(n_f)
   lx_rem_5= fltarr(n_f)
   lx_rem_err_5= fltarr(n_f)
   tx_5= fltarr(n_f)
   tx_err_5= fltarr(n_f)
   sigma_5= fltarr(n_f)
   sigma_err_5= fltarr(n_f)
   lb_5= fltarr(n_f)
   lb_err_5= fltarr(n_f)

   for i=0,n_f-1 do begin
        ;vars[0]= lx_max1[i]
        ;vars[1]= lx_max2[i]
        ;vars[2]= lx_max3[i]
        ;lx_max_5[i]= mean(vars)
        ;lx_max_err_5[i]= sqrt(variance(vars))

        vars[0]= lx_rem1[i]
        vars[1]= lx_rem2[i]
        vars[2]= lx_rem3[i]
        vars[3]= lx_rem4[i]
        lx_rem_5[i]= mean(vars)
        lx_rem_err_5[i]= sqrt(variance(vars))

        vars[0]= tx1[i]
        vars[1]= tx2[i]
        vars[2]= tx3[i]
        vars[3]= tx4[i]
        tx_5[i]= mean(vars)
        tx_err_5[i]= sqrt(variance(vars))

        vars[0]= xvar1[i]
        vars[1]= xvar2[i]
        vars[2]= xvar3[i]
        vars[3]= xvar4[i]
        sigma_5[i]= mean(vars)
        sigma_err_5[i]= sqrt(variance(vars))

   endfor
endif


	usersym,[-1,-1,1,1,-1],[-1,1,1,-1,-1],thick=4.0,/fill

        oplot, tx_5, sigma_5, psym=-8, color=100, thick= 3.0, symsize= 1.5
        oploterror, tx_5, sigma_5, sigma_err_5, psym=-3, errcolor=100, color=100, thick= 3.0, errthick= 3.0
   



; -----------------------------------------------------------------------


print, '20 percent gas '

;
;  20 % Gas
;----------------
xvar1= [50.0, 80.0, 115.0, 160.0, 225.0, 320.0, 500.0]
xvar2= xvar1
xvar3= xvar1
xvar4= xvar1

; e
;
;frun= ['/raid4/tcox/vcs/vc0vc0e', '/raid4/tcox/vcs/vc1vc1e', '/raid4/tcox/vcs/vc2vc2e', '/raid4/tcox/vcs/vc3vc3_e', $
;                '/raid4/tcox/vcs/vc4vc4e', '/raid4/tcox/vcs/vc5vc5e', '/raid4/tcox/vcs/vc6vc6e']
frun= ['/raid4/tcox/es/e0e', '/raid4/tcox/es/e1e', '/raid4/tcox/es/e2e', '/raid4/tcox/es/e3e', $
                '/raid4/tcox/es/e4e', '/raid4/tcox/es/e5e', '/raid4/tcox/es/e6e']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
fload_fruns_sigma_e, frun, xvar1



; h
;
;frun= ['/raid4/tcox/vcs/vc0vc0h', '/raid4/tcox/vcs/vc1vc1', '/raid4/tcox/vcs/vc2vc2', '/raid4/tcox/vcs/vc3vc3', $
;                '/raid4/tcox/vcs/vc4vc4a', '/raid4/tcox/vcs/vc5vc5a', '/raid4/tcox/vcs/vc6vc6a']
frun= ['/raid4/tcox/es/e0h', '/raid4/tcox/es/e1h', '/raid4/tcox/es/e2h', '/raid4/tcox/es/e3h', $
                '/raid4/tcox/es/e4h', '/raid4/tcox/es/e5h', '/raid4/tcox/es/e6h']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
fload_fruns_sigma_e, frun, xvar2



; f
;
frun= ['/raid4/tcox/es/e0f', '/raid4/tcox/es/e1f', '/raid4/tcox/es/e2f', '/raid4/tcox/es/e3f', $
                '/raid4/tcox/es/e4f', '/raid4/tcox/es/e5f', '/raid4/tcox/es/e6f']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
fload_fruns_sigma_e, frun, xvar3



; k
;
frun= ['/raid4/tcox/es/e0k', '/raid4/tcox/es/e1k', '/raid4/tcox/es/e2k', '/raid4/tcox/es/e3k', $
                '/raid4/tcox/es/e4k', '/raid4/tcox/es/e5k', '/raid4/tcox/es/e6k']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
fload_fruns_sigma_e, frun, xvar4





; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(frun)
   ;lx_max_20= fltarr(n_f)
   ;lx_max_err_20= fltarr(n_f)
   lx_rem_20= fltarr(n_f)
   lx_rem_err_20= fltarr(n_f)
   tx_20= fltarr(n_f)
   tx_err_20= fltarr(n_f)
   sigma_20= fltarr(n_f)
   sigma_err_20= fltarr(n_f)
   lb_20= fltarr(n_f)
   lb_err_20= fltarr(n_f)

   for i=0,n_f-1 do begin
        ;vars[0]= lx_max1[i]
        ;vars[1]= lx_max2[i]
        ;vars[2]= lx_max3[i]
        ;vars[3]= lx_max4[i]
        ;lx_max_20[i]= mean(vars)
        ;lx_max_err_20[i]= sqrt(variance(vars))

        vars[0]= lx_rem1[i]
        vars[1]= lx_rem2[i]
        vars[2]= lx_rem3[i]
        vars[3]= lx_rem4[i]
        lx_rem_20[i]= mean(vars)
        lx_rem_err_20[i]= sqrt(variance(vars))

        vars[0]= tx1[i]
        vars[1]= tx2[i]
        vars[2]= tx3[i]
        vars[3]= tx4[i]
        tx_20[i]= mean(vars)
        tx_err_20[i]= sqrt(variance(vars))

        vars[0]= xvar1[i]
        vars[1]= xvar2[i]
        vars[2]= xvar3[i]
        vars[3]= xvar4[i]
        sigma_20[i]= mean(vars)
        sigma_err_20[i]= sqrt(variance(vars))

        ;vars[0]= lb1[i]
        ;vars[1]= lb2[i]
        ;vars[2]= lb3[i]
        ;lb_20[i]= mean(vars)
        ;lb_err_20[i]= sqrt(variance(vars))

   endfor
endif


        oplot, tx_20, sigma_20, psym=-2, color=50, thick= 3.0, symsize= 1.5
        oploterror, tx_20, sigma_20, sigma_err_20, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0
   


; -----------------------------------------------------------------------



print, '40 percent gas '

;
;  40 % Gas
;----------------

; e
;
frun= ['/raid4/tcox/ds/d0e2_q', '/raid4/tcox/ds/d1e2_q', '/raid4/tcox/ds/d2e2_q', '/raid4/tcox/ds/d3e7', $
                '/raid4/tcox/ds/d4e2_q', '/raid4/tcox/ds/d5e2_q', '/raid4/tcox/ds/d6e2_q']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
fload_fruns_sigma_e, frun, xvar1



; h
;
frun= ['/raid4/tcox/ds/d0h2_q', '/raid4/tcox/ds/d1h2_q', '/raid4/tcox/ds/d2h2_q', '/raid4/tcox/ds/d3h7', $
                '/raid4/tcox/ds/d4h2_q', '/raid4/tcox/ds/d5h2_q', '/raid4/tcox/ds/d6h2_q']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
fload_fruns_sigma_e, frun, xvar2



; k
;
frun= ['/raid4/tcox/ds/d0k2_q', '/raid4/tcox/ds/d1k2_q', '/raid4/tcox/ds/d2k2_q', '/raid4/tcox/ds/d3k7', $
                '/raid4/tcox/ds/d4k2_q', '/raid4/tcox/ds/d5k2_q', '/raid4/tcox/ds/d6k2_q']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
fload_fruns_sigma_e, frun, xvar3



; f
;
frun= ['/raid4/tcox/ds/d0f2_q', '/raid4/tcox/ds/d1f2_q', '/raid4/tcox/ds/d2f2_q', '/raid4/tcox/ds/d3f7', $
                '/raid4/tcox/ds/d4f2_q', '/raid4/tcox/ds/d5f2_q', '/raid4/tcox/ds/d6f2_q']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
fload_fruns_sigma_e, frun, xvar4





; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(xvar1)
   ;lx_max_40= fltarr(n_f)
   ;lx_max_err_40= fltarr(n_f)
   lx_rem_40= fltarr(n_f)
   lx_rem_err_40= fltarr(n_f)
   tx_40= fltarr(n_f)
   tx_err_40= fltarr(n_f)
   sigma_40= fltarr(n_f)
   sigma_err_40= fltarr(n_f)
   lb_40= fltarr(n_f)
   lb_err_40= fltarr(n_f)

   for i=0,n_f-1 do begin
        ;vars[0]= lx_max1[i]
        ;vars[1]= lx_max2[i]
        ;vars[2]= lx_max3[i]
        ;lx_max_40[i]= mean(vars)
        ;lx_max_err_40[i]= sqrt(variance(vars))

        vars[0]= lx_rem1[i]
        vars[1]= lx_rem2[i]
        vars[2]= lx_rem3[i]
        vars[3]= lx_rem4[i]
        lx_rem_40[i]= mean(vars)
        lx_rem_err_40[i]= sqrt(variance(vars))

        vars[0]= tx1[i]
        vars[1]= tx2[i]
        vars[2]= tx3[i]
        vars[3]= tx4[i]
        tx_40[i]= mean(vars)
        tx_err_40[i]= sqrt(variance(vars))

        vars[0]= xvar1[i]
        vars[1]= xvar2[i]
        vars[2]= xvar3[i]
        vars[3]= xvar4[i]
        sigma_40[i]= mean(vars)
        sigma_err_40[i]= sqrt(variance(vars))

        ;vars[0]= lb1[i]
        ;vars[1]= lb2[i]
        ;vars[2]= lb3[i]
        ;lb_40[i]= mean(vars)
        ;lb_err_40[i]= sqrt(variance(vars))

   endfor
endif


        f = (2.0*!pi/16.0)*findgen(17)
        usersym,cos(f),sin(f),/fill

        oplot, tx_40, sigma_40, psym=-8, color=150, thick= 3.0, symsize= 1.5
        oploterror, tx_40, sigma_40, sigma_err_40, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0



; -------------------------------------------------------------------------



;print, '100 percent gas '
print, '80 percent gas '

;
;  100 % Gas
;----------------

; e
;
;frun= ['/raid4/tcox/As/A0e', '/raid4/tcox/As/A1e', '/raid4/tcox/As/A2e', '/raid4/tcox/As/A3e', $
;                '/raid4/tcox/As/A4e', '/raid4/tcox/As/A5e', '/raid4/tcox/As/A6e']
frun= ['/raid4/tcox/bs/b0e', '/raid4/tcox/bs/b1e', '/raid4/tcox/bs/b2e', '/raid4/tcox/bs/b3e', $
                '/raid4/tcox/bs/b4e', '/raid4/tcox/bs/b5e', '/raid4/tcox/bs/b6e']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
;fload_fruns_sigma_e, frun, xvar1
;tx1[2]= 0.15
;tx1[3]= 0.25



; h
;
;frun= ['/raid4/tcox/As/A0', '/raid4/tcox/As/A1', '/raid4/tcox/As/A2', '/raid4/tcox/As/A3', $
;                '/raid4/tcox/As/A4', '/raid4/tcox/As/A5', '/raid4/tcox/As/A6']
frun= ['/raid4/tcox/bs/b0h', '/raid4/tcox/bs/b1h', '/raid4/tcox/bs/b2h', '/raid4/tcox/bs/b3h', $
                '/raid4/tcox/bs/b4h', '/raid4/tcox/bs/b5h', '/raid4/tcox/bs/b6h']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
;fload_fruns_sigma_e, frun, xvar2
;tx1[2]= 0.25
;tx1[3]= 0.35



; f
;
frun= ['/raid4/tcox/bs/b0f', '/raid4/tcox/bs/b1f', '/raid4/tcox/bs/b2f', '/raid4/tcox/bs/b3f', $
                '/raid4/tcox/bs/b4f', '/raid4/tcox/bs/b5f', '/raid4/tcox/bs/b6f']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
;fload_fruns_sigma_e, frun, xvar3



; k
;
frun= ['/raid4/tcox/bs/b0k', '/raid4/tcox/bs/b1k', '/raid4/tcox/bs/b2k', '/raid4/tcox/bs/b3k', $
                '/raid4/tcox/bs/b4k', '/raid4/tcox/bs/b5k', '/raid4/tcox/bs/b6k']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
;fload_fruns_sigma_e, frun, xvar4




; average the runs
; ------------------
;average_three_runs = 0
average_three_runs = 1
if average_three_runs eq 1 then begin
   vars= fltarr(4)
   n_f= n_elements(frun)
   ;lx_max_100= fltarr(n_f)
   ;lx_max_err_100= fltarr(n_f)
   lx_rem_100= fltarr(n_f)
   lx_rem_err_100= fltarr(n_f)
   tx_100= fltarr(n_f)
   tx_err_100= fltarr(n_f)
   sigma_100= fltarr(n_f)
   sigma_err_100= fltarr(n_f)
   lb_100= fltarr(n_f)
   lb_err_100= fltarr(n_f)

   for i=0,n_f-1 do begin
        ;vars[0]= lx_max1[i]
        ;vars[1]= lx_max2[i]
        ;vars[2]= lx_max3[i]
        ;vars[3]= lx_max4[i]
        ;lx_max_100[i]= mean(vars)
        ;lx_max_err_100[i]= sqrt(variance(vars))

        vars[0]= lx_rem1[i]
        vars[1]= lx_rem2[i]
        vars[2]= lx_rem3[i]
        vars[3]= lx_rem4[i]
        lx_rem_100[i]= mean(vars)
        lx_rem_err_100[i]= sqrt(variance(vars))

        vars[0]= tx1[i]
        vars[1]= tx2[i]
        vars[2]= tx3[i]
        vars[3]= tx4[i]
        tx_100[i]= mean(vars)
        tx_err_100[i]= sqrt(variance(vars))

        vars[0]= xvar1[i]
        vars[1]= xvar2[i]
        vars[2]= xvar3[i]
        vars[3]= xvar4[i]
        sigma_100[i]= mean(vars)
        sigma_err_100[i]= sqrt(variance(vars))

   endfor
endif


        oplot, tx_100, sigma_100, psym=-6, color=200, thick= 3.0, symsize= 1.5
        oploterror, tx_100, sigma_100, sigma_err_100, psym=-3, errcolor=200, color=200, thick= 3.0, errthick= 3.0
   




; -------------------------------------------------------------------------




; extras
; ---------
;xyouts, 0.2, 0.92, 'bulge', /normal, charthick=3.0, size=1.33, color= 0


; Put actual data on there
; -------------------------
read_osullivan_03, loglx, loglb, tempx
read_osullivan_03b, sigma
;symsize= 0.2
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;oplot, tempx, loglx, psym=8, color= 0
oplot, tempx, sigma, psym=7, color= 0   ;, symsize=0.5


; a little key
x0=0.031  &  x1= 0.15
y0=295.0  &  y1= 345.0
;xyouts, 0.7, 0.28, 'OFP (2001)', color= 0, charthick=2.0, /normal
;oplot, [x0+0.13], [y0+0.75], psym=8, color= 0
xyouts, 0.31, 0.83, 'OPC (2003)', color= 0, charthick=2.0, /normal
oplot, [x0+0.007], [y0+25.0], psym=7, color= 0
oplot, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], psym=-3, color=0



; Sigma_T Beta_spec=1 (I think?)
; -------------------
x= [0.001,10.0]
y= 10^(0.480426*alog10(x) + alog10(300.0))

oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0



; Data key
; ----------
x0= 0.9  &  y0= 85.0   & dy= 10.0
xn= 0.8  &  yn= 0.39
usersym,[-1,-1,1,1,-1],[-1,1,1,-1,-1],thick=4.0,/fill
oplot, [x0], [y0], psym=-8, color=100, thick= 3.0, symsize= 1.5
xyouts, xn, yn, '5%', color=100, charthick=2.0, /normal, size= 1.2
;----
oplot, [x0], [y0-dy], psym=-2, color=50, thick= 3.0, symsize= 1.5
xyouts, xn, yn-0.05, '20%', color=50, charthick=2.0, /normal, size= 1.2
;----
f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill
oplot, [x0], [y0-2*dy], psym=-8, color=150, thick= 3.0, symsize= 1.5
xyouts, xn, yn-0.10, '40%', color=150, charthick=2.0, /normal, size= 1.2
;----
oplot, [x0], [y0-3*dy], psym=-6, color=200, thick= 3.0, symsize= 1.5
xyouts, xn, yn-0.15, '80%', color=200, charthick=2.0, /normal, size= 1.2

; enclosing box
oplot, [x0-0.15,x0-0.15,x0+1.05,x0+1.05,x0-0.15], $
	[y0-35,y0+6,y0+6,y0-35,y0-35], psym=-3, color=0




; helpful info
; --------------

; bh seed mass
;xyouts, 0.13, 42.6, 'BH seed mass', color= 0, charthick=1.2, /data
;arrow, 0.13, 42.5, 0.13, 42.0, COLOR=0, THICK=3.0, hthick=3.0, /data


;xyouts, 0.22, 0.22, 'Without discrete source contribution', /normal, color= 0, size=1.4


device, /close




end









;====================================================================================











; genereric plotting program
; ----------------------------
pro readandplot_lx_and_else, frun, pointselection, msg, y0, $
				tx=tx, $
				lb=lb, $
				sig=sig, $
				yevolfac=yevolfac

    thispsym= -3
    thiscolor= 150   ;red 

    ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
    read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

    ; get last index
    lstidx= n_elements(time)-1  ; take last

    ; default y-axis
    ;yval= xray[lstidx]
    yval= xray_rs_s[lstidx]

    ; T_x
    ; ----
    if keyword_set(tx) or keyword_set(sig) then begin
	xval= temp_keV_X[lstidx]
    endif

    ; L_b
    ; -----
    read_colors_file, frun, time,  mags
    lstidx= n_elements(time)-1
    Mb= mags[2,*]
    Lumb= -0.4*(Mb-5.51)

    if keyword_set(lb) then begin
	xval= Lumb[lstidx]
    endif

    ; Sigma
    ; ------
    if keyword_set(sig) then begin
	read_sigma_file, frun, time, Asigxy, Asigxz, Asigyz, Asigavg, Asigerr, $
                                Bsigxy, Bsigxz, Bsigyz, Bsigavg, Bsigerr, $
                                gsigavg, gsigerr
	slstidx= n_elements(time)-1
	yval= Asigavg[slstidx]
    endif

    ; add point source component to L_x
    ;L_x_discrete = 29.5 + Lumb[lstidx]
    ;yval= yval - 38.0
    ;L_x_discrete= L_x_discrete - 38.0
    ;yval= (10^(yval)) + (10^(L_x_discrete))
    ;yval= alog10(yval) + 38.0


;  open (or filled) square
; --------------------------
if pointselection eq 1 then begin
        symsize= 1.5
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=6.0
        symsel= 8
        symcolor= 50
endif

;  x's
; -----
if pointselection eq 2 then begin
        symsel= 7
        symcolor= 100
endif

;  open circle
; --------------
if pointselection eq 3 then begin
        symsize= 1.5
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        ;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
        symsel= 8
        symcolor= 150
endif

;  open triangle
; ----------------
if pointselection eq 4 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
        usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
        symsel= 8
        symcolor= 0
endif

;  asterisk
; ----------
if pointselection eq 5 then begin
        symsel= 2
        symcolor= 200
endif

;  filled square
; --------------------------
if pointselection eq 6 then begin
        symsize= 1.6
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 8
        symcolor= 100
endif

;  diamond
; ----------
if pointselection eq 7 then begin
        symsel= 4 
        symcolor= 120
endif

;  triangle (open)
; -----------------
if pointselection eq 8 then begin
        symsel= 5 
        symcolor= 180
endif

;  plus sign
; ----------
if pointselection eq 9 then begin
        symsel= 1 
        symcolor= 170
endif


print, xval, yval

    oplot, [xval], [yval], thick=3.0, psym=symsel, color= symcolor

    ;oplot, [2.2, 2.4], [41.48,41.48], thick=3.0, psym=-2, color= 50
    if not keyword_set(msg) then msg= ' '
    if not keyword_set(y0) then y0= 0.0
    xyouts, 0.73, y0, msg, /normal, charthick=3.0, size=1.33, color=thiscolor

if not keyword_set(sig) and yevolfac gt 0 then begin
   ; draw arrow for pt. source contribution
   L_x_discrete = 29.5 + Lumb[lstidx]
   yval_wpts= yval - 38.0
   L_x_discrete= L_x_discrete - 38.0
   yval_wpts= (10^(yval_wpts)) + (10^(L_x_discrete))
   yval_wpts= alog10(yval_wpts) + 38.0
   arrow, [xval], [yval], [xval], [yval_wpts], COLOR=100, THICK=3.0, hthick=3.0, /data

   yval=yval_wpts

   ; draw arrow for 5 Gyr evolution
   xevolfac=0.0
   if not keyword_set(yevolfac) then yevolfac= 0.0
   yevolfac=alog10(1+yevolfac)
   if keyword_set(lb) then xevolfac=-1.0*alog10(2.0)
   arrow, [xval], [yval], [xval+xevolfac], [yval+yevolfac], COLOR=symcolor, THICK=3.0, hthick=3.0, /data

endif


end




; genereric plotting program
;  - only this one does time
;    evolution
; ----------------------------
pro readandplot_lx_and_else_time, frun, msg, y0, $
				tx=tx, $
				lb=lb

    thispsym= -3
    thiscolor= 150   ;red 

    ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
    read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h


    ;yval= xray
    yval= xray_rs_s

    ; T_x
    ; ----
    if keyword_set(tx) then begin
	xval= temp_keV_X
    endif

    ; L_b
    ; -----
    read_colors_file, frun, time,  mags
    lstidx= n_elements(time)-1
    Mb= mags[2,*]
    Lumb= -0.4*(Mb-5.51)

    if keyword_set(lb) then begin
	xval= Lumb
    endif

    ; add point source component to L_x
    ; ------------------------------------
    ;L_x_discrete = 29.5 + Lumb
    ;yval= yval - 38.0
    ;L_x_discrete= L_x_discrete - 38.0
    ;yval= (10^(yval)) + (10^(L_x_discrete))
    ;yval= alog10(yval) + 38.0


    oplot, xval, yval, thick=3.0, psym=thispsym, color= thiscolor

    ;oplot, [2.2, 2.4], [41.48,41.48], thick=3.0, psym=-2, color= 50
    if not keyword_set(msg) then msg= ' '
    if not keyword_set(y0) then y0= 0.0
    xyouts, 0.73, y0, msg, /normal, charthick=3.0, size=1.33, color=thiscolor

end



; -----------------------------------------------------------------------------------




;=========================================================
;
;  Use Raymond & Smith ('77, I think) to determine
;  the X-ray luminosity.
;
;
;=========================================================

pro load_raymondsmith_xraylum, hard_xray_lum, soft_xray_lum, $
					zero_metallicity=zero_metallicity


ngas= long(fload_npart(0))
;Ttime= float(fload_time(1))
Redshift= double(0.0)

GasMetallicity = (1/0.02)*fload_gas_metallicity(1)
GasNe = fload_gas_nume(1)
GasMass = 1.0e+10*fload_gas_mass(1)  ; in solar masses
GasTemp = float(fload_gas_temperature(1))   ; in K (for some reason this comes back in double)
GasHsml = fload_gas_hsml(1)


rho= fload_gas_rho(1)

; units
;hub= 0.7
hub= 1.0
ProtonMass=         1.6726e-24
UnitLength_in_cm=   3.085678d+21
UnitMass_in_g=      1.989d+43
UnitDensity_in_cgs= UnitMass_in_g/(UnitLength_in_cm^3)
print, UnitDensity_in_cgs

; hydrogen number density (nH cm-3)
GasHIRho= rho*0.76*(hub*hub)*UnitDensity_in_cgs/ProtonMass
GasHIRho= float(GasHIRho)

; electron number density (ne cm-3)
GasNeRho= float(GasNe*GasHIRho)

idx=where((rho lt 0.000854924) and (GasTemp ge 1.0e+5))
if idx(0) ne -1 then begin
    print, n_elements(idx), " out of ",ngas," particles will have diffuse X-ray emission."
    ngas= n_elements(idx)
    GasMetallicity= GasMetallicity(idx)
    GasNeRho= GasNeRho(idx)
    GasHIRho= GasHIRho(idx)
    GasMass= GasMass(idx)
    GasTemp= GasTemp(idx)
    GasHsml= GasHsml(idx)
endif else begin
    ngas= 0
endelse


; for testing purposes
;    ngas= 10L
;    GasMetallicity= GasMetallicity(5600:5609)
;    GasNeRho= GasNeRho(5600:5609)
;    GasHIRho= GasHIRho(5600:5609)
;    GasMass= GasMass(5600:5609)
;    GasTemp= GasTemp(5600:5609)
;    GasHsml= GasHsml(5600:5609)
;

if keyword_set(zero_metallicity) then begin
	GasMetallicity(*)= 0.0
endif


;
; new incarnation of code doesn't
; like the zero metallicities
;
idx=where(GasMetallicity le 0.0)
if idx(0) ne -1 then begin
    GasMetallicity(idx)= 1.0e-5
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;       Find luminosities
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if(ngas gt 0) then begin

        soft_xray_lum = fltarr(ngas)
        hard_xray_lum = fltarr(ngas)

        ;S = CALL_EXTERNAL('/home/brant/code/idl/raymond_smith/raymond_smith', $
        S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/RaymondSmith/raymond_smith', $
                'raymond_smith', $
                Ngas, $
                Redshift, $
                GasMetallicity, $
                GasNeRho, $
		GasHIRho, $
                GasMass, $
                GasTemp, $
                GasHsml, $
                soft_xray_lum, $
                hard_xray_lum)

        ;LUMINOSITIES are in h^-1 units

endif else begin
        print,'No gas, no x-ray luminosity.'
	hard_xray_lum= [0]
	soft_xray_lum= [0]
endelse



; brant returns this in solar luminosities
; multiply by 3.989e33 to get ergs/sec


hard_xray_lum = hard_xray_lum*3.989d33
soft_xray_lum = soft_xray_lum*3.989d33
print, "Total hard= ",total(hard_xray_lum)," erg sec-1"
print, "Total soft= ",total(soft_xray_lum)," erg sec-1"


end












;========================================
;
;   Read the O'Sullivan data
;
;========================================

; This is the O'Sullivan, Forbes & Ponman
; 2001, MNRAS, 328, 461 data:
; A catalogue and analysis of L_X of early-type galaxies 
;
pro read_osullivan_01, loglb, loglx, ttype


; this first file has spaces in the name which 
; screws things up, use  the second one.
;osullivanfile= '/home/tcox/osullivan_01_all.txt'

osullivanfile= '/home/tcox/OSullivan/osullivan_01_all2.txt'


;
;  Format of this file is:
;#
;#
;#
;#  Name		D	Log LB	 	Log LX	Source	T
;# 		(Mpc)	(LB)	 	(erg s1)	 	 
;ESO10114	30.12	9.93*	<	41.02	B	3.0
;ESO1074	38.89	10.22	<	40.94	B	4.0
; etc...
;

spawn, "wc "+osullivanfile,result
lines=long(result)
datalines=lines(0)-5
loglb= fltarr(datalines)
loglx= fltarr(datalines)
ttype= fltarr(datalines)

openr, 1, osullivanfile
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
	readf, 1, junk
	;print, junk
	tempjunk= strsplit(junk,/extract,count=count)
	name= tempjunk(0)
	;distance= float(tempjunk(1))
	loglb(i)= float(tempjunk(2))   ; doesn't need a trap for asterisk
	if count eq 7 then begin
		idontknow= tempjunk(3)
		loglx(i)= float(tempjunk(4))
		source= tempjunk(5)
		ttype(i)= float(tempjunk(6))
	endif else begin
		loglx(i)= float(tempjunk(3))
		source= tempjunk(4)
		ttype(i)= float(tempjunk(5))
	endelse

endfor

close, 1



end




;========================================
;========================================
; This is the O'Sullivan, Ponman & Collins
; 2003, MNRAS, 340, 1375 data:
; X-ray scaling properties of early-type galaxies
;
pro read_osullivan_03, loglx, loglb, tempx


; first we'll read the '01 data, to load the L_B
; values
osullivanfile01= '/home/tcox/OSullivan/osullivan_01_all2.txt'

spawn, "wc "+osullivanfile01,result
lines=long(result)
datalines01=lines(0)-5
name01= strarr(datalines01)
loglb01= fltarr(datalines01)
loglx01= fltarr(datalines01)
ttype01= fltarr(datalines01)

openr, 1, osullivanfile01
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
	readf, 1, junk
	tempjunk= strsplit(junk,/extract,count=count)
	name01(i)= tempjunk(0)
	loglb01(i)= float(tempjunk(2))   ; doesn't need a trap for asterisk
	if count eq 7 then begin
		loglx01(i)= float(tempjunk(4))
		ttype01(i)= float(tempjunk(6))
	endif else begin
		loglx01(i)= float(tempjunk(3))
		ttype01(i)= float(tempjunk(5))
	endelse

endfor

close, 1





; next we'll actually read the 03 data, and fix
; the values to send back
osullivanfile03= '/home/tcox/OSullivan/osullivan_03_lxfixed.txt'

spawn, "wc "+osullivanfile03,result
lines=long(result)
datalines03=lines(0)-5
loglb= fltarr(datalines03)
loglx= fltarr(datalines03)
tempx= fltarr(datalines03)

openr, 1, osullivanfile03
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
        readf, 1, junk
        tempjunk= strsplit(junk,/extract,count=count)
        name= tempjunk(0)
	idx= where(name eq name01)
	if idx(0) ne -1 then begin
		loglb(i)= loglb01(idx)
		;nh(i)= float(tempjunk(1))                 ; 10^21 cm^2
		tempx(i)= float(tempjunk(2))               ; keV
		;metallicity(i)= float(tempjunk(3))        ; in solar
		loglx(i)= float(tempjunk(4))
	endif else begin
		print, "PROBLEM: can't find ",name
	endelse

endfor

close, 1




end




;========================================
;========================================
; This is the O'Sullivan, Ponman & Collins
; 2003, MNRAS, 340, 1375 data (informational):
; X-ray scaling properties of early-type galaxies
;
pro read_osullivan_03b, sigma


; we read the 03 data, note that this is mainly
; informational data, so sigma, R_e, etc.  L_x, T_x,
; are all enclosed above.
osullivanfile03= '/home/tcox/OSullivan/osullivan_03_infofixed.txt'

spawn, "wc "+osullivanfile03,result
lines=long(result)
datalines03=lines(0)-5
sigma= fltarr(datalines03)

openr, 1, osullivanfile03
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
        readf, 1, junk
        tempjunk= strsplit(junk,/extract,count=count)
        name= tempjunk(0)
	sigma(i)= float(tempjunk(1))                 ; km sec^-1
	;distance(i)= float(tempjunk(3))             ; Mpc
	;re(i)= float(tempjunk(4))                   ; arcmin
	;ttype(i)= float(tempjunk(5))
	;environment(i)= float(tempjunk(6))

endfor

close, 1




end




;========================================
;========================================




