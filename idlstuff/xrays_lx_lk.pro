

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

        read_file_xrays, frun[i], time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

        allxmass[i]= mass_xraygas(n_elements(time)-1)
        allxz[i]= z_X(n_elements(time)-1)

	; raymond & smith
        ;allxmax[i]= max(xray_rs_s)
        ;allxrem[i]= xray_rs_s(n_elements(time)-1)

	; reg.
        allxmax[i]= max(xray)
        allxrem[i]= xray(n_elements(time)-1)
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

        read_file_sigma, frun[i], time, Asigxy, Asigxz, Asigyz, Asigavg, Asigerr, $
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
    read_file_colors, frun[i], time,  mags
    lstidx= n_elements(time)-1
    Mb= mags[2,*]
    alllb[i]= -0.4*(Mb[lstidx]-5.51)
   endfor

   lb= alllb

end




;   blue luminosity
; --------------------------
pro fload_lkinfo, frun, lk

   ns= n_elements(frun)
   alllk= fltarr(ns)

   for i=0,ns-1 do begin

    ; L_b
    ; -----
    read_file_colors, frun[i], time,  mags
    lstidx= n_elements(time)-1
    Mb= mags[8,*]
    alllk[i]= -0.4*(Mb[lstidx]-3.28)
   endfor

   lk= alllk

end





;=====================================================================================
;=====================================================================================
;=====================================================================================
;=====================================================================================
;=====================================================================================
;=====================================================================================









;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------
;
;       L_x   -   L_k
;
;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------

pro plot_Lx_Lk, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_Lx_Lk, junk"
	print, " "
	print, " "
	return
endif

filename='LxLk.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=18.0, newysize=12.0


;yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
yaxistitle = "!6Log L!DX!N (ergs s!E-1!N)"
xaxistitle = "!6Log L!DK!N (L!D!9n!6!N)"
xmax = 12.1
xmin = 10.5
ymax = 42.5
ymin = 38.0


;---------------------------

!p.position= [0.10, 0.14, 0.98, 0.98]

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
frun= ['data/zs/z0e', 'data/zs/z1e', 'data/zs/z2e', 'data/zs/z3e', $
                'data/zs/z4e', 'data/zs/z5e', 'data/zs/z6e']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
fload_lkinfo, frun, lb1


; h
;
frun= ['data/zs/z0h', 'data/zs/z1h', 'data/zs/z2h', 'data/zs/z3h', $
                'data/zs/z4h', 'data/zs/z5h', 'data/zs/z6h']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
fload_lkinfo, frun, lb2


; f
;
frun= ['data/zs/z0f', 'data/zs/z1f', 'data/zs/z2f', 'data/zs/z3f', $
                'data/zs/z4f', 'data/zs/z5f', 'data/zs/z6f']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
fload_lkinfo, frun, lb3


; k
;
frun= ['data/zs/z0k', 'data/zs/z1k', 'data/zs/z2k', 'data/zs/z3k', $
                'data/zs/z4k', 'data/zs/z5k', 'data/zs/z6k']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
fload_lkinfo, frun, lb4




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

        ;oplot, lb_5, lx_rem_5, psym=-8, color=100, thick= 3.0, symsize= 1.5
        ;oploterror, lb_5, lx_rem_5, lx_rem_err_5, psym=-3, errcolor=100, color=100, thick= 3.0, errthick= 3.0
   



; -----------------------------------------------------------------------


print, '20 percent gas '

;
;  20 % Gas
;----------------

; e
;
;frun= ['data/vcs/vc0vc0e', 'data/vcs/vc1vc1e', 'data/vcs/vc2vc2e', 'data/vcs/vc3vc3_e', $
;                'data/vcs/vc4vc4e', 'data/vcs/vc5vc5e', 'data/vcs/vc6vc6e']
;frun= ['data/es/e0e', 'data/es/e1e', 'data/es/e2e', 'data/es/e3e', $
;                'data/es/e4e', 'data/es/e5e', 'data/es/e6e']
frun= ['data/es/e2e', 'data/es/e3e', 'data/es/e4e', 'data/es/e5e']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
fload_lkinfo, frun, lb1



; h
;
;frun= ['data/vcs/vc0vc0h', 'data/vcs/vc1vc1', 'data/vcs/vc2vc2', 'data/vcs/vc3vc3', $
;                'data/vcs/vc4vc4a', 'data/vcs/vc5vc5a', 'data/vcs/vc6vc6a']
;frun= ['data/es/e0h', 'data/es/e1h', 'data/es/e2h', 'data/es/e3h', $
;                'data/es/e4h', 'data/es/e5h', 'data/es/e6h']
frun= ['data/es/e2h', 'data/es/e3h', 'data/es/e4h', 'data/es/e5h']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
fload_lkinfo, frun, lb2



; f
;
;frun= ['data/es/e0f', 'data/es/e1f', 'data/es/e2f', 'data/es/e3f', $
;                'data/es/e4f', 'data/es/e5f', 'data/es/e6f']
frun= ['data/es/e2f', 'data/es/e3f', 'data/es/e4f', 'data/es/e5f']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
fload_lkinfo, frun, lb3



; k
;
;frun= ['data/es/e0k', 'data/es/e1k', 'data/es/e2k', 'data/es/e3k', $
;                'data/es/e4k', 'data/es/e5k', 'data/es/e6k']
frun= ['data/es/e2k', 'data/es/e3k', 'data/es/e4k', 'data/es/e5k']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
fload_lkinfo, frun, lb4



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


        ;oplot, lb_20, lx_rem_20, psym=-2, color=50, thick= 3.0, symsize= 1.5
        ;oploterror, lb_20, lx_rem_20, lx_rem_err_20, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0
   
        oplot, [lb1, lb2, lb3, lb4],  [lx_rem1, lx_rem2, lx_rem3, lx_rem4], psym=2, color=150, thick= 5.0, symsize= 2.5


; -----------------------------------------------------------------------



print, '40 percent gas '

;
;  40 % Gas
;----------------

; e
;
;frun= ['data/ds/d0e2_q', 'data/ds/d1e2_q', 'data/ds/d2e2_q', 'data/ds/d3e7', $
;                'data/ds/d4e2_q', 'data/ds/d5e2_q', 'data/ds/d6e2_q']
frun= ['data/ds/d2e2_q', 'data/ds/d3e7', 'data/ds/d4e2_q', 'data/ds/d5e2_q']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
fload_lkinfo, frun, lb1



; h
;
;frun= ['data/ds/d0h2_q', 'data/ds/d1h2_q', 'data/ds/d2h2_q', 'data/ds/d3h7', $
;                'data/ds/d4h2_q', 'data/ds/d5h2_q', 'data/ds/d6h2_q']
frun= ['data/ds/d2h2_q', 'data/ds/d3h7', 'data/ds/d4h2_q', 'data/ds/d5h2_q']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
fload_lkinfo, frun, lb2



; k
;
;frun= ['data/ds/d0k2_q', 'data/ds/d1k2_q', 'data/ds/d2k2_q', 'data/ds/d3k7', $
;                'data/ds/d4k2_q', 'data/ds/d5k2_q', 'data/ds/d6k2_q']
frun= ['data/ds/d2k2_q', 'data/ds/d3k7', 'data/ds/d4k2_q', 'data/ds/d5k2_q']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
fload_lkinfo, frun, lb3



; f
;
;frun= ['data/ds/d0f2_q', 'data/ds/d1f2_q', 'data/ds/d2f2_q', 'data/ds/d3f7', $
;                'data/ds/d4f2_q', 'data/ds/d5f2_q', 'data/ds/d6f2_q']
frun= ['data/ds/d2f2_q', 'data/ds/d3f7', 'data/ds/d4f2_q', 'data/ds/d5f2_q']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
fload_lkinfo, frun, lb4





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

        ;oplot, lb_40, lx_rem_40, psym=-8, color=150, thick= 3.0, symsize= 1.5
        ;oploterror, lb_40, lx_rem_40, lx_rem_err_40, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0


        oplot, [lb1, lb2, lb3, lb4],  [lx_rem1, lx_rem2, lx_rem3, lx_rem4], psym=2, color=150, thick= 5.0, symsize= 2.5

; -----------------------------------------------------------------------


;print, '100 percent gas '
print, '80 percent gas '

;
;  100 % Gas
;----------------

; e
;
;frun= ['data/As/A0e', 'data/As/A1e', 'data/As/A2e', 'data/As/A3e', $
;                'data/As/A4e', 'data/As/A5e', 'data/As/A6e']
frun= ['data/bs/b0e', 'data/bs/b1e', 'data/bs/b2e', 'data/bs/b3e', $
                'data/bs/b4e', 'data/bs/b5e', 'data/bs/b6e']

fload_xrayinfo, frun, xraygasmass1, xraygasz1, lx_max1, lx_rem1, tx1
fload_lkinfo, frun, lb1



; h
;
;frun= ['data/As/A0', 'data/As/A1', 'data/As/A2', 'data/As/A3', $
;                'data/As/A4', 'data/As/A5', 'data/As/A6']
frun= ['data/bs/b0h', 'data/bs/b1h', 'data/bs/b2h', 'data/bs/b3h', $
                'data/bs/b4h', 'data/bs/b5h', 'data/bs/b6h']

fload_xrayinfo, frun, xraygasmass2, xraygasz2, lx_max2, lx_rem2, tx2
fload_lkinfo, frun, lb2



; f
;
frun= ['data/bs/b0f', 'data/bs/b1f', 'data/bs/b2f', 'data/bs/b3f', $
                'data/bs/b4f', 'data/bs/b5f', 'data/bs/b6f']

fload_xrayinfo, frun, xraygasmass3, xraygasz3, lx_max3, lx_rem3, tx3
fload_lkinfo, frun, lb3



; k
;
frun= ['data/bs/b0k', 'data/bs/b1k', 'data/bs/b2k', 'data/bs/b3k', $
                'data/bs/b4k', 'data/bs/b5k', 'data/bs/b6k']

fload_xrayinfo, frun, xraygasmass4, xraygasz4, lx_max4, lx_rem4, tx4
fload_lkinfo, frun, lb4







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


        ;oplot, lb_100, lx_rem_100, psym=-6, color=200, thick= 3.0, symsize= 1.5
        ;oploterror, lb_100, lx_rem_100, lx_rem_err_100, psym=-3, errcolor=200, color=200, thick= 3.0, errthick= 3.0
   




; -------------------------------------------------------------------------




; extras
; ---------
;xyouts, 0.2, 0.92, 'bulge', /normal, charthick=3.0, size=1.33, color= 0


; Put actual data on there
; -------------------------
;read_data_osullivan_01, loglb, loglx, ttype
;symsize= 0.2
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;oplot, loglb, loglx, psym=8, color= 0
;oplot, loglb, loglx, psym=7, color= 0, symsize= 0.5

;read_data_osullivan_03, loglx, loglb, tempx
;oplot, loglb, loglx, psym=7, color= 0


; a little key
;x0=10.9  &  x1= 12.15
;y0=37.4  &  y1= 38.4
;xyouts, 0.7, 0.28, 'OFP (2001)', color= 0, charthick=2.0, /normal
;oplot, [x0+0.13], [y0+0.75], psym=8, color= 0
;xyouts, 0.7, 0.22, 'OPC (2003)', color= 0, charthick=2.0, /normal
;oplot, [x0+0.13], [y0+0.25], psym=7, color= 0
;oplot, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], psym=-3, color=0


;
;  Mulchaey & Jeltema 2010 field sample 
read_data_mulchaey_10, lk, lx, upperlimit
a= 3.0
ht= a * sin(60.*!PI/180.)
usersym,[-a/2.,0,a/2.,-a/2.],[-ht/2.,ht/2.,-ht/2.,-ht/2.],thick=2.0
oplot, lk, lx, psym=8, color= 0

idx= where(upperlimit eq 1)
if idx(0) ne -1 then begin
  for i=0, n_elements(idx)-1 do begin
	x= lk(idx(i))
	y= lx(idx(i))
	arrow, x, y, x, y-0.25, /data, color= 0, thick=2.0, hsize=-0.5, hthick=2.0
  endfor
endif

x= [1.0,10.0,100.0]
y= 38.90 + 3.92 * (x-11.0)
oplot, x, y, psym=-3, color= 0, thick=2.0, linestyle= 1
xyouts, 0.25, 0.92, 'Field Early-types (   )', color= 0, charthick=2.0, /normal
oplot, [11.15], [42.23], psym=8, color= 0
oplot, [10.57,10.73], [42.22,42.22], psym=-3, color= 0, linestyle= 1, thick= 4.0



;
;  Jeltema 2008 group/cluster sample 
;
read_data_jeltema_08, lk, lx
usersym,[-a/2.,0,a/2.,-a/2.],[-ht/2.,ht/2.,-ht/2.,-ht/2.],thick=2.0,/fill
oplot, lk, lx, psym=8, color= 0

x= [1.0,10.0,100.0]
y= 39.67 + 1.86 * (x-11.0)
oplot, x, y, psym=-3, color= 0, thick=2.0, linestyle= 0
xyouts, 0.25, 0.86, 'Group/Cluster Early-types (   )', color= 0, charthick=2.0, /normal
oplot, [11.325], [41.91], psym=8, color= 0
oplot, [10.57,10.73], [41.91,41.91], psym=-3, color= 0, linestyle= 0, thick= 3.0





; Lx_Lb fit
; ----------
; fit for all E's (OFP '01)
;x= [1.0,10.0,100.0]
;y= 17.98 + 2.17*x
;oplot, x, y, psym=-3, linestyle= 3, thick=2.0, color= 0

; bright E galaxy slope (OPC '03)
;x= [1.0,10.0,100.0]
;y= 11.9 + 2.7*x     ; this is the OPC fit, but B lums are for h=0.5
;y= 12.5 + 2.7*x     ; corrected for this plot because we use B lums from OFP, h=0.75
;oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0





; discrete source Lx  (Ciotti et al. 1991)
; -----------------------
;x= [1.0,10.0,100.0]
;y= 29.45 + x
;oplot, x, y, psym=-3, linestyle= 1, thick=2.0, color= 0




; Data key
; ----------
x0= 10.15  &  y0= 40.95   & dy= 0.255
xn= 0.20  &  yn= 0.68
;usersym,[-1,-1,1,1,-1],[-1,1,1,-1,-1],thick=4.0,/fill
;oplot, [x0], [y0], psym=-8, color=100, thick= 3.0, symsize= 1.5
;xyouts, xn, yn, '5%', color=100, charthick=2.0, /normal, size= 1.2
;----
;oplot, [x0], [y0-dy], psym=-2, color=50, thick= 3.0, symsize= 1.5
;xyouts, xn, yn-0.05, '20%', color=50, charthick=2.0, /normal, size= 1.2
;----
;f = (2.0*!pi/16.0)*findgen(17)
;usersym,cos(f),sin(f),/fill
;oplot, [x0], [y0-2*dy], psym=-8, color=150, thick= 3.0, symsize= 1.5
;xyouts, xn, yn-0.10, '40%', color=150, charthick=2.0, /normal, size= 1.2
;----
;oplot, [x0], [y0-3*dy], psym=-6, color=200, thick= 3.0, symsize= 1.5
;xyouts, xn, yn-0.15, '80%', color=200, charthick=2.0, /normal, size= 1.2

; enclosing box
;oplot, [x0-0.15,x0-0.15,x0+0.65,x0+0.65,x0-0.15], $
;	[y0-1.45,y0+0.2,y0+0.2,y0-1.45,y0-1.45], psym=-3, color=0


xyouts, 0.25, 0.80, 'Simulated Merger Remnants (   )', color= 0, charthick=2.0, /normal
oplot, [11.35], [41.58], psym=2, color= 150, symsize= 1.5, thick=4.0



; helpful info
; --------------

; bh seed mass
;xyouts, 8.8, 42.6, 'BH seed mass', color= 0, charthick=1.2, /data
;arrow, 8.8, 42.5, 8.8, 42.0, COLOR=0, THICK=3.0, hthick=3.0, /data


;xyouts, 0.22, 0.22, 'Without discrete source contribution', /normal, color= 0, size=1.4

device, /close


end








;=====================================================================================













