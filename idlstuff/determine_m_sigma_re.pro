;==================================================================
;==================================================================
;
;      Velocity Dispersion Plotting Routines
;
;
;==================================================================
;==================================================================



; need
;
;   .run time_m_sigma_re
;   .run bh_details
;
;

pro do_all_plots, junk

;time_m_sigma_re, frun
;generate_bh_mass_file_frombh, frun
plot_m_sigma, frun
plot_sigma, frun
plot_re, frun

end








;=====================================================================







;-------------------------------------------------
;
;-------------------------------------------------
pro plot_m_sigma, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_m_sigma, junk"
	print, " "
	print, " "
	return
endif

filename='msigma.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4




xaxistitle = "!7r!6!N(km s!E-1!N)"
yaxistitle = "!6M!DBH!N (M!D!9n!6!N)"
xmax = 600
xmin = 30
ymax = 8.0e+9
ymin = 2.0e+5


;---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
	/xlog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; -----------------------------------


;do_one_msig, "/raid4/tcox/vc3vc3e_2"


; -----------------------------------

;do_one_msig, "data1/priya/v3b1v3b1_e", clr= 0, /show_remnant
;do_one_msig, "data1/priya/v3b1v3b1_f", clr= 0, /show_remnant
;do_one_msig, "data1/priya/v3b1v3b1_k", clr= 0, /show_remnant



do_clr1= 1
;do_clr1= 0
if do_clr1 eq 1 then begin
	do_one_msig, "data1/priya/v3b1v3b1_e", clr= 0, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b1v3b1_f", clr= 0, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b1v3b1_k", clr= 0, /show_remnant, /show_to_from

	do_one_msig, "data1/priya/v3b1v3b4_e", clr= 50, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b1v3b4_f", clr= 50, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b1v3b4_k", clr= 50, /show_remnant, /show_to_from

	do_one_msig, "data1/priya/v3b2v3b1_e", clr= 50, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b2v3b1_f", clr= 50, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b2v3b1_k", clr= 50, /show_remnant, /show_to_from

	do_one_msig, "data1/priya/v3b2v3b2_e", clr= 100, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b2v3b2_f", clr= 100, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b2v3b2_k", clr= 100, /show_remnant, /show_to_from

	do_one_msig, "data1/priya/v3b2v3b4_e", clr= 100, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b2v3b4_f", clr= 100, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b2v3b4_k", clr= 100, /show_remnant, /show_to_from

	do_one_msig, "data1/priya/v3b3v3b1_e", clr= 150, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b3v3b1_f", clr= 150, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b3v3b1_k", clr= 150, /show_remnant, /show_to_from

	do_one_msig, "data1/priya/v3b3v3b2_e", clr= 150, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b3v3b2_f", clr= 150, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b3v3b2_k", clr= 150, /show_remnant, /show_to_from

	do_one_msig, "data1/priya/v3b3v3b3_e", clr= 220, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b3v3b3_f", clr= 220, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b3v3b3_k", clr= 220, /show_remnant, /show_to_from

	do_one_msig, "data1/priya/v3b3v3b4_e", clr= 150, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b3v3b4_f", clr= 150, /show_remnant, /show_to_from
	do_one_msig, "data1/priya/v3b3v3b4_k", clr= 150, /show_remnant, /show_to_from
endif


; -----------------------------------


do_manual= 0
;do_manual= 1
if do_manual eq 1 then begin
	; do some manually

	; vc3vc3e_2
	sig= [156.3]
	bh= [0.00556257]
	oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-2, color= 0
	xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "std", /data, charthick=3.0, size=1, color= 150

	; sb8BH
	sig= [125.5]
	bh= [0.00107389]
	oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-5, color= 0
	;xyouts, sig+10, bh * 1e+10 / 0.7, "sb8", /data, charthick=3.0, size=1, color= 150
	xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=2.0, !8v!6=837", /data, charthick=3.0, size=1, color= 50

	; sb10BH
	sig= [142.7]
	bh= [0.00285036]
	oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-5, color= 0
	;xyouts, sig+10, bh * 1e+10 / 0.7, "sb10", /data, charthick=3.0, size=1, color= 150
	xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=0.5, !8v!6=837", /data, charthick=3.0, size=1, color= 50

	;sb10BHtr1
	sig= [141.0]
	bh= [0.00388664]
	oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-5, color= 0

	;sb10BHtr2
	sig= [142.3]
	bh= [0.00433478]
	oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-5, color= 0

	; sb13BH
	sig= [158.8]
	bh= [0.158933]
	oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-5, color= 0
	;xyouts, sig+10, bh * 1e+10 / 0.7, "sb13", /data, charthick=3.0, size=1, color= 150
	xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=2.0, !8v!6=209", /data, charthick=3.0, size=1, color= 50
endif


; -----------------------------------

mblackhole_sigma_relation, 1

;xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0
;if bhmsg NE '' then xyouts, 0.25, 0.80, bhmsg, /normal, charthick=3.0, size=1.33, color=0




; -----------------------------------

device, /close


end




; -----------------------------------




;   m-sigma relation
; -----------------------------------
pro mblackhole_sigma_relation, junk

sigm=[0.1, 1.0, 10.0, 100.0, 1000.0]
mbh= (1.2e+8)*((sigm/200.0)^(3.75))
oplot, sigm, mbh, psym=-3, linestyle=1, thick=3.0, color=0

end



;===================================================================================



pro do_one_msig, frun, $
		clr=clr, $
		show_merger=show_merger, $
		show_remnant=show_remnant, $
		show_track=show_track, $
		show_pre=show_pre, $
		show_post=show_post, $
		show_to_from=show_to_from

print, "processing: ", frun
print, "clr: ", clr

cmd= "ls "+frun+"/bh_*.txt"
spawn, cmd, result

; assume there are two BHs
read_bh_file, frun, 1, time1, bhm1, pm1, mdot1, mdotedd1, bhfile=result[0], mergertime=mergertime1
idx= where(bhm1 gt 0.0)
N1= n_elements(idx)
time1= time1(idx)
bhm1= bhm1(idx) * 1.0e+10 / 0.7
if mergertime1 gt 0 then bhmergertime= mergertime1 else bhmergertime= 1.25
if n_elements(result) gt 1 then begin
	read_bh_file, frun, 1, time2, bhm2, pm2, mdot2, mdotedd2, bhfile=result[1], mergertime=mergertime2
	idx= where(bhm2 gt 0.0)
	N2= n_elements(idx)
	time2= time2(idx)
	bhm2= bhm2(idx) * 1.0e+10 / 0.7 
	if mergertime2 gt 0 then bhmergertime= mergertime2
endif else begin
	N2= 0
end

;
; we'd want #1 to be the final BH
;
if N2 gt N1 then begin
	N_temp= N2
	time_temp= time2
	bhm_temp= bhm2
	N2= N1
	time2= time1
	bhm2= bhm1
	N1= N_temp
	time1= time_temp
	bhm1= bhm_temp
endif

read_sigma_file, frun, time, Asigxy, Asigxz, Asigyz, Asigavg, Asigerr, $
                                Bsigxy, Bsigxz, Bsigyz, Bsigavg, Bsigerr, $
                                gsigavg, gsigerr


;
;
idx= where(Asigavg gt 0.0)
if idx(0) ne -1 then begin
	if n_elements(idx) eq N1 then sig1= Asigavg(idx) else sig2= Asigavg(idx)
endif

;
;
idx= where(Bsigavg gt 0.0)
if idx(0) ne -1 then begin
	if n_elements(idx) eq N1 then sig1= Bsigavg(idx) else sig2= Bsigavg(idx)
endif



;
;
;
;if n_elements(time1) ne n_elements(time) then begin
;	print, " "
;	print, " n(time)= ", n_elements(time)
;	print, " n(time1)= ", n_elements(time1)
;	print, " PROBLEM: sigma.txt and bh_######.txt don't have identical elements"
;	print, " "
;endif
;

print, " n(time)= ", n_elements(time)
print, " n(Asigavg)= ", n_elements(Asigavg)
print, " n(Bsigavg)= ", n_elements(Bsigavg)
print, " n(time1)= ", n_elements(time1)
print, " n(bhm1)= ", n_elements(bhm1)
print, " n(sig1)= ", n_elements(sig1)
print, " n(time2)= ", n_elements(time2)
print, " n(bhm2)= ", n_elements(bhm2)
print, " n(sig2)= ", n_elements(sig2)


;
; first BH
;
if keyword_set(show_track) then oplot, sig1, bhm1, thick= 1.0, color=clr, psym=-3

if keyword_set(show_pre) then begin
	idx= where(time1 lt bhmergertime)
	oplot, sig1(idx), bhm1(idx), thick=1.0, color= clr, psym=-3
endif

if keyword_set(show_post) then begin
	idx= where(time1 ge bhmergertime)
	if idx(0) ne -1 then oplot, sig1(idx), bhm1(idx), thick=1.0, color= clr, psym=-3
endif

if keyword_set(show_merger) then begin
	idx= where(time1 ge bhmergertime)
	if idx(0) ne -1 then oplot, [sig1(idx[0])], [bhm1(idx[0])], thick=2.0, color= clr, psym=-5
endif

if keyword_set(show_to_from) then oplot, [sig1[0],sig1(N1-1)], [bhm1[0],bhm1(N1-1)], thick=1.0, color= clr, psym=-6
if keyword_set(show_remnant) then oplot, [sig1(N1-1)], [bhm1(N1-1)], thick=5.0, color= clr, psym=-2

print, "bhm1[0]= ", bhm1[0]



;
; second BH (this is assumed to never be the final BH)
;
if N2 gt 0 then begin
    if keyword_set(show_track) then oplot, sig2, bhm2, thick= 1.0, color=clr, psym=-3

    if keyword_set(show_pre) then begin
	idx= where(time1 lt bhmergertime)
	oplot, sig2(idx), bhm2(idx), thick=1.0, color= clr, psym=-3
    endif

    ;if keyword_set(show_post) then begin
    ;	idx= where(time1 ge bhmergertime)
    ;	if idx(0) ne -1 then oplot, sig2(idx), bhm2(idx), thick=1.0, color= clr, psym=-3
    ;endif

    ;if keyword_set(show_merger) then begin
    ;	idx= where(time1 ge bhmergertime)
    ;	if idx(0) ne -1 then oplot, [sig2(idx[0])], [bhm2(idx[0])], thick=2.0, color= clr, psym=-5
    ;endif

    ;if keyword_set(show_to_from) then oplot, [sig2[0],sig1(N1-1)], [bhm2[0],bhm1(N1-1)], thick=1.0, color= clr, psym=-6
    ;if keyword_set(show_to_from) then oplot, [sig2[0],sig2(N2-1)], [bhm2[0],bhm1(N1-1)], thick=1.0, color= clr, psym=-6
    ;if keyword_set(show_remnant) then oplot, [sig2(N2-1)], [bhm2(N2-1)], thick=5.0, color= clr, psym=-2
    print, "bhm2[0]= ", bhm2[0]
endif


;arrow, 132.0, 6.0e+5, 135.0, 2.0e+6, COLOR=0, THICK=3.0, hthick=3.0, /data
;arrow, 160.0, 8.0e+6, 210.0, 1.3e+7, COLOR=0, THICK=3.0, hthick=3.0, /data



end








;===================================================================================







pro plot_sigma, frun


if not keyword_set(frun) then begin
	print, " "
	print, " plot_sigma, frun"
	print, " "
	print, " "
	return
endif

filename=frun+'/sigma.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


read_sigma_file, frun, time, Asigxy, Asigxz, Asigyz, Asigavg, Asigerr, $
				Bsigxy, Bsigxz, Bsigyz, Bsigavg, Bsigerr, gsigavg, gsigerr

time= time/0.7    ; h=0.7

yaxistitle = "!7r!3!N(km s!E-1!N)"
xaxistitle = "Time (Gyr)"
xmax = max(time)
xmin = 0
ymax = 400
ymin = 0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; center 1
thispsym= 6
thiscolor= 70
oplot, time, Asigavg, thick=4.0, psym=-thispsym, color=thiscolor
oplot, time, Asigavg+Asigerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
oplot, time, Asigavg-Asigerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1

;center 2
show_gal2= 0
if show_gal2 eq 1 then begin
	thispsym=2
	thiscolor= 70
	oplot, time, Bsigavg, thick=4.0, psym=-thispsym, color=thiscolor
	oplot, time, Bsigavg+Bsigerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
	oplot, time, Bsigavg-Bsigerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
endif


xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0
;if bhmsg NE '' then xyouts, 0.25, 0.80, bhmsg, /normal, charthick=3.0, size=1.33, color=0

show_gas_too= 1
if show_gas_too eq 1 then begin
	thispsym= 3
	thiscolor= 150
	oplot, time, gsigavg, thick=4.0, psym=-thispsym, color=thiscolor
	oplot, time, gsigavg+gsigerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
	oplot, time, gsigavg-gsigerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
endif

;if timemerge gt 0 then begin
;        arrow, timemerge, 1.0e-4, timemerge, 10^(-3.5), COLOR=0, THICK=3.0, hthick=3.0, /data
;endif


device, /close




end





;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro plot_sigma_comparison, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_sigma_comparison, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='sigma.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



yaxistitle = "!7r!6!N (km s!E-1!N)"
xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
;xmax = max(time)
xmax = 2.0
;xmax = 1.5
;xmin = 1.0
xmin = 0.5

;ymax = 400
ymax = 325
;ymax = 275
ymin = 75
;ymin = 0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; Galaxy 1
;----------
;frun="pool/vc3vc3"
;frun="pool/vc3bvc3b"
frun="/raid4/tcox/vc3vc3e_2"
read_sigma_file, frun, time, Asigxy, Asigxz, Asigyz, Asigavg, Asigerr, $
				Bsigxy, Bsigxz, Bsigyz, Bsigavg, Bsigerr, gsigavg, gsigerr

;oplot, time, Asigavg, thick=3.0, psym=-2, color= 150
oplot, time, Asigavg, thick=3.0, psym=-5, color= 150
oploterror, time, Asigavg, Asigerr, psym=-3, errcolor=150, color=150, errthick= 3.0, thick= 3.0


xyouts, 0.7, 0.86, 'stars', /normal, charthick=2, size=1.33, color=150

xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0




; gas
;----------
;frun="pool/vc3vc3_wBH"
;frun="pool/vc3bvc3b_wBH"
;frun="/raid4/tcox/vc3vc3e_no"
;read_sigma_file, frun, time, Asigxy, Asigxz, Asigyz, Asigavg, Asigerr, $
;				Bsigxy, Bsigxz, Bsigyz, Bsigavg, Bsigerr, gsigavg, gsigerr, $
;				sigfile='sigma_gas.txt'

;oplot, time, Asigavg, thick=3.0, psym=-6, color= 100
;oploterror, time, Asigavg, Asigerr, psym=-3, errcolor=100, color=100, errthick= 3.0, thick= 3.0

;xyouts, 0.7, 0.86, 'H!7a!6', /normal, charthick=2, size=1.33, color=100





;------------------
overplot_desika= 0
;overplot_desika= 1
if overplot_desika eq 1 then begin
	read_desika_FWHM_file, "COData/e2_z0.1.dat", time, FWHMxy, FWHMxz, FWHMyz
	FWHMxy= FWHMxy / 2.3548
	FWHMxz= FWHMxz / 2.3548
	FWHMyz= FWHMyz / 2.3548
		n_fwhm= n_elements(FWHMxy)
		FWHM_avg= fltarr(n_fwhm)
		FWHM_sig= fltarr(n_fwhm)
		for i=0,n_fwhm-1 do begin
		  ddata= [FWHMxy[i], FWHMxz[i], FWHMyz[i]]
		  FWHM_avg[i]= mean(ddata)
		  FWHM_sig[i]= sqrt(variance(ddata))
		endfor

	oplot, time, FWHM_avg, thick=3.0, psym=-2, color= 50
	oploterror, time, FWHM_avg, FWHM_sig, psym=-3, errcolor=50, color=50, errthick= 3.0, thick= 3.0

	xyouts, 0.7, 0.80, 'CO (2-1)', /normal, charthick=2, size=1.33, color=50
endif


;------------------
device, /close




end




;--------------------------------------
;   at various radii
;----------------------------------------
pro plot_many_sigma, frun


if not keyword_set(frun) then begin
	print, " "
	print, " plot_many_sigma, frun"
	print, " "
	print, " "
	return
endif

filename=frun+'/manysigma.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


sigfile= frun+'/sigma_gas.txt'
;sigfile= frun+'/sigma_stars.txt'
read_many_sigma_file, sigfile, time, sig1, sig2, sig3, sig4, sig5, sig6

;time= time/0.7    ; h=0.7

yaxistitle = "!7r!3!N(km s!E-1!N)"
;xaxistitle = "Time (Gyr)"
xaxistitle = "Time (Gyr h!E-1!N)"
xmax = max(time)
xmin = 0
ymax = 225
ymin = 25


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



thispsym= 3

oplot, time, sig1, thick=4.0, psym=-thispsym, color=50
oplot, time, sig2, thick=4.0, psym=-thispsym, color=70
oplot, time, sig3, thick=4.0, psym=-thispsym, color=90
oplot, time, sig4, thick=4.0, psym=-thispsym, color=110
oplot, time, sig5, thick=4.0, psym=-thispsym, color=150
oplot, time, sig6, thick=4.0, psym=-thispsym, color=200


;xyouts, 0.72, 0.82, 'Stars', /normal, charthick=4.0, color= 0, size=1.6
;xyouts, 0.58, 0.46, 'Aperture Size (kpc)', /normal, charthick=3.0, color= 0, size=1.1
;xyouts, 0.74, 0.42, '1.0', /normal, charthick=3.0, color= 50, size=1.1
;xyouts, 0.74, 0.38, '1.0', /normal, charthick=3.0, color= 50, size=1.1
;xyouts, 0.74, 0.34, '2.0', /normal, charthick=3.0, color= 70, size=1.1
;xyouts, 0.74, 0.30, '3.0', /normal, charthick=3.0, color= 90, size=1.1
;xyouts, 0.74, 0.26, '4.0', /normal, charthick=3.0, color= 110, size=1.1
;xyouts, 0.74, 0.22, '5.0', /normal, charthick=3.0, color= 150, size=1.1
;xyouts, 0.74, 0.18, '6.0', /normal, charthick=3.0, color= 200, size=1.1

xyouts, 0.22, 0.85, 'Gas', /normal, charthick=4.0, color= 0, size=1.6
xyouts, 0.58, 0.86, 'Aperture Size (kpc)', /normal, charthick=3.0, color= 0, size=1.1
xyouts, 0.74, 0.82, '1.0', /normal, charthick=3.0, color= 50, size=1.1
xyouts, 0.74, 0.78, '1.0', /normal, charthick=3.0, color= 50, size=1.1
xyouts, 0.74, 0.74, '2.0', /normal, charthick=3.0, color= 70, size=1.1
xyouts, 0.74, 0.70, '3.0', /normal, charthick=3.0, color= 90, size=1.1
xyouts, 0.74, 0.66, '4.0', /normal, charthick=3.0, color= 110, size=1.1
xyouts, 0.74, 0.62, '5.0', /normal, charthick=3.0, color= 150, size=1.1


device, /close




end















;--------------------------------------
;
;----------------------------------------
pro plot_re, frun


if not keyword_set(frun) then begin
	print, " "
	print, " plot_re, frun"
	print, " "
	print, " "
	return
endif

filename=frun+'/re.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


read_Re_file, frun, time, re, reerr, re2, re2err, gre, greerr

yaxistitle = "R!De!N(kpc)"
xaxistitle = "Time (Gyr)"
xmax = max(time)
xmin = 0
ymax = 10
ymin = 0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; center 1
thispsym= 6
thiscolor= 70
oplot, time, re, thick=3.0, psym=-thispsym, color=thiscolor
;errplot, time, re-reerr, re+reerr      ; produces black and white
oplot, time, re+reerr, thick=1.0, psym=-3, linestyle= 1, color=thiscolor
oplot, time, re-reerr, thick=1.0, psym=-3, linestyle= 1, color=thiscolor

; gas
; -----
show_gas_too= 1
if show_gas_too eq 1 then begin
	thispsym= 3
	thiscolor= 150
	oplot, time, gre, thick=4.0, psym=-thispsym, color=thiscolor
	oplot, time, gre+greerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
	oplot, time, gre-greerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
endif

; extras
; -------
xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0

; done
; ------
device, /close


end






;==================================================================
;==================================================================









;========================================
;========================================
;
;
; procedure reads phil's list of 
; directories to process
;
;
pro read_phils_list, philslist


; we read the 03 data, note that this is mainly
; informational data, so sigma, R_e, etc.  L_x, T_x,
; are all enclosed above.
;philsfile= '/home/tcox/sigmas_to_dump.txt'
;philsfile= '/home/tcox/sigmas_1.txt'
philsfile= '/home/tcox/sigmas_2.txt'
;philsfile= '/home/tcox/sigmas_1a.txt'

spawn, "wc "+philsfile,result
lines=long(result)
datalines03=lines(0)
sigma= fltarr(datalines03)

openr, 1, philsfile
junk=''

philslist= strarr(lines(0))

; read the data
for i=0,lines(0)-1 do begin

        readf, 1, junk

        ;tempjunk= strsplit(junk,/extract,count=count)
        ;philslist(i)= tempjunk(0)

	philslist(i)= strcompress(junk,/remove_all)

endfor

close, 1




end




;========================================
;========================================

pro do_phils_list, junk

	read_phils_list, philslist
	n_philslist= n_elements(philslist)

	for i=0,n_philslist-1 do begin
	   print, "--------------------------------------------------------------"
	   print, ""
	   print, "   xxxxxx  processing file:", philslist[i]
	   print, ""
	   print, "--------------------------------------------------------------"
	   time_m_sigma_re, philslist[i]
	   ;print, philslist[i]
	endfor

	print, "--------------------------------------------------------------"
	print, ""
	print, ""
	print, ""
	print, "--------------------------------------------------------------"

end




;========================================
;========================================




