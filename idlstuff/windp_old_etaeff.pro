
; =================================================
;
;  This section analyzes the properties of the
;  hot, outflowing gas.
;
;
;
;
; ==================================================
pro time_etaeff, junk


if not keyword_set(junk) then begin
        print, " "
        print, " time_etaeff, junk"
        print, " "
        print, " "
        print, " "
        print, " "
        return
endif


frun= junk
;frun= "/raid4/tcox/vc3vc3e_2"
;frun= "/raid4/tcox/sbw/sb8"
;frun= "/raid4/tcox/sbw/sb10"
;frun= "/raid4/tcox/sbw/sb20"
;frun= "/raid4/tcox/sbw/sb10BH"
;frun= "/raid4/tcox/sbw/sb22BH"


; this assumes they are ordered
;  0 through


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])



; --------------
; Parameters
; --------------

;r_wind= 150.0      ; used in wind_gtx
r_wind= 50.0      ; used in wind_gtx

; std vc3c mergers
;g1_sid= 1L
;g1_npart= 200001L
;g2_sid= 200002L
;g2_npart= 200001L

; std vc3c mergers
g1_sid= 1L
g1_npart= 550001L
g2_sid= 550002L
g2_npart= 550001L



; ---------------
; Arrays
; ---------------
time= fltarr(nsnaps)

center_1= fltarr(nsnaps,3)
center_2= fltarr(nsnaps,3)

mass_gtx_1=  fltarr(nsnaps)
mass_gtx_2=  fltarr(nsnaps)
z_gtx_1=  fltarr(nsnaps)
z_gtx_2=  fltarr(nsnaps)

mass_gtx=  fltarr(nsnaps)
z_gtx=  fltarr(nsnaps)
outflowr= fltarr(nsnaps)
etaeff= fltarr(nsnaps)

sfr_inst= fltarr(nsnaps)
sfr_avg= fltarr(nsnaps)



;-------------------------------------------
;   Get SFR rate from txt - for each file
;-------------------------------------------
open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass




; ----------------------------------------
; This part loops through the snapshots
; and compiles various information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"


        ; open snapshot
        ;ok=fload_snapshot(frun,snapnum)
        ok=fload_snapshot_bh(frun,i)


        ; what time is it?
        time[i]= fload_time(1)
        TTime= float(time[i])


       ; instead use BH positions
        ; ------------------------
        use_bh_positions= 1
        if use_bh_positions eq 1 then begin

           ; get blackhole id's
           if i eq 0 then begin
                bhid= fload_blackhole_id(1)
                bhid1= bhid[0]
                bhid2= bhid[1]
                print, "Blackhole ID's: ", bhid1, bhid2
           endif

           if fload_npart(5) gt 1 then begin
                center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
                center_2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
                ;center_1= fload_blackhole_xyz('xyz',idtofollow=bhid1)
                ;center_2= fload_blackhole_xyz('xyz',idtofollow=bhid2)
           endif else begin
                bhid= fload_blackhole_id(1)
                center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
                ;center_1= fload_blackhole_xyz('xyz',idtofollow=bhid)
                center_2= center_1
           endelse
        endif



	; gas properties
	; -----------------

	;  1
	; ===
	r1= fload_1gal_gas_xyz('r', g1_sid, g1_npart, center=center_1)
	gm1= fload_1gal_gas_mass(1, g1_sid, g1_npart)
	gz1= fload_1gal_gas_metallicity(1, g1_sid, g1_npart)

	; greater than r_wind
	gtr_idx= where(r1 gt r_wind)
	if gtr_idx(0) ne -1 then m1_gtx= total(gm1(gtr_idx)) else m1_gtx= 0
	if gtr_idx(0) ne -1 then z1_gtx= mean(gz1(gtr_idx)) else z1_gtx= 0


	;  2
	; ===
	r2= fload_1gal_gas_xyz('r', g2_sid, g2_npart, center=center_2)
	gm2= fload_1gal_gas_mass(1, g2_sid, g2_npart)
	gz2= fload_1gal_gas_metallicity(1, g2_sid, g2_npart)

	; greater than r_wind
	gtr_idx= where(r2 gt r_wind)
	if gtr_idx(0) ne -1 then m2_gtx= total(gm2(gtr_idx)) else m2_gtx= 0
	if gtr_idx(0) ne -1 then z2_gtx= mean(gz2(gtr_idx)) else z2_gtx= 0


	mass_gtx_1[i]= m1_gtx
	mass_gtx_2[i]= m2_gtx
	z_gtx_1[i]= z1_gtx
	z_gtx_2[i]= z2_gtx

	mass_gtx[i]= m1_gtx + m2_gtx
	z_gtx[i]= z1_gtx + z2_gtx



	;  SFR
	; ----------
        idx_gtr_snaptime= where(sfrtime ge time[i])
        idx_curr_snaptime= idx_gtr_snaptime[0]

        ; instantaneous sfr
        ;sfr_inst[i]= sfrsfr(idx_curr_snaptime)

        ; set time window
        dt= 0.01
        idx_window= where((sfrtime ge (time[i]-dt)) and (sfrtime le (time[i]+dt)))
        sfr_dt= sfrsfr(idx_window)
        sfr_avg[i]= mean(sfr_dt)


endfor




t_last= 0 & t_next= 0 & t_i= 0
m_last= 0 & m_next= 0 & m_i= 0

for i=0,nsnaps-1 do begin

	if i eq 0 then t_last= -1 else t_last= time[i-1]
	t_i= time[i]
	if i eq (nsnaps-1) then t_next= -1 else t_next= time[i+1]

	if i eq 0 then m_last= mass_gtx[i+1] else m_last= mass_gtx[i-1]
	m_i= mass_gtx[i]
	if i eq (nsnaps-1) then m_next= mass_gtx[i-1] else m_next= mass_gtx[i+1]

	;dm_wind= 0.25 * (m_next + 2.0*m_i + m_last)
	dm_wind= 0.5 * (m_next - m_last)
	dt= 0.5 * (t_next - t_last)

	outflowr[i]= 10.0*dm_wind/dt    ; in 10^10 msolar / Gyr   = 10 Msolar/Yr



	if sfr_avg[i] le 0.1 then sfr_avg[i]= 0.1
	print, "i= ",i, " outflow/sfr= ", outflowr[i], sfr_avg[i]
	etaeff[i]= outflowr[i]/sfr_avg[i]

endfor




; ----------------------------------------
; write to file

;openw, 1, frun+'/sbwind_150.txt', ERROR=err
openw, 1, frun+'/sbwind.txt', ERROR=err
        
printf, 1, "#   sbwind.txt"     
printf, 1, "#    (all units are gadget units),   r_wind=",r_wind
printf, 1, "#               "
printf, 1, "#  time      gtx#1    gtx#2   totgtx   wind z  outflow  sfr_avg      etaeff "               
printf, 1, "#(Gyr/h)       (m)      (m)      (m)      (z)     rate (msun/yr)  "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '("  ",F6.3,"  ", 6(F8.4," "),"   ",F8.4)', $
                time[i], mass_gtx_1[i], mass_gtx_2[i], mass_gtx[i], $
                        z_gtx[i], outflowr[i], sfr_avg[i], etaeff[i]
endfor      
close, 1    
            
            
        

; ----------------------------------------
        
print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"
        




end






;========================================================================================







; ----------------------------
;  Read wind.txt file
; ----------------------------
pro read_sbwind_file, frun, time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 

filename= frun+'/sbwind.txt'

spawn, "wc "+filename,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then filedata= fltarr(8,lines)

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
mass_gtx_1= filedata[1,*]
mass_gtx_2= filedata[2,*]
mass_gtx= filedata[3,*]
z_gtx= filedata[4,*]
outflowr= filedata[5,*]
sfr_avg= filedata[6,*]
etaeff= filedata[7,*]


end







;========================================================================================






;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro windc, junk


if not keyword_set(junk) then begin
	print, " "
	print, " windc, junk"
	print, " "
	print, " "
	return
endif

filename='etaeffwind.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;yaxistitle = "!6Log !7g!6!Deff!N"
yaxistitle = "!7g!6!Deff!N"
;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"

;xmax = max(time)
xmax = 4.25
;xmax = 3.0
;xmax = 2.0
xmin = 0

;ymax = 150.0
ymax = 25.0
;ymax = 15.0
;ymax = 10.0
;ymin = 7.0
;ymin = 6.0
ymin = 2.5e-2


;---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata





;read_sbwind_file, "/raid4/tcox/vc3vc3e_2", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
;read_sbwind_file, "/raid4/tcox/sbw/sb8", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
;time= time/0.7
;idx=where(etaeff gt 0)
;oplot, time(idx), etaeff(idx), thick=6.0, color=100, psym=-3, linestyle=0

read_sbwind_file, "/raid4/tcox/sbw/sb10", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
; ---
idx= where(time lt 2.0)
print, 'mean etaeff=', mean(etaeff)
print, 'mean etaeff (t<2)=', mean(etaeff(idx))
; ---
idx=where(etaeff gt 0)
oplot, time(idx), etaeff(idx), thick=6.0, color=150, psym=-5, linestyle=0
oplot, [0.3,0.6], [11,11], thick=6.0, color=150, psym=-5, linestyle=0
xyouts, 0.33, 0.87, 'standard', color=150, charthick=3.0, size=1.3, /normal

read_sbwind_file, "/raid4/tcox/sbw/sb10BH", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
; ---
idx= where(time lt 2.0)
print, 'mean etaeff=', mean(etaeff)
print, 'mean etaeff (t<2)=', mean(etaeff(idx))
; ---
idx=where(etaeff gt 0)
oplot, time(idx), etaeff(idx), thick=6.0, color=50, psym=-3, linestyle=0
oplot, [0.3,0.6], [6,6], thick=6.0, color=50, psym=-5, linestyle=0
xyouts, 0.33, 0.80, 'with BH', color=50, charthick=3.0, size=1.3, /normal




; extras

x=[xmin,xmax]
y=[0.5,0.5]
oplot, x, y, psym=-3, linestyle=2, color=0


device, /close




end







;==========================================================================






;===================================
;
;  Outflow Information
;
;===================================


;
pro outflow_info, junk


if not keyword_set(junk) then begin
	print, " "
	print, " outflow_info, junk"
	print, " "
	print, " "
	return
endif

filename='outflowinfo.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=12.0, newysize= 20.0


;---------------------------



xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
;xaxistitle = "!6Time (Gyr)"
;xmax = 4.25
xmax = 3.0
xmin = 0.0

yaxistitle = "!6Log Outflow Mass (M!D!9n!6!N)"
ymax = 10.5
ymin = 9.1



;---------------------------

x0= 0.18
x1= 0.98

y0= 0.08
y1= 0.53
y2= 0.98

;---------------------------


!p.position= [x0, y1, x1, y2]


plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        ;xtitle=xaxistitle, $
	xtickformat='(a1)', $
        ytitle=yaxistitle, $
        /nodata



; ----------------------------


;read_sbwind_file, "/raid4/tcox/sbw/sb8", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
;mass_gtx= alog10(mass_gtx * 1.0d+10 / 0.7)
;oplot, time, mass_gtx, thick=4.0, color=0, psym=-3, linestyle=0

read_sbwind_file, "/raid4/tcox/sbw/sb10", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
mass_gtx= alog10(mass_gtx * 1.0d+10 / 0.7)
oplot, time, mass_gtx, thick=4.0, color=50, psym=-3, linestyle=0

read_sbwind_file, "/raid4/tcox/sbw/sb10BH", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
mass_gtx= alog10(mass_gtx * 1.0d+10 / 0.7)
oplot, time, mass_gtx, thick=6.0, color=150, psym=-3, linestyle=0

read_sbwind_file, "/raid4/tcox/vc3vc3e_2", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
mass_gtx= alog10(mass_gtx * 1.0d+10 / 0.7)
oplot, time, mass_gtx, thick=6.0, color=100, psym=-3, linestyle=0


; ----------------------------


yaxistitle = "!6Outflow Rate (M!D!9n!6!N Yr!E-1!N)"
ymax = 10.0
ymin = -1.0


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, /noerase, $
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




; ----------------------------

;read_sbwind_file, "/raid4/tcox/sbw/sb8", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
;oplot, time, outflowr, thick=4.0, color=0, psym=-3, linestyle=0

read_sbwind_file, "/raid4/tcox/sbw/sb10", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
oplot, time, outflowr, thick=4.0, color=50, psym=-3, linestyle=0

read_sbwind_file, "/raid4/tcox/sbw/sb10BH", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
oplot, time, outflowr, thick=8.0, color=150, psym=-3, linestyle=0

read_sbwind_file, "/raid4/tcox/vc3vc3e_2", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
oplot, time, outflowr, thick=8.0, color=100, psym=-3, linestyle=0



; ----------------------------
device, /close

end








;==========================================================================






;===================================
;
;  Figure for Wind Paper
;
;===================================


;
pro windfig, junk


if not keyword_set(junk) then begin
	print, " "
	print, " windfig, junk"
	print, " "
	print, " "
	return
endif

filename='windfig.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=12.0, newysize= 20.0


;---------------------------

;     ----------------
;     |              |
;     |              |
;     |   Outflow    |
;     |    Rate      |
;     |              |
;     |              |
;     ----------------
;     |              |
;     |              |
;     |   Eta_eff    |
;     |              |
;     |              |
;     |              |
;     ----------------



;---------------------------



;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"
xmax = 4.25
;xmax = 3.0
xmin = 0.0


yaxistitle = "!6Outflow Rate (M!D!9n!6!N Yr!E-1!N)"
ymax = 50.0
ymin = 3.0e-2



;---------------------------

x0= 0.18
x1= 0.98

y0= 0.08
y1= 0.53
y2= 0.98

;---------------------------


!p.position= [x0, y1, x1, y2]


plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        ;xtitle=xaxistitle, $
	xtickformat='(a1)', $
        ytitle=yaxistitle, $
        /nodata



; ----------------------------


read_sbwind_file, "/raid4/tcox/sbw/sb10", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
idx= where(outflowr gt 0.0)
oplot, time(idx), outflowr(idx), thick=6.0, color=150, psym=-5, linestyle=0

; ----------------------------

read_sbwind_file, "/raid4/tcox/sbw/sb10BH", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
idx= where(outflowr gt 0.0)
oplot, time(idx), outflowr(idx), thick=6.0, color=50, psym=-3, linestyle=0


; ----------------------------


yaxistitle = "!7g!6!Deff!N"
ymax = 25.0
ymin = 2.5e-2


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, /noerase, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; ----------------------------


read_sbwind_file, "/raid4/tcox/sbw/sb10", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
; ---
idx= where(time lt 2.0)
print, 'mean etaeff=', mean(etaeff)
print, 'mean etaeff (t<2)=', mean(etaeff(idx))
; ---
idx=where(etaeff gt 0)
oplot, time(idx), etaeff(idx), thick=6.0, color=150, psym=-5, linestyle=0
oplot, [0.3,0.6], [11,11], thick=6.0, color=150, psym=-5, linestyle=0
xyouts, 0.33, 0.87, 'standard', color=150, charthick=3.0, size=1.3, /normal


; ----------------------------


read_sbwind_file, "/raid4/tcox/sbw/sb10BH", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff  
time= time/0.7
; ---
idx= where(time lt 2.0)
print, 'mean etaeff=', mean(etaeff)
print, 'mean etaeff (t<2)=', mean(etaeff(idx))
; ---
idx=where(etaeff gt 0)
oplot, time(idx), etaeff(idx), thick=6.0, color=50, psym=-3, linestyle=0
oplot, [0.3,0.6], [6,6], thick=6.0, color=50, psym=-5, linestyle=0
xyouts, 0.33, 0.80, 'with BH', color=50, charthick=3.0, size=1.3, /normal




; extras
; ---------
x=[xmin,xmax]
y=[0.5,0.5]
oplot, x, y, psym=-3, linestyle=2, color=0





; ----------------------------
device, /close

end








;==========================================================================






;===================================
;
;  Figure for Wind Paper
;
;===================================


;
pro windfig2, junk


if not keyword_set(junk) then begin
	print, " "
	print, " windfig2, junk"
	print, " "
	print, " "
	return
endif

filename='windfig.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=12.0, newysize= 22.0


;---------------------------

;     ----------------
;     |              |
;     |              |
;     |   SFR        |
;     |              |
;     |              |
;     |              |
;     ----------------
;     |              |
;     |              |
;     |   Outflow    |
;     |    Rate      |
;     |              |
;     |              |
;     ----------------
;     |              |
;     |              |
;     |   Eta_eff    |
;     |              |
;     |              |
;     |              |
;     ----------------



;---------------------------




;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"
xmax = 4.25
;xmax = 3.0
xmin = 0.0


yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
ymax = 300.0
ymin = 7.0e-3



;---------------------------

x0= 0.21
x1= 0.98

y0= 0.09
y1= 0.39
y2= 0.69
y3= 0.99

;---------------------------


!p.position= [x0, y2, x1, y3]


plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        /ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        ;xtitle=xaxistitle, $
        xtickformat='(a1)', $
        ytitle=yaxistitle, $
        /nodata



;---------------------------



process_one_sfr, "sbw/sb10", lcolor=150, lpsym=-5, lthick= 2.0, msg=' ', x0= 0.90, y0= 0.90, h=0.7
oplot, [2.3,2.5], [44,44], thick=6.0, color=150, psym=-5, linestyle=0
xyouts, 0.69, 0.93, 'sb winds', color=150, charthick=3.0, size=1.3, /normal



process_one_sfr, "sbw/sb10BH", lcolor=50, lthick= 2.0, msg=' ', x0= 0.90, y0= 0.90, h=0.7
oplot, [2.25,2.55], [11,11], thick=6.0, color=50, psym=-3, linestyle=0
xyouts, 0.69, 0.89, 'sb winds+BH', color=50, charthick=3.0, size=1.3, /normal



;---------------------------


yaxistitle = "!6Outflow Rate (M!D!9n!6!N Yr!E-1!N)"
ymax = 50.0
ymin = 3.0e-2





;---------------------------

!p.position= [x0, y1, x1, y2]


plot, [1.0],[1.0], psym=-3, /noerase, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        ;xtitle=xaxistitle, $
	xtickformat='(a1)', $
        ytitle=yaxistitle, $
        /nodata



; ----------------------------


read_sbwind_file, "/raid4/tcox/sbw/sb10", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
idx= where(outflowr gt 0.0)
oplot, time(idx), outflowr(idx), thick=6.0, color=150, psym=-5, linestyle=0

; ----------------------------

read_sbwind_file, "/raid4/tcox/sbw/sb10BH", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
idx= where(outflowr gt 0.0)
oplot, time(idx), outflowr(idx), thick=6.0, color=50, psym=-3, linestyle=0


; ----------------------------


yaxistitle = "!7g!6!Deff!N"
ymax = 25.0
ymin = 2.5e-2


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, /noerase, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; ----------------------------


read_sbwind_file, "/raid4/tcox/sbw/sb10", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
; ---
idx= where(time lt 2.0)
print, 'mean etaeff=', mean(etaeff)
print, 'mean etaeff (t<2)=', mean(etaeff(idx))
; ---
idx=where(etaeff gt 0)
oplot, time(idx), etaeff(idx), thick=6.0, color=150, psym=-5, linestyle=0
;oplot, [0.3,0.6], [11,11], thick=6.0, color=150, psym=-5, linestyle=0
;xyouts, 0.33, 0.87, 'standard', color=150, charthick=3.0, size=1.3, /normal


; ----------------------------


read_sbwind_file, "/raid4/tcox/sbw/sb10BH", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff  
time= time/0.7
; ---
idx= where(time lt 2.0)
print, 'mean etaeff=', mean(etaeff)
print, 'mean etaeff (t<2)=', mean(etaeff(idx))
; ---
idx=where(etaeff gt 0)
oplot, time(idx), etaeff(idx), thick=6.0, color=50, psym=-3, linestyle=0
;oplot, [0.3,0.6], [6,6], thick=6.0, color=50, psym=-5, linestyle=0
;xyouts, 0.33, 0.80, 'with BH', color=50, charthick=3.0, size=1.3, /normal




; extras
; ---------
x=[xmin,xmax]
y=[0.5,0.5]
oplot, x, y, psym=-3, linestyle=2, color=0





; ----------------------------
device, /close

end






;==========================================================================






;===================================
;
;  Figure for Wind Paper (v3)
;
;===================================


;
pro windfig3, junk


if not keyword_set(junk) then begin
	print, " "
	print, " windfig3, junk"
	print, " "
	print, " "
	return
endif

filename='windfig3.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=12.0, newysize= 22.0


;---------------------------

;     ----------------
;     |              |
;     |              |
;     |   SFR        |
;     |              |
;     |              |
;     |              |
;     ----------------
;     |              |
;     |              |
;     |   Outflow    |
;     |    Rate      |
;     |              |
;     |              |
;     ----------------
;     |              |
;     |              |
;     |   Eta_eff    |
;     |              |
;     |              |
;     |              |
;     ----------------



;---------------------------




;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"
xmax = 4.25
;xmax = 3.0
xmin = 0.0


yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
ymax = 300.0
ymin = 7.0e-3



;---------------------------

x0= 0.21
x1= 0.98

y0= 0.09
y1= 0.39
y2= 0.69
y3= 0.99

;---------------------------


!p.position= [x0, y2, x1, y3]


plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        /ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        ;xtitle=xaxistitle, $
        xtickformat='(a1)', $
        ytitle=yaxistitle, $
        /nodata



;---------------------------



process_one_sfr, "vc3vc3e_no", lcolor=0, lstyle= 1, lthick= 2.0, msg=' ', x0= 0.90, y0= 0.90, h=0.7
oplot, [2.25,2.55], [85,85], thick=6.0, color=0, linestyle=1
xyouts, 0.69, 0.95, 'no winds', color=0, charthick=3.0, size=1.3, /normal



process_one_sfr, "sbw/sb10", lcolor=150, lpsym=-5, lthick= 2.0, msg=' ', x0= 0.90, y0= 0.90, h=0.7
oplot, [2.3,2.5], [31,31], thick=6.0, color=150, psym=-5, linestyle=0
xyouts, 0.69, 0.92, 'sb winds', color=150, charthick=3.0, size=1.3, /normal



process_one_sfr, "sbw/sb10BH", lcolor=50, lthick= 2.0, msg=' ', x0= 0.90, y0= 0.90, h=0.7
oplot, [2.25,2.55], [11,11], thick=6.0, color=50, psym=-3, linestyle=0
xyouts, 0.69, 0.89, 'sb winds+BH', color=50, charthick=3.0, size=1.3, /normal



;---------------------------


yaxistitle = "!6Outflow Rate (M!D!9n!6!N Yr!E-1!N)"
ymax = 50.0
ymin = 3.0e-2





;---------------------------

!p.position= [x0, y1, x1, y2]


plot, [1.0],[1.0], psym=-3, /noerase, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        ;xtitle=xaxistitle, $
	xtickformat='(a1)', $
        ytitle=yaxistitle, $
        /nodata



; ----------------------------

read_sbwind_file, "/raid4/tcox/vc3vc3e_no", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
;idx= where(outflowr gt 0.0)
;oplot, time(idx), outflowr(idx), thick=6.0, color=0, psym=-3, linestyle=1
;idx= where(outflowr le 0.0)
;if idx(0) ne -1 then outflowr(idx)= 6.0e-2
;oplot, time, outflowr, thick=6.0, color=0, psym=-3, linestyle=1

oplot, time[2:16], abs(outflowr[2:16]), thick=6.0, color=0, psym=-3, linestyle=1
oplot, time[35:82], abs(outflowr[35:82]), thick=6.0, color=0, psym=-3, linestyle=1
oplot, time[88:102], abs(outflowr[88:102]), thick=6.0, color=0, psym=-3, linestyle=1
oplot, time[105:107], abs(outflowr[105:107]), thick=6.0, color=0, psym=-3, linestyle=1

;oplot, time[16:35], -outflowr[16:35], thick=6.0, color=0, psym=-3, linestyle=0
;oplot, time[82:88], -outflowr[82:88], thick=6.0, color=0, psym=-3, linestyle=0
;oplot, time[102:105], -outflowr[102:105], thick=6.0, color=0, psym=-3, linestyle=0

; ----------------------------

read_sbwind_file, "/raid4/tcox/sbw/sb10", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
idx= where(outflowr gt 0.0)
oplot, time(idx), outflowr(idx), thick=6.0, color=150, psym=-5, linestyle=0
;idx= where(outflowr le 0.0)
;if idx(0) ne -1 then outflowr(idx)= 1.0e-5
;oplot, time, outflowr, thick=6.0, color=150, psym=-5, linestyle=0

; ----------------------------

read_sbwind_file, "/raid4/tcox/sbw/sb10BH", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
idx= where(outflowr gt 0.0)
oplot, time(idx), outflowr(idx), thick=6.0, color=50, psym=-3, linestyle=0
;idx= where(outflowr le 0.0)
;if idx(0) ne -1 then outflowr(idx)= 1.0e-5
;oplot, time, outflowr, thick=6.0, color=50, psym=-3, linestyle=0


; ----------------------------


yaxistitle = "!7g!6!Deff!N"
ymax = 25.0
ymin = 2.5e-2


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, /noerase, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; ----------------------------


read_sbwind_file, "/raid4/tcox/vc3vc3e_no", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
; ---
idx= where(time lt 2.0)
print, 'mean etaeff=', mean(etaeff)
print, 'mean etaeff (t<2)=', mean(etaeff(idx))
print, "mean etaeff (outflow wt'd)=", total(etaeff*outflowr)/total(outflowr)
; ---
;idx=where(etaeff gt 0)
;oplot, time(idx), etaeff(idx), thick=6.0, color=0, psym=-3, linestyle=1
;oplot, [0.3,0.6], [11,11], thick=6.0, color=150, psym=-5, linestyle=0
;xyouts, 0.33, 0.87, 'standard', color=150, charthick=3.0, size=1.3, /normal

oplot, time[2:16], abs(etaeff[2:16]), thick=6.0, color=0, psym=-3, linestyle=1
oplot, time[35:82], abs(etaeff[35:82]), thick=6.0, color=0, psym=-3, linestyle=1
oplot, time[88:102], abs(etaeff[88:102]), thick=6.0, color=0, psym=-3, linestyle=1
oplot, time[105:107], abs(etaeff[105:107]), thick=6.0, color=0, psym=-3, linestyle=1


; ----------------------------


read_sbwind_file, "/raid4/tcox/sbw/sb10", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff 
time= time/0.7
; ---
idx= where(time lt 2.0)
print, 'mean etaeff=', mean(etaeff)
print, 'mean etaeff (t<2)=', mean(etaeff(idx))
print, "mean etaeff (outflow wt'd)=", total(etaeff*outflowr)/total(outflowr)
; ---
idx=where(etaeff gt 0)
oplot, time(idx), etaeff(idx), thick=6.0, color=150, psym=-5, linestyle=0
;oplot, [0.3,0.6], [11,11], thick=6.0, color=150, psym=-5, linestyle=0
;xyouts, 0.33, 0.87, 'standard', color=150, charthick=3.0, size=1.3, /normal


; ----------------------------


read_sbwind_file, "/raid4/tcox/sbw/sb10BH", time, mass_gtx_1, mass_gtx_2, mass_gtx, z_gtx, outflowr, sfr_avg, etaeff  
time= time/0.7
; ---
idx= where(time lt 2.0)
print, 'mean etaeff=', mean(etaeff)
print, 'mean etaeff (t<2)=', mean(etaeff(idx))
print, "mean etaeff (outflow wt'd)=", total(etaeff*outflowr)/total(outflowr)
; ---
idx=where(etaeff gt 0)
oplot, time(idx), etaeff(idx), thick=6.0, color=50, psym=-3, linestyle=0
;oplot, [0.3,0.6], [6,6], thick=6.0, color=50, psym=-5, linestyle=0
;xyouts, 0.33, 0.80, 'with BH', color=50, charthick=3.0, size=1.3, /normal




; extras
; ---------
x=[xmin,xmax]
y=[0.5,0.5]
oplot, x, y, psym=-3, linestyle=2, color=0





; ----------------------------
device, /close

end





