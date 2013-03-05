

;-------------------------------------------------------------------


pro process_one_sfr, frun, lcolor=lcolor, lstyle=lstyle, lthick=lthick, $
				ctab=ctab, msg=msg, y0=y0, x0=x0, h=h, $
				gasmass=gasmass, normalize=normalize


	if not keyword_set(lthick) then lthick= 1.0
	if not keyword_set(lstyle) then lstyle= 0
	if not keyword_set(lcolor) then lcolor= 0


	if keyword_set(ctab) then begin
                        loadct, ctab
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0
	endif

	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass


	; physical units
	if keyword_set(h) then begin
		sfrtime = sfrtime / h
		sfrgasmass = sfrgasmass /h
	endif

	if keyword_set(normalize) then sfrgasmass= sfrgasmass / sfrgasmass(0)
	if lcolor eq 150 then lpsym=-5 else lpsym=-3

	    lthick= 4.0 * lthick

	    if keyword_set(cumulative) then begin
		oplot, sfrtime, sfrmfs, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick
	    endif else begin
		if keyword_set(gasmass) then begin
			oplot, sfrtime, sfrgasmass, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick
		endif else begin
			oplot, sfrtime, sfrsfr, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick
		endelse
	    endelse



if not keyword_set(y0) then y0= 0.92
if not keyword_set(x0) then x0= 0.75

if keyword_set(msg) then begin
    xyouts, x0, y0, msg, /normal, charthick=3.0, size=1.33, color= lcolor
endif else begin
    xyouts, x0, y0, frun, /normal, charthick=3.0, size=1.33, color= lcolor
endelse





end






;======================================================================





pro process_one_wind, frun, lcolor=lcolor, lstyle=lstyle, lthick=lthick, $
				ctab=ctab, msg=msg, y0=y0, x0=x0, h=h, $
				gasmass=gasmass, normalize=normalize


	if not keyword_set(lthick) then lthick= 1.0
	if not keyword_set(lstyle) then lstyle= 0
	if not keyword_set(lcolor) then lcolor= 0


	if keyword_set(ctab) then begin
                        loadct, ctab
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0
	endif

	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass


        ;------------------
        ;  get wind info
        ;------------------
        wfrun= '/raid4/tcox/'+frun
        read_wind_file, wfrun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
        windm= mass_egt


        if keyword_set(h) then begin
		time=time/h
		windm= windm/h
		sfrgasmass= sfrgasmass/h
	endif

	if keyword_set(normalize) then windm= windm / sfrgasmass(0)

	lthick= 4.0 * lthick

	if lcolor eq 150 then lpsym=-5 else lpsym=-3

	oplot, time, windm, color= lcolor, linestyle= lstyle, thick= lthick, psym=lpsym


end






;======================================================================






pro open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass

	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	    ; get sfr data (from file)
	    ; -------------------------
	    loopnum= 0
	    repeat begin
	    ;sfrfile= '/data/tcox/sfr/'+fload_getid(fruns[i])+'.sfr'        ; twopiter
	    ;sfrfile= '/home/tcox/data/sfr/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/raid4/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/raid2/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/data6/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard

		case loopnum of
		  0: datadir='/raid4'
		  1: datadir='/home'
		  2: datadir='/raid2'
		  3: datadir='/data'
		  4: datadir='/data6'
		  5: datadir='/data7'
		  else: break
		endcase

		;sfrfile= datadir+'/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'

		nfrun= fload_getid(frun)
		spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/*.sfr ',result

		sfrfile=strcompress(result[0],/remove_all)



		get_lun, unit
		openr, unit, sfrfile, ERROR=err
		close, unit
		free_lun, unit

		if (err NE 0) then begin
		   ERR= 0
		   foundit= 0
		endif else begin
		   foundit= 1
		endelse

		loopnum= loopnum+1

	    endrep until (foundit eq 1)

	    if foundit eq 0 then begin
		msg= 'Problem: sfr file no found  file=[data]/tcox/'+fload_getid(frun)+'/'+fload_getid(frun)+'.sfr'
	        print, "  "
	        print, msg
	        print, "  "
		return
	    endif

	    ; read to open
	    print, "opening: ",sfrfile
	    sfrdata= read_ascii(sfrfile)
	    sfrtime= sfrdata.field1[0,*]
	    sfrsfr= sfrdata.field1[1,*]
	    sfrmfs= sfrdata.field1[3,*]    ; amount of new stars
	    n_cols= n_elements(sfrdata.field1[*,0])
	    n_rows= n_elements(sfrdata.field1[0,*])
	    finalnewstarmass= sfrmfs[n_rows-1]
	    ; gas mass
	    if n_cols ge 6 then sfrgasmass= sfrdata.field1[6,*] else sfrgasmass=[20.0,20.0-finalnewstarmass]
	    sfrmsg= ''


            t1per= 0
            idx=where(sfrtime ge 1.0)
            if idx(0) ne -1 and n_cols gt 5 then begin
                gm= sfrgasmass[0]-sfrgasmass(idx(0))
                t1per= 100.0*gm/sfrgasmass[0]
            endif

            t4per= 0
            idx=where(sfrtime ge 4.0)
            if idx(0) ne -1 and n_cols gt 5 then begin
                gm= sfrgasmass[0]-sfrgasmass(idx(0))
                t4per= 100.0*gm/sfrgasmass[0]
            endif

            maxsfr= max(sfrsfr)
            idx= where(sfrsfr eq maxsfr)
            maxsfrtime= sfrtime(idx)

            sfraftertmerg= 0.0
            smaftertmerg= 0.0
            taftertmerg= 0.0
            tmerger= 0.0
            if n_elements(tmerg) gt 1 then begin
                idx= where(sfrtime ge tmerg[i]+0.2)
                if idx(0) eq -1 then begin
                        sfraftertmerg= sfrsfr[n_rows-1]
                        smaftertmerg= sfrmfs[n_rows-1]
                        tmerger= sfrtime[n_rows-1]
                        taftertmerg= tmerger
                endif else begin
                        sfraftertmerg= sfrsfr[idx(0)]
                        smaftertmerg= sfrmfs[idx(0)]
                        tmerger= tmerg[i]
                        taftertmerg= sfrtime[idx(0)]
                endelse
            endif

            print, "-------------------------------------"
            print, "maximum sfr= ", maxsfr, " occurs at ", maxsfrtime
            print, "average sfr= ", total(sfrsfr)/n_rows
            print, "  "
            print, "      tmerg= ",tmerger,'  and 200 Myr later= ',taftertmerg
            print, "        sfr= ", sfraftertmerg,'  200 Myr after t_merger'
            print, "         sm= ", smaftertmerg,'  200 Myr after t_merger'
            n_mfs= n_elements(sfrgasmass)
            print, "-------------------------------------"
            print, "original gas mass   =", sfrgasmass[0]
            print, "remnant gas mass    =", sfrgasmass[n_mfs-1]
            print, "gas consumed        =", sfrgasmass[0]-sfrgasmass[n_mfs-1]
            print,"                      ", 100.0*(sfrgasmass[0]-sfrgasmass[n_mfs-1])/sfrgasmass[0],'  %'
            print,"                      ", t1per,' % after 1 Gyr'
            print,"                      ", t4per,' % after 4 Gyr'
            print, "-------------------------------------"


end






;==============================================================================






pro wp, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "wp, junk"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='wp.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=18.0, newysize=16.0


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
yaxistitle="!6Mass (10!e10!n h!E-1!NM!D!9n!6!N)"
yaxistitle="!6Gas Mass (10!e10!n h!E-1!NM!D!9n!6!N)"

xaxistitle = "!6Time (Gyr)"  & h= 0.7

xmax = 4.25
xmin = 0

ymax = 300
ymin= 0.002



;-----------------------------------------
;
;   Load the panels
;
;-----------------------------------------
;
;
;        vc1         vc3        vc5
;
;     |---------------------------------|
;     |          |          |           |
; S   |          |          |           |
; F   |    1     |    2     |    3      |
; R   |          |          |           |
;     |          |          |           |
;     |---------------------------------|
;     |          |          |           |
; M   |          |          |           |
; a   |    4     |    5     |    6      |
; s   |          |          |           |
; s   |          |          |           |
;     |---------------------------------|
;     |          |          |           |
; M   |          |          |           |
; a   |    7     |    8     |    9      |
; s   |          |          |           |
; s   |          |          |           |
;     |---------------------------------|
;
;

x0= 0.12
xs= (1./3.)*(0.98-x0)
x1= x0+xs
x2= x0+xs+xs
x3= x0+xs+xs+xs

y0= 0.10
ys= (1./3.)*(0.98-y0)
y1= y0+ys
y2= y0+ys+ys
y3= y0+ys+ys+ys




; top row
;======================

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"

xaxistitle = "!6Time (Gyr)"  & h= 0.7

xmax = 4.25
xmin = 0

ymax = 700
ymin= 0.002




;  1
; ---
!p.position= [x0,y2,x1,y3]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, $
	xtickformat='(a1)', ytitle=yaxistitle, ytickformat='exp_label'

process_one_sfr, "ds/d1e2_q", lcolor=100, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sb10_mass/d1e", lcolor=50, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sb10_mass/d1e_no", lcolor=150, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
xyouts, x1-0.1, y3-0.06, '10!E11!N M!D!9n!6!N', /normal, charthick=3, size=1.33, color=0



oplot, [0.3,0.7], [0.08,0.08], psym=-3, color=100, thick= 8.0
xyouts, 0.9, 0.07, 'BH', /data, charthick=3, size=1.33, color=100

oplot, [0.3,0.7], [0.025,0.025], psym=-3, color=50, thick= 8.0
xyouts, 0.9, 0.020, 'eta=0.5, !8v!6!DW!N=837 km s!E-1!N', /data, charthick=3, size=1.2, color=50

oplot, [0.3,0.7], [0.01,0.01], psym=-3, color=150, thick= 8.0
;xyouts, 0.9, 0.020, 'eta=0.5, !8v!6!DW!N=837 km s!E-1!N, +BH', /data, charthick=3, size=1.2, color=150
xyouts, 0.9, 0.008, 'sb winds + BH', /data, charthick=3, size=1.2, color=150



;  2
; ---
!p.position= [x1,y2,x2,y3]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

process_one_sfr, "ds/d3e7", lcolor=100, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sb10_mass/d3e", lcolor=50, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sb10_mass/d3e_no", lcolor=150, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x2-0.1, y3-0.06, '10!E12!N M!D!9n!6!N', /normal, charthick=3, size=1.33, color=0



;  3
; ---
!p.position= [x2,y2,x3,y3]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

process_one_sfr, "ds/d5e2_q", lcolor=100, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sb10_mass/d5e", lcolor=50, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sb10_mass/d5e_no", lcolor=150, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x3-0.1, y3-0.06, '10!E13!N M!D!9n!6!N', /normal, charthick=3, size=1.33, color=0





; middle row
;======================

yaxistitle="!6Gas Mass / M!D0!N"

ymax = 1.0
ymin = 0.0



;  4
; ---
!p.position= [x0,y1,x1,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtickformat='(a1)', ytitle=yaxistitle

process_one_sfr, "ds/d1e2_q", lcolor=100, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass, /normalize
process_one_sfr, "sb10_mass/d1e", lcolor=50, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass, /normalize
process_one_sfr, "sb10_mass/d1e_no", lcolor=150, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass, /normalize




;  5
; ---
!p.position= [x1,y1,x2,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

process_one_sfr, "ds/d3e7", lcolor=100, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass, /normalize
process_one_sfr, "sb10_mass/d3e", lcolor=50, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass, /normalize
process_one_sfr, "sb10_mass/d3e_no", lcolor=150, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass, /normalize




;  6
; ---
!p.position= [x2,y1,x3,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

process_one_sfr, "ds/d5e2_q", lcolor=100, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass, /normalize
process_one_sfr, "sb10_mass/d5e", lcolor=50, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass, /normalize
process_one_sfr, "sb10_mass/d5e_no", lcolor=150, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass, /normalize






; bottom row
;======================

yaxistitle="!6Wind Mass / M!D0!N"

ymax = 0.75
ymin= 0.0



;  7
; ---
!p.position= [x0,y0,x1,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytitle=yaxistitle

process_one_wind, "ds/d1e2_q", lcolor=100, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /normalize
process_one_wind, "sb10_mass/d1e", lcolor=50, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /normalize
process_one_wind, "sb10_mass/d1e_no", lcolor=150, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /normalize




;  8
; ---
!p.position= [x1,y0,x2,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytickformat='(a1)'

process_one_wind, "ds/d3e7", lcolor=100, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /normalize
process_one_wind, "sb10_mass/d3e", lcolor=50, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /normalize
process_one_wind, "sb10_mass/d3e_no", lcolor=150, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /normalize




;  9
; ---
!p.position= [x2,y0,x3,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytickformat='(a1)'

process_one_wind, "ds/d5e2_q", lcolor=100, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /normalize
process_one_wind, "sb10_mass/d5e", lcolor=50, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /normalize
process_one_wind, "sb10_mass/d5e_no", lcolor=150, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h, /normalize









;--------------------------------------
device, /close


end













;======================================================================





pro wp2, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "wp2, junk"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='wp2.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=13.0, newysize=15.0


;--------------------------------------
;  Print the Shit
;--------------------------------------

xaxistitle = "!6Time (Gyr)"  & h= 0.7
xmax = 4.25
xmin = 0


;---------------------------

;
;  BIG - this is how we fix which
;        plot to do
;
do_sbwind= 1
;do_sbwind= 0

if do_sbwind eq 1 then do_nowind= 0 else do_nowind= 1




;-----------------------------------------
;
;   Load the panels
;
;-----------------------------------------
;
;
;
;     |-----------
;     |          |
;     |          |
;     |    1     |
;     |          |
;     |          |
;     |-----------
;     |          |
;     |          |
;     |    2     |
;     |          |
;     |          |
;     |-----------
;
;

x0= 0.18
x1= 0.98

y0= 0.10
ys= (1./2.)*(0.98-y0)
y1= y0+ys
y2= y0+ys+ys




; top row
;======================

xaxistitle = "!6Time (Gyr)"  & h= 0.7
xmax = 4.25
xmin = 0

yaxistitle="!6Gas Mass (10!E10!N M!D!9n!6!N)"
ymin = 0.0
ymax = 4.8



; ---
!p.position= [x0,y1,x1,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtickformat='(a1)', ytitle=yaxistitle

;process_one_sfr, "ds/d3e7", lcolor=100, lthick= 2.0, lstyle= 0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass
;process_one_sfr, "sb10_mass/d3e", lcolor=50, lthick= 2.0, lstyle= 1, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass
;process_one_sfr, "sb10_mass/d3e_no", lcolor=150, lthick= 2.0, lstyle= 2, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass


if do_sbwind eq 1 then begin
	process_one_sfr, "sbw/sb10", lcolor=150, lthick= 2.0, lstyle= 0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass
	process_one_sfr, "sbw/sb10BH", lcolor=50, lthick= 2.0, lstyle= 0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass
endif

if do_nowind eq 1 then begin
	process_one_sfr, "vc3vc3e_no", lcolor=150, lthick= 2.0, lstyle= 0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass
	process_one_sfr, "vc3vc3e_2", lcolor=50, lthick= 2.0, lstyle= 0, msg=' ', x0= 0.0, y0= 0.0, h=h, /gasmass
endif




; bottom row
;======================

yaxistitle="!6Wind Mass (10!E10!N M!D!9n!6!N)"
ymax = 2.8
ymin= 0.0



; ---
!p.position= [x0,y0,x1,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
	;yticks= 5, $
	ytickname=['0',' ','1',' ','2',' '], $
        xtitle=xaxistitle, ytitle=yaxistitle

;process_one_wind, "ds/d3e7", lcolor=100, lthick= 2.0, lstyle= 0, msg=' ', x0= 0.0, y0= 0.0, h=h
;process_one_wind, "sb10_mass/d3e", lcolor=50, lthick= 2.0, lstyle= 1, msg=' ', x0= 0.0, y0= 0.0, h=h
;process_one_wind, "sb10_mass/d3e_no", lcolor=150, lthick= 2.0, lstyle= 2, msg=' ', x0= 0.0, y0= 0.0, h=h

if do_sbwind eq 1 then begin
	process_one_wind, "sbw/sb10", lcolor=150, lthick= 2.0, lstyle= 0, msg=' ', x0= 0.0, y0= 0.0, h=h
	process_one_wind, "sbw/sb10BH", lcolor=50, lthick= 2.0, lstyle= 0, msg=' ', x0= 0.0, y0= 0.0, h=h
endif

if do_nowind eq 1 then begin
	process_one_wind, "vc3vc3e_no", lcolor=150, lthick= 2.0, lstyle= 0, msg=' ', x0= 0.0, y0= 0.0, h=h
	process_one_wind, "vc3vc3e_2", lcolor=50, lthick= 2.0, lstyle= 0, msg=' ', x0= 0.0, y0= 0.0, h=h
endif


; labels and such
;

if do_nowind eq 1 then begin
        xyouts, 0.62, 0.83, 'no SB winds', color=0, charthick=3.0, size=1.8, /normal
endif

if do_sbwind eq 1 then begin
        xyouts, 0.67, 0.83, 'SB winds', color=0, charthick=3.0, size=1.8, /normal
endif


oplot, [0.4,0.6], [2.4,2.4], thick=6.0, color=150, psym=-5, linestyle=0
xyouts, 0.33, 0.47, 'without BH', color=150, charthick=3.0, size=1.3, /normal
oplot, [0.35,0.65], [2.1,2.1], thick=6.0, color=50, psym=-3, linestyle=0
xyouts, 0.33, 0.43, 'BH', color=50, charthick=3.0, size=1.3, /normal


;oplot, [0.3,0.7], [3.0,3.0], psym=-3, color=100, thick= 8.0, linestyle= 0
;xyouts, 0.9, 2.95, 'BH', /data, charthick=3, size=1.33, color=100
;
;oplot, [0.3,0.7], [2.4,2.4], psym=-3, color=50, thick= 8.0, linestyle= 1
;xyouts, 0.9, 2.35, 'sb winds, no BH', /data, charthick=3, size=1.2, color=50
;
;oplot, [0.3,0.7], [1.8,1.8], psym=-3, color=150, thick= 8.0, linestyle= 2
;xyouts, 0.9, 1.75, 'sb winds + BH', /data, charthick=3, size=1.2, color=150


;--------------------------------------
device, /close


end













;======================================================================







