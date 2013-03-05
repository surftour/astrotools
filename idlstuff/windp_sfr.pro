

;-------------------------------------------------------------------


pro process_one_sfr, frun, lcolor=lcolor, lstyle=lstyle, lthick=lthick, $
				ctab=ctab, msg=msg, y0=y0, x0=x0, h=h


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
	if keyword_set(h) then sfrtime = sfrtime / h
	;if i eq 1 then sfrtime = sfrtime/(0.7)    ; did this for comparing h and noh


	    lthick= 4.0 * lthick

	    if keyword_set(cumulative) then begin
		oplot, sfrtime, sfrmfs, psym=-3, color= lcolor, linestyle= lstyle, thick= lthick
	    endif else begin
		if keyword_set(gasmass) then begin
			oplot, sfrtime, sfrgasmass, psym=-3, color= lcolor, linestyle= lstyle, thick= lthick
		endif else begin
			oplot, sfrtime, sfrsfr, psym=-3, color= lcolor, linestyle= lstyle, thick= lthick
		endelse
	    endelse



; cooling
;xyouts, 0.7, 0.86, "Cooling", /normal, charthick=1, size=1.33, color=0
;oplot, [1.65,2.2],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.7, 0.80, "Explicit", /normal, charthick=1, size=1.33, color=50
;xyouts, 0.7, 0.75, "mzero.cie", /normal, charthick=1, size=1.33, color=100
;xyouts, 0.7, 0.70, "m-00.cie", /normal, charthick=1, size=1.33, color=150

; sfr law
;;xyouts, 0.7, 0.86, "SF Law", /normal, charthick=1, size=1.33, color=0
;xyouts, 0.74, 0.86, '!7s!3!N!D!20E!3!N', /normal, charthick=1, size=1.33, color=0
;oplot, [1.65,2.2],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.7, 0.80, "constant", /normal, charthick=1, size=1.33, color=50
;xyouts, 0.7, 0.75, "dynamical", /normal, charthick=1, size=1.33, color=100
;xyouts, 0.7, 0.70, '1/!7q!3!N', /normal, charthick=1, size=1.33, color=150

; resolution
;xyouts, 0.7, 0.86, "Resolution", /normal, charthick=1, size=1.33, color=0
;oplot, [1.65,2.2],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.7, 0.80, '4', /normal, charthick=1, size=1.33, color=50
;xyouts, 0.7, 0.76, '2', /normal, charthick=1, size=1.33, color=75
;xyouts, 0.7, 0.72, '1', /normal, charthick=1, size=1.33, color=100
;xyouts, 0.7, 0.68, '1/2', /normal, charthick=1, size=1.33, color=125
;xyouts, 0.7, 0.64, '1/4', /normal, charthick=1, size=1.33, color=150
;xyouts, 0.7, 0.60, '1/8', /normal, charthick=1, size=1.33, color=175

; feedback - within on n
;xyouts, 0.65, 0.86, 'n=2 Feedback', /normal, charthick=1, size=1.33, color=0
;;xyouts, 0.65, 0.86, '!7b!3 for n=2', /normal, charthick=1, size=1.33, color=0
;oplot, [1.50,2.35],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.65, 0.80, 'high', /normal, charthick=1, size=1.33, color=50
;xyouts, 0.65, 0.75, 'medium high', /normal, charthick=1, size=1.33, color=87
;xyouts, 0.65, 0.70, 'medium low', /normal, charthick=1, size=1.33, color=124
;xyouts, 0.65, 0.65, 'low', /normal, charthick=1, size=1.33, color=161

; alternate feedback n's
;xyouts, 0.55, 0.86, 'Alternate Feedback n', /normal, charthick=1, size=1.33, color=0
;oplot, [1.20,2.35],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.6, 0.80, 'n=1 high fb', /normal, charthick=1, size=1.33, color=50
;xyouts, 0.6, 0.76, 'n=1 medium fb', /normal, charthick=1, size=1.33, color=75
;xyouts, 0.6, 0.72, 'n=1 low fb', /normal, charthick=1, size=1.33, color=100
;xyouts, 0.6, 0.68, 'n=0 high fb', /normal, charthick=1, size=1.33, color=125
;xyouts, 0.6, 0.64, 'n=0 medium fb', /normal, charthick=1, size=1.33, color=150
;xyouts, 0.6, 0.60, 'n=0 low fb', /normal, charthick=1, size=1.33, color=175

; mh comparison
;xyouts, 0.6, 0.85, "Isothermal Gas", /normal, charthick=1, size=1.33, color=50
;xyouts, 0.6, 0.80, "Full Gas Treatment", /normal, charthick=1, size=1.33, color=150


; springel comparison
;xyouts, 0.3, 0.80, "Conservative Entropy", /normal, charthick=1, size=1.33, color=50
;xyouts, 0.3, 0.75, "Arithmetic", /normal, charthick=1, size=1.33, color=100
;xyouts, 0.3, 0.70, "Geometric", /normal, charthick=1, size=1.33, color=150


; minor paper - iso figure
; ------------------------------
;xyouts, 0.80, 0.84, 'G3', /normal, charthick=3, size=1.33, color=zcolor_orig+deltacolor+deltacolor+deltacolor
;xyouts, 0.80, 0.38, 'G2', /normal, charthick=3, size=1.33, color=zcolor_orig+deltacolor+deltacolor
;xyouts, 0.80, 0.22, 'G1', /normal, charthick=3, size=1.33, color=zcolor_orig+deltacolor
;xyouts, 0.61, 0.36, 'G0', /normal, charthick=3, size=1.33, color=zcolor_orig
;toppos= 0.32
;botpos= 0.01
;xarrow= 0.60
;yarrow= 0.71
;arrow, xarrow, toppos, yarrow, botpos, /data, COLOR=zcolor_orig, THICK=2.0, hthick=3.0



if not keyword_set(y0) then y0= 0.92
if not keyword_set(x0) then x0= 0.75

if keyword_set(msg) then begin
    xyouts, x0, y0, msg, /normal, charthick=3.0, size=1.33, color= lcolor
endif else begin
    xyouts, x0, y0, frun, /normal, charthick=3.0, size=1.33, color= lcolor
endelse





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






pro sfr_winds, junk, $
		filename=filename, $
		cumulative=cumulative, $
		gasmass=gasmass, $
		h=h


if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_4, junk, filename=filename, /h"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=18.0, newysize=16.0
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
if keyword_set(cumulative) then yaxistitle="!6Mass (10!e10!n h!E-1!NM!D!9n!6!N)"
if keyword_set(gasmass) then yaxistitle="!6Gas Mass (10!e10!n h!E-1!NM!D!9n!6!N)"

;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

;xmax = 14.0
;xmax = 12.0
;xmax = 10.0
;xmax = 7.5
;xmax = 6.0
;xmax = 5.0
xmax = 4.25
;xmax = 4.0
;xmax = 3.0
;xmax = 2.8
;xmax = 2.0
;xmax = 1.5
;xmax = 1.3
;xmax = 1.0
;xmax = 0.5
xmin = 0

;ymax = 8000
;ymax = 2500
;ymax = 1750
;ymax = 1500
;ymax = 400
ymax = 300
;ymax = 250.0
;ymax = 180
;ymax = 120
;ymax = 100
;ymax = 90
;ymax = 80
;ymax = 60
;ymax = 50.0
;ymax = 20.0
;ymax = 16.0
;ymax = 15
;ymax = 13
;ymax = 12
;ymax = 10.0
;ymax = 8.0
;ymax = 6.0
;ymax = 5.0
;ymax = 3.5
;ymax = 2.0
;ymax = 1.5
;ymax = 1.3
;ymax = 1.2
;ymax = 0.75
;ymax = 0.6
;ymax = 0.5
;ymin = 0 
;ymin= 0.01
ymin= 0.002
;ymin= 0.001
;ymin= 0.0001
;ymin= 0.00001


; physical units
;if keyword_set(h) then begin
;	xaxistitle = "Time (Gyr)"
;	h = fload_cosmology('h')
;	;xmax = xmax / h
;	;xmin = xmin / h
;endif


;-----------------------------------------
;
;   Load the panels
;
;-----------------------------------------
;
;     |---------------------|
;     |          |          |
;     |          |          |
;     |    1     |    2     |
;     |          |          |
;     |          |          |
;     |---------------------|
;     |          |          |
;     |          |          |
;     |    3     |    4     |
;     |          |          |
;     |          |          |
;     |---------------------|
;
;

x0= 0.12
xs= 0.5*(0.98-x0)
x1= x0+xs
x2= x0+xs+xs

y0= 0.10
ys= 0.5*(0.98-y0)
y1= y0+ys
y2= y0+ys+ys



;  1
; ---
!p.position= [x0,y1,x1,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, $
	xtickformat='(a1)', ytitle=yaxistitle, ytickformat='exp_label'

fload_newcolortable, 4
;process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0, h=h
;process_one_sfr, "sbw/sb6", lcolor=220, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
;process_one_sfr, "sbw/sb2", lcolor=180, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
;process_one_sfr, "sbw/sb1", lcolor=140, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
;process_one_sfr, "sbw/sb3", lcolor=100, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
;process_one_sfr, "sbw/sb4", lcolor=60, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
;process_one_sfr, "sbw/sb5", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
;xyouts, x1-0.15, y2-0.06, '!7g!6=0.005', /normal, charthick=3, size=1.33, color=0
process_one_sfr, "vc3vc3e_no", lcolor=230, lthick=2.0, msg=' ', x0= 0.0, y0=0.0, h=h
process_one_sfr, "sbw/sb9", lcolor=30, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb18", lcolor=80, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb19", lcolor=130, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb7", lcolor=180, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
xyouts, x1-0.15, y2-0.06, '!7g!6=0.05', /normal, charthick=3, size=1.33, color=0

fload_newcolortable, 4
oplot, [0.3,0.7], [0.1,0.1], psym=-3, color=230, thick= 8.0
xyouts, 0.9, 0.09, 'no SB winds', /data, charthick=3, size=1.33, color=230

fload_newcolortable, 1
oplot, [0.3,0.7], [0.05,0.05], psym=-3, color=180, thick= 8.0
xyouts, 0.9, 0.04, '105 km s!E-1!N', /data, charthick=3, size=1.2, color=180

oplot, [0.3,0.7], [0.025,0.025], psym=-3, color=130, thick= 8.0
xyouts, 0.9, 0.020, '209 km s!E-1!N', /data, charthick=3, size=1.2, color=130

oplot, [0.3,0.7], [0.011,0.011], psym=-3, color=80, thick= 8.0
xyouts, 0.9, 0.009, '418 km s!E-1!N', /data, charthick=3, size=1.2, color=80

oplot, [0.3,0.7], [0.005,0.005], psym=-3, color=30, thick= 8.0
xyouts, 0.9, 0.0042, '837 km s!E-1!N', /data, charthick=3, size=1.2, color=30


;oplot, [0.3,0.7], [0.05,0.05], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.04, '837 km s!E-1!N', /data, charthick=3, size=1.2, color=200
;
;oplot, [0.3,0.7], [0.025,0.025], psym=-3, color=140, thick= 8.0
;xyouts, 0.9, 0.020, '418 km s!E-1!N', /data, charthick=3, size=1.2, color=140
;
;oplot, [0.3,0.7], [0.011,0.011], psym=-3, color=80, thick= 8.0
;xyouts, 0.9, 0.009, '209 km s!E-1!N', /data, charthick=3, size=1.2, color=80
;
;oplot, [0.3,0.7], [0.005,0.005], psym=-3, color=20, thick= 8.0
;xyouts, 0.9, 0.0042, '105 km s!E-1!N', /data, charthick=3, size=1.2, color=20



;  2
; ---
!p.position= [x1,y1,x2,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

fload_newcolortable, 4
process_one_sfr, "vc3vc3e_no", lcolor= 230, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb10", lcolor=30, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb15", lcolor=80, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb16", lcolor=130, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb17", lcolor=180, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x2-0.15, y2-0.06, '!7g!6=0.5', /normal, charthick=3, size=1.33, color=0



;  3
; ---
!p.position= [x0,y0,x1,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytitle=yaxistitle, ytickformat='exp_label'

fload_newcolortable, 4
process_one_sfr, "vc3vc3e_no", lcolor=230, lthick=2.0, msg=' ', x0= 0.0, y0=0.0, h=h
process_one_sfr, "sbw/sb8", lcolor=30, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb12", lcolor=80, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb13", lcolor=130, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb14", lcolor=180, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x1-0.15, y1-0.06, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0




;  4
; ---
!p.position= [x1,y0,x2,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytickformat='(a1)'

fload_newcolortable, 4
process_one_sfr, "vc3vc3e_no", lcolor=230, lthick=2.0, msg=' ', x0= 0.0, y0=0.0, h=h
process_one_sfr, "sbw/sb11", lcolor=30, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb20", lcolor=80, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb21", lcolor=130, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb22", lcolor=180, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x2-0.15, y1-0.06, '!7g!6=5.0', /normal, charthick=3, size=1.33, color=0




;--------------------------------------
device, /close


end













;======================================================================





pro sfr_winds_BH, junk, $
		filename=filename, $
		cumulative=cumulative, $
		gasmass=gasmass, $
		h=h


if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_4, junk, filename=filename, /h"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=18.0, newysize=16.0
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
if keyword_set(cumulative) then yaxistitle="!6Mass (10!e10!n h!E-1!NM!D!9n!6!N)"
if keyword_set(gasmass) then yaxistitle="!6Gas Mass (10!e10!n h!E-1!NM!D!9n!6!N)"

;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

;xmax = 14.0
;xmax = 12.0
;xmax = 10.0
;xmax = 7.5
;xmax = 6.0
;xmax = 5.0
xmax = 4.25
;xmax = 4.0
;xmax = 3.0
;xmax = 2.8
;xmax = 2.0
;xmax = 1.5
;xmax = 1.3
;xmax = 1.0
;xmax = 0.5
xmin = 0

;ymax = 8000
;ymax = 2500
;ymax = 1750
;ymax = 1500
;ymax = 400
ymax = 300
;ymax = 250.0
;ymax = 180
;ymax = 120
;ymax = 100
;ymax = 90
;ymax = 80
;ymax = 60
;ymax = 50.0
;ymax = 20.0
;ymax = 16.0
;ymax = 15
;ymax = 13
;ymax = 12
;ymax = 10.0
;ymax = 8.0
;ymax = 6.0
;ymax = 5.0
;ymax = 3.5
;ymax = 2.0
;ymax = 1.5
;ymax = 1.3
;ymax = 1.2
;ymax = 0.75
;ymax = 0.6
;ymax = 0.5
;ymin = 0 
;ymin= 0.01
ymin= 0.002
;ymin= 0.001
;ymin= 0.0001
;ymin= 0.00001


; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;-----------------------------------------
;
;   Load the panels
;
;-----------------------------------------
;
;     |---------------------|
;     |          |          |
;     |          |          |
;     |    1     |    2     |
;     |          |          |
;     |          |          |
;     |---------------------|
;     |          |          |
;     |          |          |
;     |    3     |    4     |
;     |          |          |
;     |          |          |
;     |---------------------|
;
;

x0= 0.12
xs= 0.5*(0.98-x0)
x1= x0+xs
x2= x0+xs+xs

y0= 0.10
ys= 0.5*(0.98-y0)
y1= y0+ys
y2= y0+ys+ys



;  1
; ---
!p.position= [x0,y1,x1,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, $
	xtickformat='(a1)', ytitle=yaxistitle, ytickformat='exp_label'

fload_newcolortable, 4
process_one_sfr, "vc3vc3e_2", lcolor= 150, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x1-0.20, y2-0.06, 'no SB winds', /normal, charthick=3, size=1.33, color=0

fload_newcolortable, 4
oplot, [0.3,0.7], [0.1,0.1], psym=-3, color=200, thick= 8.0
xyouts, 0.9, 0.09, 'no BH', /data, charthick=3, size=1.33, color=200

oplot, [0.3,0.7], [0.05,0.05], psym=-3, color=150, thick= 8.0
xyouts, 0.9, 0.04, 'BH', /data, charthick=3, size=1.2, color=150



;  2
; ---
!p.position= [x1,y1,x2,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

fload_newcolortable, 4
;process_one_sfr, "vc3vc3e_no", lcolor= 200, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb10", lcolor=200, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb10BH", lcolor=20, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x2-0.15, y2-0.06, '!7g!6=0.5', /normal, charthick=3, size=1.33, color=0
xyouts, x2-0.15, y2-0.11, '!8v!6!DW!N=837', /normal, charthick=3, size=1.33, color=0

;xyouts, 0.20, 0.37, 'WindSpeed', /normal, charthick=3, size=1.33, color=0
;xyouts, 0.22, 0.33, '(km/s)', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "sbw/sb13", lcolor=200, ctab= 1, lthick= 2.0, msg='209', x0= 0.20, y0= 0.28
;process_one_sfr, "sbw/sb13BH", lcolor=20, ctab= 1, lthick= 2.0, msg='209 (BH)', x0= 0.20, y0= 0.24
;xyouts, 0.20, 0.18, 'WindEfficiency=2.0', /normal, charthick=3, size=1.33, color=0

;fload_newcolortable, 4
;oplot, [0.3,0.7], [0.1,0.1], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.09, 'no SB winds', /data, charthick=3, size=1.2, color=200

fload_newcolortable, 1
oplot, [0.3,0.7], [0.05,0.05], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.04, '837 km s!E-1!N', /data, charthick=3, size=1.2, color=200
xyouts, 0.9, 0.04, 'no BH', /data, charthick=3, size=1.2, color=200

oplot, [0.3,0.7], [0.025,0.025], psym=-3, color=20, thick= 8.0
;xyouts, 0.9, 0.020, '837 km s!E-1!N, with BH', /data, charthick=3, size=1.2, color=20
xyouts, 0.9, 0.020, 'BH', /data, charthick=3, size=1.2, color=20





;  3
; ---
!p.position= [x0,y0,x1,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytitle=yaxistitle, ytickformat='exp_label'

fload_newcolortable, 4
;process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0, h=h
process_one_sfr, "sbw/sb8", lcolor=200, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb8BH", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x1-0.15, y1-0.06, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0
xyouts, x1-0.15, y1-0.11, '!8v!6!DW!N=837', /normal, charthick=3, size=1.33, color=0

;fload_newcolortable, 4
;oplot, [0.3,0.7], [0.1,0.1], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.09, 'no SB winds', /data, charthick=3, size=1.2, color=200

fload_newcolortable, 1
oplot, [0.3,0.7], [0.05,0.05], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.04, '837 km s!E-1!N', /data, charthick=3, size=1.2, color=200
xyouts, 0.9, 0.04, 'no BH', /data, charthick=3, size=1.2, color=200

oplot, [0.3,0.7], [0.025,0.025], psym=-3, color=20, thick= 8.0
;xyouts, 0.9, 0.020, '837 km s!E-1!N, with BH', /data, charthick=3, size=1.2, color=20
xyouts, 0.9, 0.020, 'BH', /data, charthick=3, size=1.2, color=20





;  4
; ---
!p.position= [x1,y0,x2,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytickformat='(a1)'

fload_newcolortable, 4
;process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0, h=h
process_one_sfr, "sbw/sb13", lcolor=200, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb13BH", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x2-0.15, y1-0.06, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0
xyouts, x2-0.15, y1-0.11, '!8v!6!DW!N=209', /normal, charthick=3, size=1.33, color=0

;fload_newcolortable, 4
;oplot, [0.3,0.7], [0.1,0.1], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.09, 'no SB winds', /data, charthick=3, size=1.2, color=200

fload_newcolortable, 1
oplot, [0.3,0.7], [0.05,0.05], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.04, '209 km s!E-1!N', /data, charthick=3, size=1.2, color=200
xyouts, 0.9, 0.04, 'no BH', /data, charthick=3, size=1.2, color=200

oplot, [0.3,0.7], [0.025,0.025], psym=-3, color=20, thick= 8.0
;xyouts, 0.9, 0.020, '209 km s!E-1!N, with BH', /data, charthick=3, size=1.2, color=20
xyouts, 0.9, 0.020, 'BH', /data, charthick=3, size=1.2, color=20





;--------------------------------------
device, /close


end













;======================================================================








pro sfr_bigp, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_bigp, junk "
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=18.0, newysize=16.0
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"

xaxistitle = "!6Time (Gyr)"  & h= 0.7

;xmax = 14.0
;xmax = 12.0
;xmax = 10.0
;xmax = 7.5
;xmax = 6.0
;xmax = 5.0
xmax = 4.25
;xmax = 4.0
;xmax = 3.0
;xmax = 2.8
;xmax = 2.0
;xmax = 1.5
;xmax = 1.3
;xmax = 1.0
;xmax = 0.5
xmin = 0

;ymax = 8000
;ymax = 2500
;ymax = 1750
;ymax = 1500
;ymax = 400
ymax = 300
;ymax = 250.0
;ymax = 180
;ymax = 120
;ymax = 100
;ymax = 90
;ymax = 80
;ymax = 60
;ymax = 50.0
;ymax = 20.0
;ymax = 16.0
;ymax = 15
;ymax = 13
;ymax = 12
;ymax = 10.0
;ymax = 8.0
;ymax = 6.0
;ymax = 5.0
;ymax = 3.5
;ymax = 2.0
;ymax = 1.5
;ymax = 1.3
;ymax = 1.2
;ymax = 0.75
;ymax = 0.6
;ymax = 0.5
;ymin = 0 
;ymin= 0.01
ymin= 0.002
;ymin= 0.001
;ymin= 0.0001
;ymin= 0.00001


; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;-----------------------------------------
;
;   Load the panels
;
;-----------------------------------------
;
;     |---------------------|
;     |          |          |
;     |          |          |
;     |    1     |    2     |
;     |          |          |
;     |          |          |
;     |---------------------|
;     |          |          |
;     |          |          |
;     |    3     |    4     |
;     |          |          |
;     |          |          |
;     |---------------------|
;
;

x0= 0.12
xs= 0.5*(0.98-x0)
x1= x0+xs
x2= x0+xs+xs

y0= 0.10
ys= 0.5*(0.98-y0)
y1= y0+ys
y2= y0+ys+ys



;  1
; ---
!p.position= [x0,y1,x1,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, $
	xtickformat='(a1)', ytitle=yaxistitle, ytickformat='exp_label'

fload_newcolortable, 4
process_one_sfr, "vc3vc3e_2", lcolor= 150, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x1-0.20, y2-0.06, 'no SB winds', /normal, charthick=3, size=1.33, color=0

fload_newcolortable, 4
oplot, [0.3,0.7], [0.1,0.1], psym=-3, color=200, thick= 8.0
xyouts, 0.9, 0.09, 'no BH', /data, charthick=3, size=1.33, color=200

oplot, [0.3,0.7], [0.05,0.05], psym=-3, color=150, thick= 8.0
xyouts, 0.9, 0.04, 'BH', /data, charthick=3, size=1.2, color=150



;  2
; ---
!p.position= [x1,y1,x2,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

fload_newcolortable, 4
;process_one_sfr, "vc3vc3e_no", lcolor= 200, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb10", lcolor=200, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb10BH", lcolor=20, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x2-0.15, y2-0.06, '!7g!6=0.5', /normal, charthick=3, size=1.33, color=0
xyouts, x2-0.15, y2-0.11, '!8v!6!DW!N=837', /normal, charthick=3, size=1.33, color=0

;xyouts, 0.20, 0.37, 'WindSpeed', /normal, charthick=3, size=1.33, color=0
;xyouts, 0.22, 0.33, '(km/s)', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "sbw/sb13", lcolor=200, ctab= 1, lthick= 2.0, msg='209', x0= 0.20, y0= 0.28
;process_one_sfr, "sbw/sb13BH", lcolor=20, ctab= 1, lthick= 2.0, msg='209 (BH)', x0= 0.20, y0= 0.24
;xyouts, 0.20, 0.18, 'WindEfficiency=2.0', /normal, charthick=3, size=1.33, color=0

;fload_newcolortable, 4
;oplot, [0.3,0.7], [0.1,0.1], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.09, 'no SB winds', /data, charthick=3, size=1.2, color=200

fload_newcolortable, 1
oplot, [0.3,0.7], [0.05,0.05], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.04, '837 km s!E-1!N', /data, charthick=3, size=1.2, color=200
xyouts, 0.9, 0.04, 'no BH', /data, charthick=3, size=1.2, color=200

oplot, [0.3,0.7], [0.025,0.025], psym=-3, color=20, thick= 8.0
;xyouts, 0.9, 0.020, '837 km s!E-1!N, with BH', /data, charthick=3, size=1.2, color=20
xyouts, 0.9, 0.020, 'BH', /data, charthick=3, size=1.2, color=20





;  3
; ---
!p.position= [x0,y0,x1,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytitle=yaxistitle, ytickformat='exp_label'

fload_newcolortable, 4
;process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0, h=h
process_one_sfr, "sbw/sb8", lcolor=200, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb8BH", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x1-0.15, y1-0.06, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0
xyouts, x1-0.15, y1-0.11, '!8v!6!DW!N=837', /normal, charthick=3, size=1.33, color=0

;fload_newcolortable, 4
;oplot, [0.3,0.7], [0.1,0.1], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.09, 'no SB winds', /data, charthick=3, size=1.2, color=200

fload_newcolortable, 1
oplot, [0.3,0.7], [0.05,0.05], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.04, '837 km s!E-1!N', /data, charthick=3, size=1.2, color=200
xyouts, 0.9, 0.04, 'no BH', /data, charthick=3, size=1.2, color=200

oplot, [0.3,0.7], [0.025,0.025], psym=-3, color=20, thick= 8.0
;xyouts, 0.9, 0.020, '837 km s!E-1!N, with BH', /data, charthick=3, size=1.2, color=20
xyouts, 0.9, 0.020, 'BH', /data, charthick=3, size=1.2, color=20





;  4
; ---
!p.position= [x1,y0,x2,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytickformat='(a1)'

fload_newcolortable, 4
;process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0, h=h
process_one_sfr, "sbw/sb13", lcolor=200, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h
process_one_sfr, "sbw/sb13BH", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0, h=h

xyouts, x2-0.15, y1-0.06, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0
xyouts, x2-0.15, y1-0.11, '!8v!6!DW!N=209', /normal, charthick=3, size=1.33, color=0

;fload_newcolortable, 4
;oplot, [0.3,0.7], [0.1,0.1], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.09, 'no SB winds', /data, charthick=3, size=1.2, color=200

fload_newcolortable, 1
oplot, [0.3,0.7], [0.05,0.05], psym=-3, color=200, thick= 8.0
;xyouts, 0.9, 0.04, '209 km s!E-1!N', /data, charthick=3, size=1.2, color=200
xyouts, 0.9, 0.04, 'no BH', /data, charthick=3, size=1.2, color=200

oplot, [0.3,0.7], [0.025,0.025], psym=-3, color=20, thick= 8.0
;xyouts, 0.9, 0.020, '209 km s!E-1!N, with BH', /data, charthick=3, size=1.2, color=20
xyouts, 0.9, 0.020, 'BH', /data, charthick=3, size=1.2, color=20





;--------------------------------------
device, /close


end













;======================================================================





