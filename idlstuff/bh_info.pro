pro bh_one, frun, $
                xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, $
                lcolor=lcolor, lstyle=lstyle, lthick=lthick, $
                ctab=ctab, msg=msg, y0=y0, x0=x0, lpsym=lpsym, $
		filename=filename, suppressfid=suppressfid, $
		fillinregion=fillinregion, $
		h=h, msg=msg


if not keyword_set(frun) then begin
   print, "  "
   print, "bh_info, frun, filename=filename, /suppressfid,"
   print, "             /fillinregion, /h, msg=msg"
   print, "  "
   print, "  "
   return
endif



if keyword_set(filename) then begin
	if strmid(filename,strlen(filename)-3,3) eq "eps" then fname_base= strmid(filename,strlen(filename)-4,4)
	bhmass_filename=fname_base+'_bhmass.eps'
	bhacc_filename=fname_base+'_bhacc.eps'
endif else begin
	bhmass_filename='bhmass.eps'
	bhacc_filename='bhacc.eps'
endif



   
if not keyword_set(xmax) then xmax = 2.4
if not keyword_set(xmin) then xmin = 0

if not keyword_set(ymax) then ymax = 250.0
if not keyword_set(ymin) then ymin = 0.1


if not keyword_set(lcolor) then lcolor= 150
if not keyword_set(lthick) then lthick= 2.0
if not keyword_set(lstyle) then lstyle= 0

;if not keyword_set(msg) then msg=' '

if not keyword_set(x0) then x0= 0.70
if not keyword_set(y0) then y0= 0.85




;-------------------------------------
;   Get BH Info from fid.bh file
;-------------------------------------


open_blackholes_file, fruns[i], bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass, bh_mdot_edd


bhtime= bhtime / 0.7
bh_mass= 1.0e+10*bh_mass / 0.7
bh_totalmass= 1.0e+10*bh_totalmass / 0.7


ok=fload_snapshot_bh(frun,0,/nopot_in_snap)

n_bh= fload_npart(5)
bhid= fload_blackhole_id(1)


;
; bh #1
;
read_bh_file, frun, bhid1, time1, bhm1, pm1, mdot1, mdotedd1, bhfile=frun+"/bh_"+bhidlbl1+".txt", mergertime=mergertime

;
; bh #2
;
read_bh_file, frun, bhid2, time2, bhm2, pm2, mdot2, mdotedd2, bhfile=frun+"/bh_"+bhidlbl1+".txt", mergertime=mergertime






;--------------------------------
;  PLOT 1
;
;  BH Mass
;
;--------------------------------

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=bhmass_filename, colortable=4


; physical units
;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
;if keyword_set(h) then xaxistitle = "!6Time (Gyr)"
xaxistitle = "!6Time (Gyr)"



yaxistitle="!6Black Hole Mass (!8h!6!E-1!N M!D!9n!6!N)"
;ymax= 7.0e+10
;ymax= 7.0e+9
;ymax= 3.0e+9
;ymax= 3.0e+8
ymax= 1.2e+8

;ymin= 2.0e+7
ymin= 7.0e+6
;ymin= 6.0e+5
;ymin= 2.0e+5
;ymin= 7.0e+4
;ymin= 1.0e+3
;ymin= 1.0e+0



!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $ 
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata         



;---------------------------

oplot, bhtime, bh_mass, psym=-3, color= zcolor, linestyle= lstyle, thick= 4.0



; ---------------
if not keyword_set(msg) then msg= fload_fid(1)
xyouts, 0.75, 0.90, msg, /normal, charthick=3.0, size=1.5, color= 0


; the end
;
device, /close






;--------------------------------
;  PLOT 2
;
;  BH Accretion Rate
;
;--------------------------------

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=bhacc_filename, colortable=4


;yaxistitle="!6Accretion Rate (M!D!9n!6!N Yr!E-1!N)"
yaxistitle="!6Accretion Rate (Eddington)"



; acc. rate (msolar/yr)
;ymax= 5.0e+2
;ymax= 1.0e+1
;ymax= 2.0e+0
;ymax= 1.0e+0
;ymin= 1.0e-4
;ymin= 1.0e-5
;ymin= 1.0e-8
;ymin= 1.0e-9


;ymax= 10
;ymax= 7
ymax= 4
;ymin= 1.0e-3
ymin= 1.0e-4
;ymin= 1.0e-5


!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $ 
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata         



;---------------------------

;oplot, bhtime, bh_mdot_sunyr, thick=4.0, psym=-3, color= zcolor, linestyle= lstyle

oplot, bhtime, bh_mdot_edd, thick=4.0, psym=-3, color= zcolor, linestyle= lstyle


; ---------------
if not keyword_set(msg) then msg= fload_fid(1)
xyouts, 0.75, 0.90, msg, /normal, charthick=3.0, size=1.5, color= 0


; the end
;
device, /close










end








; ================================================================================












