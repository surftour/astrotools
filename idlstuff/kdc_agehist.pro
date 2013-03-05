;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Histograms of Stellar Ages
;     ------------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------


pro age_hist, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "age_hist, junk"
   print, "  "
   return
endif

filename='agehist.eps'


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; radius, linear scale
xaxistitle= "!6Time (!8h!6!E-1!NGyr)"
xmax = 3.0     ; also set these in the process_hist procedure
xmin = 0.0

; number (histogram)
;yaxistitle= "!6V!Dr!N (km s!E-1!N)"
yaxistitle= ' '
ymax = 1.05
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.05
y0= 0.15
x1= 0.95
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
	;/ylog, $
	;/xlog, $
	color= 0, $
	xrange=[xmin,xmax], $
	yrange=[ymin,ymax], $
	xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
	xtitle=xaxistitle, $
	ytickformat='(a1)', $
	;ytitle=yaxistitle, $
	/nodata   ;, /noerase


; -----------------------------------------------

;frun= "/raid4/tcox/vc3vc3i" & snapnum= 30 & mannorm= 7623.0 & kdcnormfac= 10.0
frun= "/raid4/tcox/vc3vc3m" & snapnum= 30 & mannorm= 7730.0 & kdcnormfac= 20.0

ok=fload_snapshot_bh(frun,snapnum)

as_age= fload_allstars_age(1)
as_ids= fload_allstars_ids(1)

; the whole shebang
oplotit= 50
temp= process_histogram(as_age, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=mannorm)
print, min(temp), max(temp)
xyouts, 0.70, 0.80, "all stars", size=1.5, color=oplotit, /normal, charthick=3.0



; -----------------------------------------------

; the whole shebang
oplotit= 150
idx= intarr(n_elements(as_ids))
lstid= -1
duplicates= 0
idlist_fmfile= fload_id_list(frun+'/id_list_test.txt')
print, "id's in file= ",n_elements(idlist_fmfile)
sidx=sort(idlist_fmfile)
idlist_fmfile= idlist_fmfile(sidx)
for i=0,n_elements(idlist_fmfile)-1 do begin
    if idlist_fmfile[i] eq lstid then begin
	duplicates= duplicates+1
	;print, lstid, idlist_fmfile[i]
    endif
    inlst= where(as_ids eq idlist_fmfile[i])
    if inlst(0) ne -1 then begin
	idx(inlst)= 1
    endif else begin
	;print, "couldn't find id=",idlist_fmfile[i]
    endelse
    lstid= idlist_fmfile[i]
endfor

midx= where(idx eq 1)
print, "matching id's= ",n_elements(midx)
print, "duplicate id's= ",duplicates

; now the histogram
if midx(0) ne -1 then begin
	idages= as_age(midx)
	temp= process_histogram(idages, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=mannorm/kdcnormfac)
	print, min(temp), max(temp)
	xyouts, 0.70, 0.75, "KDC stars", size=1.5, color=oplotit, /normal, charthick=3.0
endif

; -----------------------------------------------


; add star formation rate on top of this

;   Get SFR rate from txt - for each file
;-------------------------------------------
open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass


sfrsfr=sfrsfr/max(sfrsfr)
oplot, sfrtime, sfrsfr, psym=-3, color= 0, linestyle= 0, thick= 3.0



; -----------------------------------------------



; print extras
; -------------

; done
; -----
device, /close


end



;========================================================================================


