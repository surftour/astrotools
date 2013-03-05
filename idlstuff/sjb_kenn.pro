pro kenn_iso, junk, filename=filename

if not keyword_set(junk) then begin
   print, "  "
   print, "  "
   print, "methods_kenn_iso, junk"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename = 'sjb_kenn_iso.eps'

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename=filename, colortable=4
; should I try and send a different size for this?



;-------------
;  Plot axes 
;-------------

xmax = 3.16e5      ; for direct comparison to kennicutt (or volker for that matter)
;xmin = 0.316
xmin = 0.0316
ymax = 3.16e3
;ymin = 3.16e-6     ; changed from -5
ymin = 3.16e-8     ; changed from -5


logplot = 1
if logplot eq 1 then begin
        xmax = alog10(xmax)
        xmin = alog10(xmin)
        ymax = alog10(ymax)
        ymin = alog10(ymin)
        xaxistitle="Log !4R !D!3gas !N(M!D!9n!3!N pc!E-2!N)"
        yaxistitle="Log !4R !D!3sfr !N(M!D!9n!3!N yr!E-1!N kpc!E-2!N)"
endif


!p.position= [0.20, 0.15, 0.95, 0.95]

plot, [10.0], [10.0], psym=-3, linestyle=0, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


        if logplot eq 1 then ok=oplot_kenn_log('ps',/linear) else ok = oplot_kenn_log('ps')



;-----------------------
;  Put data on graph
;-----------------------2
snapnum= 10
;snapnum= 20

r_kenn= 2.0



;  n=2
; -----
symsize= 1.5
usersym,symsize*[-1,-1,1,1],symsize*[-1,1,1,-1],/fill
frun= "/raid4/tcox/ssft/iso1a/"
ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 150, thick=2.0

;frun= "/data/tcox/Sbc11i4-u38"
;frun= "/data/tcox/Sbc11i4-u17"
;frun= "execute/Sbc11i4-u43"
frun= "/raid4/tcox/ssft/iso1b/"
ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 150, thick=2.0

frun= "/raid4/tcox/ssft/iso/"
ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 150, thick=2.0


;  n=1
; -----
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0   ;,/fill
;frun= "/data/tcox/Sbc11i4-u6"
;frun= "/data/tcox/Sbc11i4-u17"
;frun= "execute/Sbc11i4-u53"
;frun= "/raid1/sbush/antennae/isossft/iso2a/"
frun= "/raid4/tcox/ssft/iso2a/"
ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 100, thick=2.0

;frun= "/data/tcox/Sbc11i4-u16"
;frun= "/data/tcox/Sbc11i4-u17"
;frun= "execute/Sbc11i4-u52"
;frun= "/raid1/sbush/antennae/isossft/iso2b/"
frun= "/raid4/tcox/ssft/iso2b/"
ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 100, thick=2.0




;  n=0
; -----
usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=3.0, /fill
;frun= "/data/tcox/Sbc11i4-u8"
;frun= "execute/Sbc11i4-u50"
;frun= "/raid1/sbush/antennae/isossft/iso3a/"
frun= "/raid4/tcox/ssft/iso3a/"
ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 50, thick=2.0

;frun= "/data/tcox/Sbc11i4-u17"
;frun= "/data/tcox/Sbc11i4-u17"
;frun= "execute/Sbc11i4-u51"
;frun= "/raid1/sbush/antennae/isossft/iso3b/"
frun= "/raid4/tcox/ssft/iso3b/"
ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 50, thick=2.0






;--------------------------------------
device, /close


end






;--------------------------------------
; Actually process the kennicutt point
;--------------------------------------
function process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log

	;ok=fload_snapshot_bh(frun,snapnum, /nopot_in_snap)

        if (fload_snapshot_bh(frun, snapnum, /nopot_in_snap)) then begin
		print, "PROBLEM: opening file"
		return, 0
	endif

	rxy = fload_gas_xyz('rxy',center=[0,0,0])

	sfr = fload_gas_sfr(1) & print, "total sfr= ", total(sfr)
	; trap for bad sfr's (NaN's and negative numbers)
        idx=where(finite(sfr) eq 0, n_i)
        if idx(0) ne -1 then sfr(idx)= 0.0
        idx=where(sfr lt 0)
        if idx(0) ne -1 then sfr(idx)= 0.0

	gas_mass= fload_gas_mass(1)
	gasmfs_mass= fload_gas_mfs(1)

	idx_within_r_kenn = where(rxy LE r_kenn)
        if idx_within_r_kenn(0) ne -1 then begin
		sfr_within_r_kenn = sfr(idx_within_r_kenn)
		mass_within_r_kenn = gas_mass(idx_within_r_kenn)-gasmfs_mass(idx_within_r_kenn)

		sd_kpc2 = !PI * r_kenn * r_kenn
		sd_pc2 = sd_kpc2 * 1e6
		gas_sd = total(mass_within_r_kenn)*(1e10)/sd_pc2    ; units Msolar/pc2
		sfr_sd = total(sfr_within_r_kenn)/sd_kpc2           ;   "   Msolar/kpc2 

		gas_sd_log= [alog10(gas_sd)]
		if sfr_sd gt 0 then begin
		   sfr_sd_log= [alog10(sfr_sd)]
		endif else begin
		   sfr_sd_log= [-6]
		   return, 1
		endelse
	endif else begin
		gas_sd_log= [-6]
		sfr_sd_log= [-6]
		return, 1
	endelse


	print, "gas_sd, sfr_sd = ", gas_sd_log, sfr_sd_log

	return, 0

end

