;-----------------------------------------------------------
;
;
;
;
;
;-----------------------------------------------------------
pro kenn_iso, junk, filename=filename

if not keyword_set(junk) then begin
   print, "  "
   print, "  "
   print, "kenn_iso, junk"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename = 'kenn_iso.eps'

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

        if logplot eq 1 then begin
		ok=oplot_kenn_log('ps',/linear)

		h = 0.70     ; 1.0 means in units of h, h-2 etc.
		; h = 1.0

		x=100000*findgen(10)+0.1
		y=(2.5e-4)*((h/0.75)^0.6)*(x^(1.5))

		x= alog10(x)
		y= alog10(y)

        	oplot, x, y, psym=-3, linestyle=2, color= 0, thick= 2.0

	endif



;-----------------------
;  Put data on graph
;-----------------------2
snapnum= 2
;snapnum= 10
;snapnum= 20

r_kenn= 2.0


;do_altsfr= 1
do_altsfr= 0

if do_altsfr eq 1 then begin
	;  N = 2
	; --------
	symsize= 1.5
	usersym,symsize*[-1,-1,1,1],symsize*[-1,1,1,-1],/fill
	frun= "/raid4/tcox/isolated/d3_N2/"
	ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
	if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
	oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 150, thick=2.0

	frun= "/raid4/tcox/isolated/d3_N2a/"
	ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
	if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
	oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 150, thick=2.0

	frun= "/raid4/tcox/isolated/d3_N2b"
	ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
	if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
	oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 150, thick=2.0


	;  N = 1.5
	; --------
	usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0   ;,/fill
	frun= "/raid4/tcox/isolated/d3_N1.5/"
	ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
	if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
	oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 100, thick=2.0

	frun= "/raid4/tcox/isolated/d3/"
	ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
	if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
	oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 100, thick=2.0




	;  N = 1
	; --------
	usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=3.0, /fill
	frun= "/raid4/tcox/isolated/d3_N0/"
	ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
	if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
	oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 50, thick=2.0

	frun= "/raid4/tcox/isolated/d3_N0a/"
	ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
	if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
	oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 50, thick=2.0

	frun= "/raid4/tcox/isolated/d3_N0b/"
	ok=process_kenn_info(frun,snapnum,r_kenn,gas_sd_log,sfr_sd_log)
	if ok then usersym,symsize*[0,0,-.5,.5,0],symsize*[1,-1,-.5,-.5,-1]
	oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 50, thick=2.0

endif




load_and_plot_seq, "/raid4/tcox/ds/vc3vc3e_2/", 1, r_kenn=r_kenn
load_and_plot_seq, "/raid4/tcox/altsf/vc3vc3e_N1/", 5, r_kenn=r_kenn
load_and_plot_seq, "/raid4/tcox/altsf/vc3vc3e_N2/", 2, r_kenn=r_kenn




;--------------------------------------
device, /close


end






;--------------------------------------
; Actually process the kennicutt point
;--------------------------------------
function process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log, center=center, $
					loadedsnap=loadedsnap


	thiscenter= center

	if not keyword_set(loadedsnap) then begin
           if (fload_snapshot_bh(frun, snapnum, /nopot_in_snap)) then begin
		print, "PROBLEM: opening file"
		return, 0
	   endif
	endif

	rxy = fload_gas_xyz('rxy',center=thiscenter)

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






;==========================================================================






pro kenn_isoline, junk


sendto= 'ps'
filename= 'kenni.eps'

initialize_plotinfo, 1
setup_plot_stuff, sendto, filename=filename, colortable=4


; must run prodir before this, because this routime calls
; process_kenn_line, which 



; ---------------------
;  Set up Parmaters
; ---------------------

xaxistitle="!6Log !7R !D!6gas !N(M!D!9n!6!N pc!E-2!N)"
yaxistitle="!6Log !7R !D!6SFR !N(M!D!9n!6!N yr!E-1!N kpc!E-2!N)"


;xmax = 3.16e5      ; for direct comparison to kennicutt (or volker for that matter)
;xmin = 0.316
;ymax = 3.16e3
;ymin = 3.16e-5

; we'll make the plot log
xmax = 4.5
;xmax = 3.5
;xmax = 2.5
;xmax = 2.3
;xmin = 0.0
;xmin = 0.0
;xmin = -0.4
;xmin = -0.5
xmin = -1.5

ymax = 2.5
;ymax = 1.5
;ymax = 0.0
;ymin = -5.0
ymin = -6.0




; ------------------
; Plot this up
; ------------------
!p.position= [0.20, 0.15, 0.95, 0.95]
;!p.font= 0

plot, [10.0], [10.0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; ---------------------------



	snapnum= 10



	; N= 2
	xyouts, 0.25, 0.85, 'KS_2/iso_b6', color= 150, charthick=4.0, size=1.50, /normal
	load_and_plot, '/n/scratch/hernquist_lab/dnarayanan/gadgetruns/z3/KS_2/iso_b6', snapnum, plotcolor=150, plotthick=4.0;, pstyle= 5
	xyouts, 0.42, 0.85, '', color= 150, charthick=2.0, size=1.10, /normal



	; N= 1.5
	xyouts, 0.25, 0.80, 'iso_b6', color= 100, charthick=4.0, size=1.50, /normal
	load_and_plot, '/n/scratch/hernquist_lab/dnarayanan/gadgetruns/z3/iso_b6', snapnum, plotcolor=100, plotthick=4.0;, pstyle= 5
	xyouts, 0.42, 0.80, '', color= 100, charthick=2.0, size=1.10, /normal



	; n= 1
	xyouts, 0.25, 0.75, 'KS_1/iso_b6', color= 50, charthick=4.0, size=1.50, /normal
	load_and_plot, '/n/scratch/hernquist_lab/dnarayanan/gadgetruns/z3/KS_1/iso_b6', snapnum, plotcolor=50, plotthick=4.0;, pstyle= 2
	xyouts, 0.42, 0.75, '', color= 50, charthick=2.0, size=1.10, /normal




; ---------------------------







ok=oplot_kenn_log(sendto,/linear)
;ok=oplot_kenn_log(sendto,/linear,/ditch_boundary)


;xyouts, 0.25, 0.85, fload_timelbl(1,2), /normal, size=1, color= 0, charthick= 2.5



; 
; Legend
; ---------
;display_sfrlaw_legend= 1
display_sfrlaw_legend= 0
if display_sfrlaw_legend eq 1 then begin
	symsize= 1.5
	oplot, [1.0], [-4.4], symsize= 1.5, psym=2, linestyle=0, color= 100, thick=2.0

	usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0, /fill
	oplot, [1.0], [-4.8], psym=8, linestyle=0, color= 150, thick=2.0

	oplot, [1.0], [-5.2], psym=5, symsize=1.5, linestyle=0, color=50, thick=2.0
	
	xyouts, 0.8, -4.0, '!6SFR !9? !7q!6!EN!N', color= 0, charthick=4.0, size=2.5

	xyouts, 1.2, -4.5, '!6N= 1.0', color= 100, charthick=4.0, size=1.53

	xyouts, 1.2, -4.9, '!6N= 1.5', color= 150, charthick=4.0, size=1.53

	xyouts, 1.2, -5.3, '!6N= 2.0', color= 50, charthick=4.0, size=1.53
endif




;display_iso_legend= 1
display_iso_legend= 0
if display_iso_legend eq 1 then begin
        oplot, [-0.1], [-1.0], symsize= 1.5, psym=2, linestyle=0, color= 0, thick=2.0
	oplot, [-0.25, 0.05], [-1.0,-1.0], psym=-3, linestyle=0, color= 0, thick=4.0

        ;oplot, [-0.1], [-1.4], psym=3, symsize=1.5, linestyle=0, color=0, thick=2.0
	oplot, [-0.25, 0.05], [-1.4,-1.4], psym=-3, linestyle=0, color= 0, thick=4.0

        oplot, [-0.1], [-1.8], psym=5, symsize=1.5, linestyle=0, color=0, thick=2.0
	oplot, [-0.25, 0.05], [-1.8,-1.8], psym=-3, linestyle=0, color= 0, thick=4.0

        xyouts, 0.1, -1.1, '!8low!6', color= 0, charthick=4.0, size=1.53

        xyouts, 0.1, -1.5, '!8med!6', color= 0, charthick=4.0, size=1.53

        xyouts, 0.1, -1.9, '!8high!6', color= 0, charthick=4.0, size=1.53
endif






if sendto eq 'ps' then device, /close




end









pro load_and_plot, frun, snapnum, plotcolor=plotcolor, $
				plotlinest=plotlinest, $
				plotthick=plotthick, $
				pstyle=pstyle

	bins = 20

	plotlinest= 0    ; didn't work that will with other linestyles

	ok=fload_snapshot_bh(frun,snapnum,/nopot_in_snap)
	a= fload_gas_xyz('rxy')
	b= fload_gas_sfr(1)
	c= fload_gas_mass(1)
	d= 0.0 + 0.0*c

	process_kenn_line, a, b, c, d, bins, mass_sd, sfr_sd

	mass_sd = alog10(mass_sd)
	sfr_sd = alog10(sfr_sd)

print, "SFR SD", max(sfr_sd), min(sfr_sd)
print, "Mass SD", max(mass_sd), min(mass_sd)

	;ppsym= -3
	ppsym= -5
	psymsize= 1.5
	if not keyword_set(pstyle) then pstyle= 0 else ppsym= -pstyle

	if pstyle eq 10 then begin
		ppsym= -8
		psymsize= 2.0
		usersym,[-1,-1,1,1,-1],[-1,1,1,-1,-1],thick=4.0,/fill
	endif


	oplot, mass_sd, sfr_sd, psym=ppsym, linestyle= plotlinest, $
				color= plotcolor, thick= plotthick, $
				symsize= psymsize



end





pro load_and_plot_seq, frun, pointt, bhid=bhid, r_kenn=r_kenn


	starti= 0
	spawn, "/bin/ls "+frun+"/snap* | wc ",result
	endi=long(result[0])-1

	print, "frun= ", frun
	print, "starti= ", starti
	print, "endi= ", endi

	select_thispoint, pointt, thispsym, thiscolor

	if not keyword_set(bhid) then begin
		;bhid= fload_blackhole_id(1)
		;bhid= bhid[0]
        	;bhid= bhid[1]
        	;bhid= 200001L  
        	;bhid= 280002L   ; used for z3/b4e
        	bhid= 400002L   ; used for ds/vc3vc3e_2
	endif


	for i=starti,endi do begin

		thisi= i

		ok=fload_snapshot_bh(frun, thisi, /nopot_in_snap)

        	center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)

        	print, "Blackhole ID: ", bhid
        	print, "Blackhole center: ", center_bh

	        ok=process_kenn_info(frun,thisi,r_kenn,gas_sd_log,sfr_sd_log,center=center_bh,/loadedsnap)

	        if ok then begin
			; upper limit
			usersym,1.2*[0,0,-.5,.5,0],1.2*[1,-1,-.5,-.5,-1]
	        	oplot, gas_sd_log, sfr_sd_log, psym=8, color= thiscolor, thick=2.0
		endif else begin
	        	oplot, gas_sd_log, sfr_sd_log, psym=thispsym, color= thiscolor, thick=2.0
		endelse

	endfor

end


