;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   Compiled Gas Profiles
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------






;========================================================================
;========================================================================





; 3D




; -----------------------
;  Gas Profile
; -----------------------
pro plot_3D_profile, frun, smoothlen=smoothlen, $
				filename=filename, $
				snapnum=snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_density_profile_gas"
   print, "  "
   return                       
endif                           

if not keyword_set(snapnum) then snapnum= 25
if not keyword_set(filename) then filename= 'gasprofile.eps'
                
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4



                
                
; ----------------------------- 
; set up constants      
; -----------------------------
                        
bins = 50
                
;xmax = 3
;xmin = -1
xmax = 50
xmin = 0.05

ymax = 1e+8    
ymin = 1e-5

xaxistitle= 'Log R (kpc)'
yaxistitle= '!7q!3!N(cm!E-3!N) '

divide_by_rhocrit= 0
;divide_by_rhocrit= 1

if divide_by_rhocrit eq 1 then begin
	yaxistitle="!7q!3!N(r) / !7q!3!Dcrit!N"
	ymax = 1e+10    
	ymin = 1e-2
endif


;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
	/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; ------------------------
;  determine gas profile
; ------------------------

processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, linecolor= 0, /density_in_cm
;processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, linecolor= 100, /density_in_cm, /profavg, /gasrho


; just hot gas
;processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, /hotgasonly, linecolor= 150, /density_in_cm


; just cold gas
processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, /coldgasonly, TurbulentFactor= 1.0, $
					linecolor= 50, /density_in_cm, /profavg
processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, /coldgasonly, TurbulentFactor= 10.0, $
					linecolor= 150, /density_in_cm, /profavg, /tabulated_profile
processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, /coldgasonly, TurbulentFactor= 100.0, $
					linecolor= 100, /density_in_cm, /profavg




; -----------------
;  Plot Extras
; -----------------

;x0= 0.65
;y0= 0.86
;xyouts, x0, y0, fload_fid(1), /normal, charthick= 2.0, size= 1.7, color= 0

;xyouts, 0.40, 0.40, 'Total', charthick=2.0, size=1.2, color=50
;xyouts, 0.55, 0.72, 'Cold', charthick=2.0, size=1.2, color=100

xyouts, 0.64, 0.87, '!6P!DTurbulent!N/P!DThermal!N', /normal, charthick= 2.0, size= 1.2, color= 0
xyouts, 0.64, 0.83, '1.0', /normal, charthick= 2.0, size= 1.2, color= 50
xyouts, 0.64, 0.79, '10.0', /normal, charthick= 2.0, size= 1.2, color= 150
xyouts, 0.64, 0.75, '100.0', /normal, charthick= 2.0, size= 1.2, color= 100


; smoothing length
ok = oplot_smoothing_length(smoothlen, /plog)


;--------------------------------------
;--------------------------------------

device, /close



end




;========================================================================
;========================================================================




; --------------------------
;  Gas Profile - comparison
; --------------------------
pro plot_3D_profile_comparison, junk, $
				smoothlen=smoothlen, $
				filename=filename, $
				snapnum=snapnum

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: plot_density_profile_comparison"
   print, "  "
   return                       
endif                           

if not keyword_set(filename) then filename='profgas.eps'
if not keyword_set(smoothlen) then smoothlen=0.1

                
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4



; ----------
;  snapshot
; ----------
if not keyword_set(snapnum) then snapnum= 25


                
                
; ----------------------------- 
; set up constants      
; -----------------------------
                        
bins = 50
                
;xmax = 3
;xmin = -1
;xmax = 300
xmax = 10
xmin = 0.1

ymax = 1e+2    
ymin = 1e-10

xaxistitle= '!6Log R (kpc)'
yaxistitle= '!7q!6!N(r) '

divide_by_rhocrit= 0
;divide_by_rhocrit= 1

if divide_by_rhocrit eq 1 then begin
	yaxistitle="!7q!6!N(r) / !7q!6!Dcrit!N"
	ymax = 1e+10    
	ymin = 1e-2
endif

density_in_cm3= 1
;density_in_cm3= 0
if density_in_cm3 eq 1 then begin
	yaxistitle= "!7q!6!N (cm!E-3!N)"
	ymax= 5.0e+3
	ymin= 1.0e-3
endif

;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
	/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
	ytickforma='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; ------------------------
;  determine gas profile
; ------------------------

;frun="/raid4/tcox/vc3vc3e"
;frun="/raid4/tcox/vc3vc3e_2"
;processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, linecolor= 150, $
;					density_in_cm3=density_in_cm3
;
;frun="/raid4/tcox/vc3vc3e_no"
;processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, linecolor= 50, $
;                                        density_in_cm3=density_in_cm3


;frun="/raid2/tcox/grpA_4"
;processandplot_onedensityprofile, frun, 0, xmin, xmax, bins, linecolor= 150, $
;                                       density_in_cm3=density_in_cm3
;processandplot_onedensityprofile, frun, 10, xmin, xmax, bins, linecolor= 120, $
;                                       density_in_cm3=density_in_cm3
;processandplot_onedensityprofile, frun, 20, xmin, xmax, bins, linecolor= 100, $
;                                       density_in_cm3=density_in_cm3
;processandplot_onedensityprofile, frun, 30, xmin, xmax, bins, linecolor=  80, $
;                                       density_in_cm3=density_in_cm3
;processandplot_onedensityprofile, frun, 40, xmin, xmax, bins, linecolor=  60, $
;                                       density_in_cm3=density_in_cm3
;processandplot_onedensityprofile, frun, 50, xmin, xmax, bins, linecolor=  40, $
;                                       density_in_cm3=density_in_cm3
;


frun="/raid4/tcox/ds/vc3vc3e_2"

processandplot_onedensityprofile, frun, 34, xmin, xmax, bins, linecolor=  50, /density_in_cm
xyouts, 0.25, 0.38, 'pre-sb', /normal, size=1.2, color= 50, charthick= 2.5

processandplot_onedensityprofile, frun, 40, xmin, xmax, bins, linecolor=  150, /density_in_cm
xyouts, 0.25, 0.34, 'peak-sb', /normal, size=1.2, color= 150, charthick= 2.5
processandplot_onedensityprofile, frun, 40, xmin, xmax, bins, linecolor=  0, /density_in_cm, /diffusegasonly

processandplot_onedensityprofile, frun, 52, xmin, xmax, bins, linecolor=  100, /density_in_cm
xyouts, 0.25, 0.30, 'post-sb', /normal, size=1.2, color= 100, charthick= 2.5





; -----------------
;  Plot Extras
; -----------------

; smoothing length
x=[smoothlen,smoothlen]
y=[ymin,ymax]
oplot, x, y, linestyle=1, color= 0





;--------------------------------------
;--------------------------------------

device, /close



end





;========================================================================
;========================================================================






;  determine gas profile
; ------------------------
pro processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, $
				divide_by_rhocrit=divide_by_rhocrit, $
				diffusegasonly=diffusegasonly, $
				hotgasonly=hotgasonly, $
				coldgasonly=coldgasonly, $
				density_in_cm=density_in_cm, $
				linecolor=linecolor, $
				profavg=profavg, $
				gasrho=gasrho, $
				TurbulentFactor=TurbulentFactor, $
				tabulated_profile=tabulated_profile, $
				xrwt=xrwt

ok=fload_snapshot_bh(frun,snapnum)

; determine center
center= [0,0,0]
if fload_npart(5) gt 1 then begin
        bhid= fload_blackhole_id(1)
        bhid1= bhid[0]
        bhid2= bhid[1]
        print, "Blackhole ID's: ", bhid1, bhid2
        center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
        center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
        center_bh= center1
        ;center_bh= center2
endif else begin
        bhid= fload_blackhole_id(1)
        center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
endelse


print, "center_bh= ", center_bh
a=fload_gas_xyz('r',center=center_bh)
;a=fload_gas_xyz('r')
a= a/0.7
c=fload_gas_mass(1)
c=c/0.7

if keyword_set(hotgasonly) then begin
	u=fload_gas_u(1)

	idx= where(u gt 150.0)
	if idx(0) ne -1 then begin
		a= a(idx)
		c= c(idx)
	endif else begin
		print, "No HOT particles!"
		return
	endelse
endif


if keyword_set(diffusegasonly) then begin
	print, "mass max/min= ", max(c), min(c)
	cf= fload_gas_coldfraction(1)
	print, "coldfraction max/min= ", max(cf), min(cf)
	c= c * (1.0 - cf)
	print, "diffuse gas mass max/min= ", max(c), min(c)
endif


if keyword_set(coldgasonly) then begin
	c= fload_gas_rho(1,/cold,TurbulentFactor=TurbulentFactor)
	print, "cold density max/min= ", max(c), min(c)
endif

if keyword_set(gasrho) then begin
	c= fload_gas_rho(1)
endif

if keyword_set(xrwt) then begin
	xray=fload_gas_xray_luminosity(1,/diffuse_hotgas)
	idx= where(xray gt 0.0)
	if idx(0) ne -1 then begin
		a= a(idx)
		c= c(idx)
	endif else begin
		print, "No HOT particles!"
		return
	endelse
endif

r_log_gas = fltarr(bins)
mass_sd_gas = fltarr(bins)

tempxmin= alog10(xmin)
tempxmax= alog10(xmax)

if not keyword_set(profavg) then begin
	a= alog10(a)
	process_prof_3drho, a, c, bins, tempxmax, tempxmin, r_log_gas, mass_sd_gas, weight, /r_is_log
endif else begin
	process_prof_avg, a, c, bins, tempxmax, tempxmin, r_log_gas, mass_sd_gas, weight
endelse

r_gas= 10^(r_log_gas)    ; non-log version for fitting


idx= where(mass_sd_gas gt 0)
if idx(0) ne -1 then begin
   mass_sd_gas= mass_sd_gas(idx)
   r_log_gas= r_log_gas(idx)
   r_gas= r_gas(idx)
   weight= weight(idx)
endif


if keyword_set(density_in_cm) then begin
	;
	; gadget unit = 6.76991x10-22 g/cm3
	;
	;   then divide by mu*m_p = 0.6 * 1.6726x10-24 g
	;
	;mass_sd_gas= mass_sd_gas * 674.59

	; this is the straight-up conversion
	mass_sd_gas= mass_sd_gas * 404.5
endif

;---------------------
; divide by rho_crit
;---------------------
if keyword_set(divide_by_rhocrit) then begin
	rho_crit= fload_rho_crit(1)
        mass_sd_gas= mass_sd_gas/rho_crit
        weight= weight/rho_crit
endif


if keyword_set(tabulated_profile) then begin
	openw, 1, frun+'/gasprofile.txt', ERROR=err

	printf, 1, "#   gasprofile.txt "
	printf, 1, "#   Profile of (spherically averaged) gas densities. "
	printf, 1, "#             "
	printf, 1, "#  radius  Log gas density"
	printf, 1, "#   (kpc)      (cm^-3)"
	for i=0,n_elements(r_gas)-1 do begin
	        printf, 1, FORMAT= '(F9.3,"     ",1(F9.4,"  "))', $
	                r_gas[i], alog10(mass_sd_gas[i])
	endfor
	close, 1

endif

; this is if it's a log x-axis
;oplot, r_log_gas, mass_sd_gas, psym=-3, linestyle=0, color= 0, thick=2.0


; this is if it's regular r, but with log coordinates
if linecolor eq 150 then begin
        symsize= 1.0
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        oplot, r_gas, mass_sd_gas, psym=-8, linestyle=0, color= linecolor, thick=3.0
endif

if linecolor eq 50 then begin
        symsize= 1.2
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=3.0
        oplot, r_gas, mass_sd_gas, psym=-8, linestyle=2, color= linecolor, thick=3.0
endif

if linecolor eq 100 then begin
        symsize= 1.0
        oplot, r_gas, mass_sd_gas, psym=-2, linestyle=1, color= linecolor, thick=3.0
endif

if linecolor eq 0 then begin
        symsize= 1.0
        oplot, r_gas, mass_sd_gas, psym=-3, linestyle=1, color= linecolor, thick=2.0
endif


end





;=================================================================================






; ---------------------------------
;  Cold Gas Volume Filling Factor
; ---------------------------------
pro plot_coldgasVFF_profile, frun, smoothlen=smoothlen, $
				filename=filename, $
				snapnum=snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "plot_coldgasVFF_profile, frun, "
   print, "  "
   return                       
endif                           


if not keyword_set(filename) then filename='coldgVFF.eps'
if not keyword_set(snapnum) then snapnum= 31
                
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4



                
                
; ----------------------------- 
; set up constants      
; -----------------------------
                        
bins = 50
                
;xmax = 3
;xmin = -1
xmax = 50
xmin = 0.05

ymax = 1e+0    
ymin = 1e-4

xaxistitle= 'Log R (kpc)'
yaxistitle= 'Cold Gas Volume Filling Factor'



;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
	/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
	ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; ------------------------
;  determine gas profile
; ------------------------

processandplot_onecoldgVFF_profile, frun, snapnum, xmin, xmax, bins, TurbulentFactor=1.0, $
					linecolor= 50, /profavg

processandplot_onecoldgVFF_profile, frun, snapnum, xmin, xmax, bins, TurbulentFactor=10.0, $
					linecolor= 150, /profavg, /tabulated_profile

processandplot_onecoldgVFF_profile, frun, snapnum, xmin, xmax, bins, TurbulentFactor=100.0, $
					linecolor= 100, /profavg



;processandplot_onecoldgVFF_profile, "/raid4/tcox/vc3vc3e_2", 54, xmin, xmax, bins, $
;					linecolor= 150, /profavg




; -----------------
;  Plot Extras
; -----------------

x0= 0.65
y0= 0.86
;xyouts, x0, y0, fload_fid(1), /normal, charthick= 2.0, size= 1.7, color= 0

xyouts, 0.24, 0.87, '!6P!DTurbulent!N/P!DThermal!N', /normal, charthick= 2.0, size= 1.2, color= 0
xyouts, 0.24, 0.83, '1.0', /normal, charthick= 2.0, size= 1.2, color= 50
xyouts, 0.24, 0.79, '10.0', /normal, charthick= 2.0, size= 1.2, color= 150
xyouts, 0.24, 0.75, '100.0', /normal, charthick= 2.0, size= 1.2, color= 100



; smoothing length
ok = oplot_smoothing_length(smoothlen, /plog)




;--------------------------------------
;--------------------------------------

device, /close



end






;========================================================================
;========================================================================





;  determine profile
; ------------------------
pro processandplot_onecoldgVFF_profile, frun, snapnum, xmin, xmax, bins, $
				linecolor=linecolor, $
				tabulated_profile=tabulated_profile, $
				profavg=profavg, TurbulentFactor=TurbulentFactor


ok=fload_snapshot_bh(frun,snapnum)

; determine center
center= [0,0,0]
if fload_npart(5) gt 1 then begin
        bhid= fload_blackhole_id(1)
        bhid1= bhid[0]
        bhid2= bhid[1]
        print, "Blackhole ID's: ", bhid1, bhid2
        center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
        center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
        center_bh= center1
        ;center_bh= center2
endif else begin
        bhid= fload_blackhole_id(1)
        center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
endelse


print, "center_bh= ", center_bh
a=fload_gas_xyz('r',center=center_bh)
;a=fload_gas_xyz('r')
a= a/0.7
c=fload_gas_coldgasfillingfactor(1,TurbulentFactor=TurbulentFactor)
c=c/0.7

print, "cold gas volume filling factor max/min= ", max(c), min(c)



r_log_gas = fltarr(bins)
mass_sd_gas = fltarr(bins)

tempxmin= alog10(xmin)
tempxmax= alog10(xmax)

if not keyword_set(profavg) then begin
	process_prof_rho, a, c, bins, tempxmax, tempxmin, r_log_gas, mass_sd_gas, weight
endif else begin
	process_prof_avg, a, c, bins, tempxmax, tempxmin, r_log_gas, mass_sd_gas, weight
endelse

r_gas= 10^(r_log_gas)    ; non-log version for fitting


idx= where(mass_sd_gas gt 0)
if idx(0) ne -1 then begin
   mass_sd_gas= mass_sd_gas(idx)
   r_log_gas= r_log_gas(idx)
   r_gas= r_gas(idx)
   weight= weight(idx)
endif


if keyword_set(tabulated_profile) then begin
	openw, 1, frun+'/coldgasVFF.txt', ERROR=err

	printf, 1, "#   coldgasVFF.txt "
	printf, 1, "#   Profile of (spherically averaged) cold gas volume filling fraction. "
	printf, 1, "#        "
	printf, 1, "#  radius  log cold gas volume"
	printf, 1, "#   (kpc)    filling fraction"
	for i=0,n_elements(r_gas)-1 do begin
	        printf, 1, FORMAT= '(F9.3,"       ",1(F9.4,"  "))', $
	                r_gas[i], alog10(mass_sd_gas[i])
	endfor
	close, 1

endif



; this is if it's regular r, but with log coordinates
if linecolor eq 150 then begin
        symsize= 1.0
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        oplot, r_gas, mass_sd_gas, psym=-8, linestyle=0, color= linecolor, thick=3.0
endif

if linecolor eq 50 then begin
        symsize= 1.2
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=3.0
        oplot, r_gas, mass_sd_gas, psym=-8, linestyle=2, color= linecolor, thick=3.0
endif

if linecolor eq 100 then begin
        symsize= 1.0
        oplot, r_gas, mass_sd_gas, psym=-2, linestyle=1, color= linecolor, thick=3.0
endif




end





;=================================================================================
;=================================================================================










; plot cumulative mass profile


pro plot_gas_mass_comparison, junk, smoothlen=smoothlen, filename=filename, $
				snapnum=snapnum
        
if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: plot_gas_profile"
   print, "  "
   return
endif   
        
if not keyword_set(filename) then filename="gas_cm.eps"
if not keyword_set(smoothlen) then smoothlen=0.1


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


; ----------
;  snapshot
; ----------
if not keyword_set(snapnum) then snapnum= 25


; ----------------------------- 
; set up constants      
; -----------------------------
                        
bins = 50
                
;xmax = 3
;xmin = -1
xmax = 300
xmin = 0.1

ymax = 1e+11
ymin = 1e+9

xaxistitle= 'R (h!E-1!N kpc)'
;yaxistitle= 'M!Dgas!N (M!D!9n!3!N)'
yaxistitle= 'M!Dstars!N (M!D!9n!3!N)'


;---------------------------
;  Print it
;---------------------------
   
!p.position= [0.2, 0.15, 0.95, 0.95]
   
plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
        /xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata

; ------------------------
;  determine gas profile
; ------------------------

;frun="/raid4/tcox/vc3vc3e"
frun="/raid4/tcox/vc3vc3e_2"
processandplot_onemassprofile, frun, snapnum, xmin, xmax, bins, linecolor= 150


frun="/raid4/tcox/vc3vc3e_no"
processandplot_onemassprofile, frun, snapnum, xmin, xmax, bins, linecolor= 50


; -----------------
;  Plot Extras
; -----------------

; smoothing length
x=[smoothlen,smoothlen]
y=[ymin,ymax]
oplot, x, y, linestyle=1, color= 0


;--------------------------------------
;--------------------------------------

device, /close



end







; actually do the work
; --------------------------
pro processandplot_onemassprofile, frun, snapnum, xmin, xmax, bins, $
					linecolor=linecolor


ok=fload_snapshot_bh(frun,snapnum)
;radius=fload_gas_xyz('r')
;mass=fload_gas_mass(1)
radius=fload_allstars_xyz('r')
mass=fload_allstars_mass(1)
print, "total mass= ", total(mass)
mass= 1.0d+10 * mass

rs = fltarr(bins)
cumulative_mass= fltarr(bins)

tempxmin= alog10(xmin)
tempxmax= alog10(xmax)

;--------------------
; process the cumulative profile

binsize = float((tempxmax-tempxmin))/bins

for i=0, bins-1 do begin
    currentr_log= tempxmin + (i+1)*binsize
    currentr= 10^(currentr_log)
    idx= where(radius le currentr)
    rs[i]= currentr
    cumulative_mass[i]= total(mass(idx))
endfor

oplot, rs, cumulative_mass, psym=-3, linestyle=0, color= linecolor, thick=4.0


end




;==================================================================================







