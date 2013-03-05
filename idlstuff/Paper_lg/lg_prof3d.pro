;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   Compiled 3D Profiles
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



; ----------------------------
;  Dark Matter Profile
; ----------------------------
pro lg_prof3d, junk


frun='/raid4/tcox/localgroup/v5'
smoothlen=0.1
filename=frun+'/prof3d.eps'
snapnum=105



if not keyword_set(smoothlen) then smoothlen=0.1
if not keyword_set(filename) then filename=frun+'/lgprof3d.eps'
if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: lg_prof3d"
   print, "  " 
   return
endif


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


                
; ----------------------------- 
; set up constants      
; -----------------------------

bins = 40
                
;xmax = 3
;xmin = -1
;xmax = 400
xmax = 800
xmin = 0.08

ymax = 1e+0
ymin = 1e-8

;xaxistitle= 'Log r (kpc)'
xaxistitle= '!6r (kpc)'
;yaxistitle= '!7q!3!N (r) '
yaxistitle= '!7q!6!N (10 M!D!9n!6!N pc!E-3!N) '

;divide_by_rhocrit= 0
divide_by_rhocrit= 1

if divide_by_rhocrit eq 1 then begin
        yaxistitle="!7q!6!N(r) / !7q!6!Dcrit!N"
        ymax = 5e+8
        ymin = 1e-2
endif


;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
        /xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.7, ycharsize=1.7, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; ------------------------
;  determine dm profile
; ------------------------

ok=fload_snapshot_bh(frun,snapnum)



processandplot_one_dm_profile, frun, snapnum, xmin, xmax, bins, linecolor= 0, ics=ics, $
				divide_by_rhocrit=divide_by_rhocrit, /nofit


processandplot_one_dm_profile, frun, snapnum, xmin, xmax, bins, linecolor= 150, ics=ics, $
				divide_by_rhocrit=divide_by_rhocrit, /nofit, /allstars


processandplot_one_dm_profile, frun, snapnum, xmin, xmax, bins, linecolor= 50, ics=ics, $
				divide_by_rhocrit=divide_by_rhocrit, /nofit, /gas





; -----------------
;  Plot Extras
; -----------------

x0= 0.40
y0= 0.26
oplot, [0.17, 0.5] , [9.0,9.0] , psym=-3, linestyle=0, color= 0, thick= 6.0
xyouts, x0, y0+0.10, '!6Dark Matter', /normal, charthick= 3.0, size= 1.7, color= 0
oplot, [0.17, 0.5] , [1.8,1.8] , psym=-3, linestyle=2, color= 150, thick= 6.0
xyouts, x0, y0+0.05, '!6All Stars', /normal, charthick= 3.0, size= 1.7, color= 150
oplot, [0.17, 0.5] , [0.4,0.4] , psym=-3, linestyle=1, color= 50, thick= 6.0
xyouts, x0, y0,      '!6Gas', /normal, charthick= 3.0, size= 1.7, color= 50



; smoothing length
;ok = oplot_smoothing_length(smoothlen, /plog)


xyouts, 0.8, 0.9, '(a)', /normal, size=1.8, color= 0, charthick=3.0




;--------------------------------------
;--------------------------------------

device, /close



end







;--------------------------------------------------------------------------




pro processandplot_one_dm_profile, frun, snapnum, xmin, xmax, bins, $
                                divide_by_rhocrit=divide_by_rhocrit, $ 
                                allstars=allstars, $
                                density_in_cm3=density_in_cm3, $
                                linecolor=linecolor, $
                                xrwt=xrwt, $
				ics=ics, nofit=nofit, gas=gas

;if keyword_set(ics) then begin
;	ok=fload_initial_conditions(frun)
;endif else begin
;	ok=fload_snapshot_bh(frun,snapnum)
;endelse

; total dark matter profile
; ---------------------------
a=fload_halo_xyz('r')
c=fload_halo_mass(1)

if keyword_set(allstars) then begin
	a= fload_allstars_xyz('r')
	c= fload_allstars_mass(1)
endif

if keyword_set(gas) then begin
	a= fload_gas_xyz('r')
	c= fload_gas_mass(1)
endif

a= a/0.7
c=c/0.7

print, "Dark Matter Properties"
print, "N_dm= ", n_elements(a)
print, "radius (min/max)= ", min(a),max(a)
print, "mass (min/max)= ",min(c), max(c)

r_log_dm= fltarr(bins)
mass_rho_dm= fltarr(bins)

tempxmin= alog10(xmin)
tempxmax= alog10(xmax)
process_rho_profile, a, c, bins, tempxmax, tempxmin, r_log_dm, mass_rho_dm, weight


r_dm= 10^(r_log_dm)    ; non-log version for fitting


idx= where(mass_rho_dm gt 0)
if idx(0) ne -1 then begin
   mass_rho_dm= mass_rho_dm(idx)
   r_log_dm= r_log_dm(idx)
   r_dm= r_dm(idx)
   weight= weight(idx)
endif



;---------------------
; divide by rho_crit
;---------------------
if keyword_set(divide_by_rhocrit) then begin
        rho_crit= fload_rho_crit(1)
        mass_rho_dm= mass_rho_dm/rho_crit
        weight= weight/rho_crit
endif




; ------------------------------------
; Do a curve fit to a NFW profile
; ------------------------------------
if not keyword_set(nofit) then begin
	fitandoverplot_nfw, r_dm, mass_rho_dm, weight, xmin=xmin, xmax=xmax, $
				divide_by_rhocrit=divide_by_rhocrit
endif


thislinestyle= 0
if linecolor eq 50 then thislinestyle= 1
if linecolor eq 150 then thislinestyle= 2


oplot, r_dm, mass_rho_dm, psym=-3, linestyle=thislinestyle, color= linecolor, thick= 6.0
;oploterror, r_dm, mass_rho_dm, weight, psym=-3, linestyle=0, color= linecolor, errcolor= linecolor, errthick=3.0, thick= 3.0

	

end








;============================================================================







