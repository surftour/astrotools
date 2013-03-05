;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   Compiled 3D Profiles
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------











;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     density profiles
;     -------------------------------------------
;  plot 1:  just dark matter density profile, eventually will want to fit NFW to it
;  plot 2:  baryonic density profile, by component
;  plot 3:  the ratio of baryonic to dark matter
;  
;  plus we define a process to do all the profiles the say matter
;  
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------






; ----------------------------
;  Dark Matter Profile
; ----------------------------
pro plot_density_profile_dm, frun, smoothlen=smoothlen, filename=filename, $
					ics=ics, snapnum=snapnum



if not keyword_set(smoothlen) then smoothlen=0.1
if not keyword_set(filename) then filename='prof.eps'
if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_density_profile_dm"
   print, "  " 
   return
endif


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


; ----------
;  snapshot
; ----------
;if not keyword_set(snapnum) then snapnum= 25
if not keyword_set(snapnum) then snapnum= 0


                
                
; ----------------------------- 
; set up constants      
; -----------------------------

bins = 50
                
;xmax = 3
;xmin = -1
;xmax = 400
xmax = 800
xmin = 0.04

ymax = 1e+0
ymin = 1e-8

;xaxistitle= 'Log r (kpc)'
xaxistitle= 'r (kpc)'
;yaxistitle= '!7q!3!N (r) '
yaxistitle= '!7q!3!N (10 M!D!9n!3!N pc!E-3!N) '

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
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; ------------------------
;  determine dm profile
; ------------------------

processandplot_one_dm_profile, frun, snapnum, xmin, xmax, bins, linecolor= 50, ics=ics, $
				divide_by_rhocrit=divide_by_rhocrit


;frun= frun2
;processandplot_one_dm_profile, frun, snapnum, xmin, xmax, bins, linecolor= 150


processandplot_one_dm_profile, frun, snapnum, xmin, xmax, bins, linecolor= 100, ics=ics, $
				divide_by_rhocrit=divide_by_rhocrit, /nofit, /allstars





; -----------------
;  Plot Extras
; -----------------

x0= 0.55
y0= 0.86
;xyouts, x0, y0, fload_fid(1), /normal, charthick= 2.0, size= 1.7, color= 0
;xyouts, x0-0.23, y0+0.04, "/home/tcox/MakeNewDisk/vc/vc3c.dat", /normal, charthick= 2.0, size= 1.1, color= 0
y0= 0.80
xyouts, x0, y0+0.05, 'Dark Matter Profile', /normal, charthick= 2.0, size= 1.2, color= 0
xyouts, x0, y0, 'Stellar Profile', /normal, charthick= 2.0, size= 1.2, color= 100



; smoothing length
;ok = oplot_smoothing_length(smoothlen, /plog)



;--------------------------------------
;--------------------------------------

device, /close



end







;--------------------------------------------------------------------------




pro dmprof_c, junk, filename=filename


if not keyword_set(filename) then filename='prof.eps'
if not keyword_set(junk) then begin
   print, "  "
   print, "dmprof_c, junk, filename=filename"
   print, "  " 
   return
endif


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


smoothlen= 0.1

                
                
; ----------------------------- 
; set up constants      
; -----------------------------

bins = 50
                
;xmax = 3
;xmin = -1
xmax = 400
xmin = 0.04

ymax = 1e+0
ymin = 1e-8

;xaxistitle= 'Log r (kpc)'
xaxistitle= 'r (!8h!3!E-1!Nkpc)'
;yaxistitle= '!7q!3!N (r) '
yaxistitle= '!7q!3!N (10 M!D!9n!3!N pc!E-3!N) '

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

!p.position= [0.2, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
        /xlog, $
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


; ------------------------
;  determine dm profile
; ------------------------

;xyouts, 0.38, 0.91, '250,000 DM particles', /normal, charthick= 2.0, size= 1.7, color= 0
;processandplot_one_dm_profile, '/raid2/tcox/btest/250k_1', 0, xmin, xmax, bins, linecolor= 50, /nofit
;xyouts, 0.7, 0.85, 'initial', /normal, charthick= 2.0, size= 1.2, color= 50
;processandplot_one_dm_profile, '/raid2/tcox/btest/250k_1', 5, xmin, xmax, bins, linecolor= 150, /nofit
;xyouts, 0.7, 0.81, '1 Gyr evolution', /normal, charthick= 2.0, size= 1.2, color= 150

; -------------------
;xyouts, 0.37, 0.91, '2,500,000 DM particles', /normal, charthick= 2.0, size= 1.7, color= 0
;processandplot_one_dm_profile, '/raid2/tcox/btest/2M_1', 0, xmin, xmax, bins, linecolor= 50, /nofit
;xyouts, 0.7, 0.85, 'initial', /normal, charthick= 2.0, size= 1.2, color= 50
;processandplot_one_dm_profile, '/raid2/tcox/btest/2M_1', 5, xmin, xmax, bins, linecolor= 150, /nofit
;xyouts, 0.7, 0.81, '1 Gyr evolution', /normal, charthick= 2.0, size= 1.2, color= 150
;xyouts, 0.4, 0.34, 'std', /normal, charthick= 2.0, size= 1.2, color= 0

; -------------------
xyouts, 0.37, 0.91, '2,500,000 DM particles', /normal, charthick= 2.0, size= 1.7, color= 0
processandplot_one_dm_profile, '/raid2/tcox/btest/2M_3', 0, xmin, xmax, bins, linecolor= 50, /nofit
xyouts, 0.67, 0.85, 'initial', /normal, charthick= 2.0, size= 1.2, color= 50
processandplot_one_dm_profile, '/raid2/tcox/btest/2M_3', 5, xmin, xmax, bins, linecolor= 150, /nofit
xyouts, 0.67, 0.81, '1 Gyr/h evolution', /normal, charthick= 2.0, size= 1.2, color= 150
processandplot_one_dm_profile, '/raid2/tcox/btest/2M_3', 10, xmin, xmax, bins, linecolor= 100, /nofit
xyouts, 0.67, 0.77, '5 Gyr/h evolution', /normal, charthick= 2.0, size= 1.2, color= 100
xyouts, 0.4, 0.34, 'lower', /normal, charthick= 2.0, size= 1.2, color= 0
xyouts, 0.4, 0.30, 'integration', /normal, charthick= 2.0, size= 1.2, color= 0
xyouts, 0.4, 0.26, 'accuracy', /normal, charthick= 2.0, size= 1.2, color= 0

; -------------------
;xyouts, 0.37, 0.91, '2,500,000 DM particles', /normal, charthick= 2.0, size= 1.7, color= 0
;processandplot_one_dm_profile, '/raid2/tcox/btest/2M_2', 0, xmin, xmax, bins, linecolor= 50, /nofit
;xyouts, 0.7, 0.85, 'initial', /normal, charthick= 2.0, size= 1.2, color= 50
;processandplot_one_dm_profile, '/raid2/tcox/btest/2M_2', 4, xmin, xmax, bins, linecolor= 150, /nofit
;processandplot_one_dm_profile, '/raid2/tcox/btest/2M_2', 5, xmin, xmax, bins, linecolor= 150, /nofit
;xyouts, 0.7, 0.81, '1 Gyr evolution', /normal, charthick= 2.0, size= 1.2, color= 150
;xyouts, 0.4, 0.34, 'higher', /normal, charthick= 2.0, size= 1.2, color= 0
;xyouts, 0.4, 0.30, 'integration', /normal, charthick= 2.0, size= 1.2, color= 0
;xyouts, 0.4, 0.26, 'accuracy', /normal, charthick= 2.0, size= 1.2, color= 0

; -------------------
;xyouts, 0.37, 0.91, 'standard vc3 disk', /normal, charthick= 2.0, size= 1.7, color= 0
;processandplot_one_dm_profile, '/raid4/tcox/isolated/vc3c', 0, xmin, xmax, bins, linecolor= 50, /nofit
;xyouts, 0.7, 0.85, 'initial', /normal, charthick= 2.0, size= 1.2, color= 50
;processandplot_one_dm_profile, '/raid4/tcox/isolated/vc3c', 30, xmin, xmax, bins, linecolor= 150, /nofit
;xyouts, 0.7, 0.81, '3 Gyr/h evolution', /normal, charthick= 2.0, size= 1.2, color= 150

; -------------------
;xyouts, 0.37, 0.91, 'std disk vs. merg. rem.', /normal, charthick= 2.0, size= 1.7, color= 0
;processandplot_one_dm_profile, '/raid4/tcox/isolated/vc3c', 0, xmin, xmax, bins, linecolor= 50, /nofit
;xyouts, 0.7, 0.85, 'initial disk', /normal, charthick= 2.0, size= 1.2, color= 50
;processandplot_one_dm_profile, '/raid4/tcox/vc3vc3e_2', 106, xmin, xmax, bins, linecolor= 150, /nofit
;xyouts, 0.7, 0.81, 'merger remnant', /normal, charthick= 2.0, size= 1.2, color= 150

; -------------------
;xyouts, 0.37, 0.91, 'rem vs. relaxed rem', /normal, charthick= 2.0, size= 1.7, color= 0
;processandplot_one_dm_profile, '/raid4/tcox/vc3vc3e_2', 71, xmin, xmax, bins, linecolor= 50, /nofit
;xyouts, 0.7, 0.85, 'rem (T=1.3)', /normal, charthick= 2.0, size= 1.2, color= 50
;processandplot_one_dm_profile, '/raid4/tcox/vc3vc3e_2', 106, xmin, xmax, bins, linecolor= 150, /nofit
;xyouts, 0.7, 0.81, 'rem (T=3.0)', /normal, charthick= 2.0, size= 1.2, color= 150




; -----------------
;  Plot Extras
; -----------------

x0= 0.55
y0= 0.86
;xyouts, x0, y0, fload_fid(1), /normal, charthick= 2.0, size= 1.7, color= 0
;xyouts, x0-0.23, y0+0.04, "/home/tcox/MakeNewDisk/vc/vc3c.dat", /normal, charthick= 2.0, size= 1.1, color= 0
y0= 0.80
;xyouts, x0, y0, 'Dark Matter Profile', /normal, charthick= 2.0, size= 1.2, color= 0



; smoothing length
;ok = oplot_smoothing_length(smoothlen, /plog)
;ok = oplot_smoothing_length(smoothlen)
x= [smoothlen,smoothlen]
y= [ymin,ymax]
oplot, x, y, psym=-3, linestyle=2, color= 0

x= x*2.3
oplot, x, y, psym=-3, linestyle=1, color= 0


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
				ics=ics, nofit=nofit

if keyword_set(ics) then begin
	ok=fload_initial_conditions(frun)
endif else begin
	ok=fload_snapshot_bh(frun,snapnum)
endelse

; total dark matter profile
; ---------------------------
a=fload_halo_xyz('r')
c=fload_halo_mass(1)

if keyword_set(allstars) then begin
	a= fload_allstars_xyz('r')
	c= fload_allstars_mass(1)
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
process_prof_rho, a, c, bins, tempxmax, tempxmin, r_log_dm, mass_rho_dm, weight


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



oplot, r_dm, mass_rho_dm, psym=-3, linestyle=0, color= linecolor, thick= 3.0
oploterror, r_dm, mass_rho_dm, weight, psym=-3, linestyle=0, color= linecolor, errcolor= linecolor, errthick=3.0, thick= 3.0

	

end








;============================================================================







; -----------------------
;  Baryonic Profile
; -----------------------
pro plot_density_profile_b, frun, sendto, smoothlen=smoothlen, filename=filename

if not keyword_set(sendto) then sendto='ps'
if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_density_profile_b"
   print, "  "
   return
endif


setup_plot_stuff, sendto, filename=filename



; -----------------------------
; set up constants
; -----------------------------

bins = 100

xmax = 3
xmin = -2
ymax = 1e+12
ymin = 1




        ; total (baryonic) profile
        ; --------------------------
        a=fload_baryon_xyz('r')
        c=fload_baryon_mass(1)

        r_log= fltarr(bins)
        mass_sd = fltarr(bins)

        process_prof_rho, a, c, bins, xmax, xmin, r_log, mass_sd, weight

        ; gas profile
        ; --------------
        if fload_npart(0) GT 0 then begin
                a=fload_gas_xyz('r')
                c=fload_gas_mass(1)

                r_log_gas = fltarr(bins)
                mass_sd_gas = fltarr(bins)

                process_prof_rho, a, c, bins, xmax, xmin, r_log_gas, mass_sd_gas, weight
        endif

        ; disk profile
        ; --------------
        if fload_npart(2) GT 0 then begin
                a=fload_disk_xyz('r')
                c=fload_disk_mass(1)

                r_log_disk = fltarr(bins)
                mass_sd_disk = fltarr(bins)

                process_prof_rho, a, c, bins, xmax, xmin, r_log_disk, mass_sd_disk, weight
        endif

        ; bulge profile
        ; --------------
        if fload_npart(3) GT 0 then begin
                a=fload_bulge_xyz('r')
                c=fload_bulge_mass(1)

                r_log_bulge = fltarr(bins)
                mass_sd_bulge = fltarr(bins)

                process_prof_rho, a, c, bins, xmax, xmin, r_log_bulge, mass_sd_bulge, weight
        endif

        ; new stars profile
        ; -------------------
        if fload_npart(4) GT 0 then begin
                a=fload_newstars_xyz('r')
                c=fload_newstars_mass(1)

                r_log_stars = fltarr(bins)
                mass_sd_stars = fltarr(bins)

                process_prof_rho, a, c, bins, xmax, xmin, r_log_stars, mass_sd_stars, weight
        endif




;---------------------
; divide by rho_crit
;---------------------
mass_sd= mass_sd/fload_rho_crit(1)
if fload_npart(0) gt 0 then mass_sd_gas= mass_sd_gas/fload_rho_crit(1)
if fload_npart(2) gt 0 then mass_sd_disk= mass_sd_disk/fload_rho_crit(1)
if fload_npart(3) gt 0 then mass_sd_bulge= mass_sd_bulge/fload_rho_crit(1)
if fload_npart(4) gt 0 then mass_sd_stars= mass_sd_stars/fload_rho_crit(1)




;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle="Log R (kpc)", $
        ytitle="!7q!3!N(r) / !7q!3!Dcrit!N", $
        /nodata


        oplot, r_log, mass_sd, psym=-3, linestyle=0, color= 0
        if fload_npart(0) GT 0 then oplot, r_log_gas, mass_sd_gas, psym=-3, linestyle=0, color= 50, thick=3.0
        if fload_npart(2) GT 0 then oplot, r_log_disk, mass_sd_disk, psym=-3, linestyle=0, color= 100, thick=3.0
        if fload_npart(3) GT 0 then oplot, r_log_bulge, mass_sd_bulge, psym=-3, linestyle=0, color= 200, thick=3.0
        if fload_npart(4) GT 0 then oplot, r_log_stars, mass_sd_stars, psym=-3, linestyle=0, color= 150, thick=3.0

        x0= 0.65
        y0= 0.9
        xyouts, x0, y0, 'Total', /normal, charthick=1.7, size= 1.7, color= 0
        if fload_npart(0) GT 0 then xyouts, x0, y0-0.04, 'Gas', /normal, charthick= 3.0, size= 1.7, color= 50
        if fload_npart(2) GT 0 then xyouts, x0, y0-0.08, 'Disk', /normal, charthick= 3.0, size= 1.7, color= 100
        if fload_npart(3) GT 0 then xyouts, x0, y0-0.16, 'Bulge', /normal, charthick= 3.0, size= 1.7, color= 200
        if fload_npart(4) GT 0 then xyouts, x0, y0-0.12, 'New Stars', /normal, charthick= 3.0, size= 1.7, color= 150


xyouts, 0.35, 0.85, fload_fid(1), size=1.7, color= 0, /normal


ok = oplot_smoothing_length(smoothlen, /plog)

;--------------------------------------
;--------------------------------------

if (sendto EQ 'ps') then device, /close



end






; -------------------------------------
;  Density Profile Ratio (baryonic/dm)
; -------------------------------------
pro plot_density_profile_ratio, frun, sendto, smoothlen=smoothlen, filename=filename

if not keyword_set(sendto) then sendto='ps'
if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_density_profile_ratio"
   print, "  "
   return
endif


setup_plot_stuff, sendto, filename=filename


bins = 100


xmin= -2
xmax= 3

ymax= 1000
ymin=0.1




        ; total (baryonic) profile
        ; --------------------------
        a=fload_baryon_xyz('r')
        c=fload_baryon_mass(1)

        r_log_b= fltarr(bins)
        mass_rho_b = fltarr(bins)

        process_prof_rho, a, c, bins, xmax, xmin, r_log_b, mass_rho_b, weight


        ; total dark matter profile
        ; ---------------------------
        a=fload_halo_xyz('r')
        c=fload_halo_mass(1)

        r_log_dm= fltarr(bins)
        mass_rho_dm= fltarr(bins)

        process_prof_rho, a, c, bins, xmax, xmin, r_log_dm, mass_rho_dm, weight



;---------------------
; divide by rho_crit
;---------------------
mass_rho_dm= mass_rho_dm/fload_rho_crit(1)
mass_rho_b= mass_rho_b/fload_rho_crit(1)

idx= where(mass_rho_dm gt 0)
if idx(0) ne -1 then begin
        ratios= mass_rho_b(idx)/mass_rho_dm(idx)
        r_ratios= r_log_dm(idx)
endif else begin
        ratios= mass_rho_b/mass_rho_dm
        r_ratios= r_log_dm
endelse




;------------------
; Set Ymax
;------------------
ratioslog= alog10(ratios)
ymax= fix(max(ratioslog)+1)
ymax= 10^(ymax)



; ------------------
; Plot this up
; ------------------
!p.position= [0.24, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle="Log R (kpc)", $
        ytitle="!7q!3!N!Dbaryon!N/!7q!3!Ddm!N", $
        /nodata


        oplot, r_ratios, ratios, psym=-3, linestyle=0, color= 150, thick=3.0


	ok = oplot_smoothing_length(smoothlen, /plog)


        ;plot 1
        x=[xmin,xmax]
	y=[1.0,1.0]
        oplot, x, y, psym=-3, color= 0, linestyle=1, thick=3.0



if sendto eq 'ps' then device, /close



end




; ---------------------------------------------------------------------------------
; ---------------------------------------------------------------------------------





; 3D




; -----------------------
;  Gas Profile
; -----------------------
pro plot_density_profile_gas, frun, sendto, smoothlen=smoothlen, $
				filename=filename, $
				snapnum=snapnum

if not keyword_set(sendto) then sendto='ps'
if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_density_profile_gas"
   print, "  "
   return                       
endif                           

                
initialize_plotinfo, 1
setup_plot_stuff, sendto, filename=filename, colortable= 4



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
xmin = 1.0

ymax = 1e+2    
ymin = 1e-10

xaxistitle= 'Log R (kpc)'
yaxistitle= '!7q!3!N(r) '

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

processandplot_oneprofile, frun, snapnum, xmin, xmax, bins, linecolor= 50


; just hot gas?
processandplot_oneprofile, frun, snapnum, xmin, xmax, bins, /hotgasonly, linecolor= 150







; -----------------
;  Plot Extras
; -----------------

x0= 0.65
y0= 0.86
xyouts, x0, y0, fload_fid(1), /normal, charthick= 2.0, size= 1.7, color= 0
y0= 0.80
xyouts, x0, y0, 'Gas Profile', /normal, charthick= 2.0, size= 1.7, color= 0



; smoothing length
ok = oplot_smoothing_length(smoothlen, /plog)




;--------------------------------------
;--------------------------------------

if (sendto EQ 'ps') then device, /close



end




; --------------------------
;  Gas Profile - comparison
; --------------------------
pro plot_density_profile_gas_comparison, junk, $
				smoothlen=smoothlen, $
				filename=filename, $
				snapnum=snapnum

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: plot_density_profile_comparison"
   print, "  "
   return                       
endif                           

if not keyword_set(filename) then filename='gas.eps'
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
;xmin = 1.0
xmin = 0.1

ymax = 1e+2    
ymin = 1e-10

xaxistitle= 'Log R (kpc)'
yaxistitle= '!7q!3!N!Dgas!N(r) '

divide_by_rhocrit= 0
;divide_by_rhocrit= 1

if divide_by_rhocrit eq 1 then begin
	yaxistitle="!7q!3!N(r) / !7q!3!Dcrit!N"
	ymax = 1e+10    
	ymin = 1e-2
endif

;density_in_cm3= 1
density_in_cm3= 0
if density_in_cm3 eq 1 then begin
	yaxistitle= "!7q!3!N (cm!E-3!N)"
	ymax= 1.0e+1
	ymin= 1.0e-7
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

;frun="/raid4/tcox/vc3vc3e"
frun="/raid4/tcox/vc3vc3e_2"
processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, linecolor= 150, $
					density_in_cm3=density_in_cm3

xyouts, 0.55, 0.85, 'vc3vc3e_2, with AGN', color=150, charthick=3.0, /normal

frun="/raid4/tcox/vc3vc3e_no"
processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, linecolor= 50, $
                                        density_in_cm3=density_in_cm3

xyouts, 0.55, 0.80, 'vc3vc3e_no, without', color=50, charthick=3.0, /normal


;frun="/raid2/tcox/grpA_4"
;processandplot_oneprofile, frun, 0, xmin, xmax, bins, linecolor= 150, $
;                                       density_in_cm3=density_in_cm3
;processandplot_oneprofile, frun, 10, xmin, xmax, bins, linecolor= 120, $
;                                       density_in_cm3=density_in_cm3
;processandplot_oneprofile, frun, 20, xmin, xmax, bins, linecolor= 100, $
;                                       density_in_cm3=density_in_cm3
;processandplot_oneprofile, frun, 30, xmin, xmax, bins, linecolor=  80, $
;                                       density_in_cm3=density_in_cm3
;processandplot_oneprofile, frun, 40, xmin, xmax, bins, linecolor=  60, $
;                                       density_in_cm3=density_in_cm3
;processandplot_oneprofile, frun, 50, xmin, xmax, bins, linecolor=  40, $
;                                       density_in_cm3=density_in_cm3








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







;  determine gas profile
; ------------------------
;pro processandplot_oneprofile, frun, snapnum, xmin, xmax, bins, $
pro processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, $
				divide_by_rhocrit=divide_by_rhocrit, $
				hotgasonly=hotgasonly, $
				density_in_cm3=density_in_cm3, $
				linecolor=linecolor, $
				xrwt=xrwt

ok=fload_snapshot_bh(frun,snapnum)
a=fload_gas_xyz('r')
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
process_prof_rho, a, c, bins, tempxmax, tempxmin, r_log_gas, mass_sd_gas, weight

r_gas= 10^(r_log_gas)    ; non-log version for fitting


idx= where(mass_sd_gas gt 0)
if idx(0) ne -1 then begin
   mass_sd_gas= mass_sd_gas(idx)
   r_log_gas= r_log_gas(idx)
   r_gas= r_gas(idx)
   weight= weight(idx)
endif


if keyword_set(density_in_cm3) then begin
	;
	; gadget unit = 6.76991x10-22 g/cm3
	;
	;   then divide by mu*m_p = 0.6 * 1.6726x10-24 g
	;
	mass_sd_gas= mass_sd_gas * 674.59
endif

;---------------------
; divide by rho_crit
;---------------------
if keyword_set(divide_by_rhocrit) then begin
	rho_crit= fload_rho_crit(1)
        mass_sd_gas= mass_sd_gas/rho_crit
        weight= weight/rho_crit
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




end





;=================================================================================
;
;
;
;
;==================================================================================


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
xmin = 1.0

ymax = 1e+11
ymin = 1e+6

xaxistitle= 'R (h!E-1!N kpc)'
yaxistitle= 'M!Dgas!N (M!D!9n!3!N)'


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
radius=fload_gas_xyz('r')
mass=fload_gas_mass(1)
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








;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Circular velocity, i.e. sqrt(GM(r<R)/R)
;     ---------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------


pro plot_circular_velocity, frun, sendto, filename=filename

if not keyword_set(sendto) then sendto='ps'
if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_circular_velocity"
   print, "  "
   return
endif



setup_plot_stuff, sendto, filename=filename




; -----------------------------
; set up constants
; -----------------------------

G= 43007.1
bins = 100

; x is radius in kpc/h
xmax = 50.0
xmin = 0.0

; y is velocity in km/sec
ymax = 350
ymin = 0


radius= fltarr(bins)
baryonic_vc= fltarr(bins)
dm_vc= fltarr(bins)
total_vc= fltarr(bins)



; load variables
br=fload_baryon_xyz('r')
bmass=fload_baryon_mass(1)


dmr= fload_halo_xyz('r')
dmmass= fload_halo_mass(1)


allr= fload_all_xyz('r')
allmass= fload_all_mass(1)



binsize= (xmax-xmin)/bins


for i=1,bins do begin

	radius[i-1]= i*binsize

	; baryons
	idx= where(br le radius[i-1])
	if idx(0) ne -1 then begin
		massinside= bmass(idx)
		baryonic_vc[i-1]= sqrt(G*total(massinside)/radius[i-1])
	endif

	; dark matter
        idx= where(dmr le radius[i-1])
        if idx(0) ne -1 then begin
                massinside= dmmass(idx)
                dm_vc[i-1]= sqrt(G*total(massinside)/radius[i-1])
        endif

	; all
        idx= where(allr le radius[i-1])
        if idx(0) ne -1 then begin
                massinside= allmass(idx)
                total_vc[i-1]= sqrt(G*total(massinside)/radius[i-1])
        endif
	

endfor



ymax= max(total_vc)+10


; ---------------------------------------
; ---------------------------------------




; ------------------
; Plot this up
; ------------------
!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle="R (kpc)", $
        ytitle="V!Dc!N (km/sec)", $
        /nodata


if sendto EQ 'x' then begin

endif else begin
        oplot, radius, baryonic_vc, psym=-3, color= 150, linestyle= 1, thick=3.0
        oplot, radius, dm_vc, psym=-3, color= 50, linestyle= 3, thick=3.0
        oplot, radius, total_vc, psym=-3, linestyle=0, color= 0, thick=3.0

        x0= 0.75
        y0= 0.35 
        xyouts, x0, y0, 'Total', /normal, charthick=3.0, size= 1.2, color= 0
        xyouts, x0, y0-0.04, 'Dark Matter', /normal, charthick= 3.0, size= 1.2, color= 50
        xyouts, x0, y0-0.08, 'Baryons', /normal, charthick= 3.0, size= 1.2, color= 150
endelse


x0= 0.55
y0= 0.35
xyouts, x0, y0-0.04, fload_fid(1), /normal, size= 1.2, color= 0, charthick=2.5
xyouts, x0, y0-0.08, fload_timelbl(1,2), /normal, size= 1.2, color= 0, charthick=2.5


maxvc= max(total_vc)
maxvclbl= strcompress(string(maxvc),/remove_all)
maxvclbl= strmid(maxvclbl,0,3)
maxvclbl= 'V!Dmax!N ='+maxvclbl+' km/s'
xyouts, x0, y0, maxvclbl, /normal, charthick=3.0, size=1.2, color= 0






if (sendto EQ 'ps') then device, /close


end















