;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   3D Gas Profiles
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------





;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     This is our 3D figure for the *revised* 
;     X-ray letter.
;     -------------------------------------------
;  
;  
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------




; --------------------------
;  The whole thing
; --------------------------
pro plot_whole_thing, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "plot_whole_think, junk"
   print, "  "
   return                       
endif                           

filename='gasprofiles.eps'
                
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=16, newysize=20



                
; ----------------------------- 
; set up constants      
; -----------------------------


frun1="/raid4/tcox/vc3vc3e_2"
frun2="/raid4/tcox/vc3vc3e_no"
snapnum= 107
                        
bins = 30
                
;xmax = 3
;xmin = -1
xmax = 300
xmin = 1.0

xaxistitle= 'R (kpc)'




;---------------------------
;  Print it
;---------------------------
;
; Produce a multi-figure
; such as below.
;
; --------   --------
; |      |   |      |
; |  #1  |   |  #2  |
; |      |   |      |
; --------   --------
; |      |   |      |
; |  #3  |   |  #4  |
; |      |   |      |
; --------   --------
; |      |   |      |
; |  #5  |   |  #6  |
; |      |   |      |
; --------   --------
;
;

x0= 0.12 & x1= 0.495 & x2= 0.615 & x3= 0.99

y0= 0.09 & y1= 0.39 & y2= 0.69 & y3= 0.99



; plot #1
; -------
; cumulative total mass
;
yaxistitle= 'M!Dtotal!N (M!D!9n!3!N)'
ymax = 4.0e+12
ymin = 7.0e+9
generate_axes, x0, y2, x1, y3, xmax, xmin, ymax, ymin, yaxistitle, /yexplabel

processandplot_onemassprofile, frun1, snapnum, xmin, xmax, bins, linecolor= 150, /totalmass

processandplot_onemassprofile, frun2, snapnum, xmin, xmax, bins, linecolor= 50, /totalmass


; plot #2
; -------
; gas mass
yaxistitle= 'M!Dgas!N (M!D!9n!3!N)'
ymax = 2.0e+10
ymin = 7.0e+5
generate_axes, x2, y2, x3, y3, xmax, xmin, ymax, ymin, yaxistitle, /yexplabel

processandplot_onemassprofile, frun1, snapnum, xmin, xmax, bins, linecolor= 150

processandplot_onemassprofile, frun2, snapnum, xmin, xmax, bins, linecolor= 50




; plot #3
; -------
; gas/total mass
;
yaxistitle= 'M!Dgas!N/M!Dtotal!N'
ymax = 0.010
ymin = 0.0
generate_axes, x0, y1, x1, y2, xmax, xmin, ymax, ymin, yaxistitle, /linear

processandplot_onemassprofile, frun1, snapnum, xmin, xmax, bins, linecolor= 150, /gas_to_tot

processandplot_onemassprofile, frun2, snapnum, xmin, xmax, bins, linecolor= 50, /gas_to_tot




; plot #4
; -------
; Entropy
;
yaxistitle= 'Entropy (keV cm!E-2!N) '
ymax = 700.0
ymin = 2.0
generate_axes, x2, y1, x3, y2, xmax, xmin, ymax, ymin, yaxistitle

processandplot_oneentropyprofile, frun1, snapnum, xmin, xmax, bins, linecolor= 150

processandplot_oneentropyprofile, frun2, snapnum, xmin, xmax, bins, linecolor= 50




; plot #5
; -------
; cooling time
;
yaxistitle= 't!Dcool!N (Gyr) '
ymax = 80.0
ymin = 0.1
generate_axes, x0, y0, x1, y1, xmax, xmin, ymax, ymin, yaxistitle, xaxistitle

processandplot_onetcoolprofile, frun1, snapnum, xmin, xmax, bins, linecolor= 150

processandplot_onetcoolprofile, frun2, snapnum, xmin, xmax, bins, linecolor= 50



; plot #6
; -------
; M/L
;
yaxistitle= 'M!Dtot!N/L!DB!N (M!D!9n!3!NL!D!9n!3!N!E-1!N) '
ymax = 70.0
ymin = 1.0
generate_axes, x2, y0, x3, y1, xmax, xmin, ymax, ymin, yaxistitle, xaxistitle

processandplot_onemlprofile, frun1, snapnum, xmin, xmax, bins, linecolor= 150

processandplot_onemlprofile, frun2, snapnum, xmin, xmax, bins, linecolor= 50




; -----------------
;  Plot Extras
; -----------------


xyouts, x0+0.03, y3-0.04, 'black hole', /normal, charthick=3.5, size= 1.4, color= 150
xyouts, x0+0.03, y3-0.07, 'no black hole', /normal, charthick=3.5, size= 1.4, color= 50


;--------------------------------------
;--------------------------------------

device, /close



end




; ------------------------------------------------------------------




; --------------------------
;  The whole thing
; --------------------------
pro plot_3_thing, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "plot_3_think, junk"
   print, "  "
   return                       
endif                           

filename='gasprofiles.eps'
                
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=10, newysize=20



                
; ----------------------------- 
; set up constants      
; -----------------------------


;frun1="/raid4/tcox/vc3vc3e_2"
;frun2="/raid4/tcox/vc3vc3e_no"
;snapnum1= 107
;snapnum2= snapnum1
frun1="/raid4/tcox/vc3vc3h_2"
frun2="/raid4/tcox/vc3vc3h_no"
snapnum1= 107
snapnum2= 30
;frun1="/raid4/tcox/vc3bvc3b"
;frun2="/raid4/tcox/vc3bvc3b_no"
;snapnum= 30
                        
bins = 20
                
;xmax = 3
;xmin = -1
xmax = 300
xmin = 1.0

xaxistitle= '!6radius (kpc)'




;---------------------------
;  Print it
;---------------------------
;
; Produce a multi-figure
; such as below.
;
; --------
; |      |
; |  #1  |
; |      |
; --------
; |      |
; |  #2  |
; |      |
; --------
; |      |
; |  #3  |
; |      |
; --------
;
;

x0= 0.19 & x1= 0.98

y0= 0.09 & y1= 0.39 & y2= 0.69 & y3= 0.99



; plot #1
; -------
; gas density
;
yaxistitle= '!7q!6 (cm!E-3!N)'
ymax = 2.0e-1
ymin = 2.0e-7
generate_axes, x0, y2, x1, y3, xmax, xmin, ymax, ymin, yaxistitle


; these are for procedure in processandplot_onedensityprofile
;processandplot_onedensityprofile, frun1, snapnum, xmin, xmax, bins, linecolor= 50, /density_in_cm3, /xrwt
;processandplot_onedensityprofile, frun2, snapnum, xmin, xmax, bins, linecolor= 150, /density_in_cm3, /xrwt

; these are for procedure in here
processandplot_onedensityprofile, frun1, snapnum1, xmin, xmax, bins, linecolor= 50, /xrwt
processandplot_onedensityprofile, frun2, snapnum2, xmin, xmax, bins, linecolor= 150, /xrwt


; want these on the top graph, not the bottom
symsize= 1.2
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
oplot, [30.8, 42.7], [4.4e-2,4.4e-2], thick=3.0, psym=-8, color= 50, linestyle=2
;xyouts, 0.67, 0.93, 'black hole (RS)', /normal, charthick=3.0, size=1.33, color=50
xyouts, 0.75, 0.95, 'BH (RS)', /normal, charthick=3.0, size=1.33, color=50


symsize= 1.0
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;oplot, [18.0, 24.4], [0.405,0.405], thick=3.0, psym=-2, color= 150
oplot, [30.8, 42.7], [1.3e-2,1.3e-2], thick=3.0, psym=-8, color= 150
;xyouts, 0.67, 0.90, 'no black hole', /normal, charthick=3.0, size=1.33, color=150
xyouts, 0.75, 0.925, 'std (RS)', /normal, charthick=3.0, size=1.33, color=150



; plot #2
; -------
; cooling time
;
yaxistitle= '!6t!Dcool!N (Gyr) '
ymax = 70.0
ymin = 0.03
generate_axes, x0, y1, x1, y2, xmax, xmin, ymax, ymin, yaxistitle

processandplot_onetcoolprofile, frun1, snapnum1, xmin, xmax, bins, linecolor= 50, /xrwt, Metallicity='-05'

processandplot_onetcoolprofile, frun2, snapnum2, xmin, xmax, bins, linecolor= 150, /xrwt, Metallicity='-15'



; plot #3
; -------
; Entropy
;
yaxistitle= '!6Entropy (keV cm!E2!N) '
ymax = 500.0
ymin = 0.7
generate_axes, x0, y0, x1, y1, xmax, xmin, ymax, ymin, yaxistitle, xaxistitle

processandplot_oneentropyprofile, frun1, snapnum1, xmin, xmax, bins, linecolor= 50, /xrwt

processandplot_oneentropyprofile, frun2, snapnum2, xmin, xmax, bins, linecolor= 150, /xrwt






; -----------------
;  Plot Extras
; -----------------



;--------------------------------------
;--------------------------------------

device, /close



end




; ------------------------------------------------------------------




; ----------------------
; generate axes and such
; ----------------------
pro generate_axes, x0, y0, x1, y1, $
		xmax, xmin, ymax, ymin, $
		yaxistitle, xaxistitle, $
		linear=linear, $
		yexplabel=yexplabel

!p.position= [x0, y0, x1, y1]

if keyword_set(linear) then begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
	return
endif

if keyword_set(yexplabel) then begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, /ylog, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytickformat='exp_label', ytitle=yaxistitle, /nodata, /noerase
	return
endif

if keyword_set(xaxistitle) then begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, /xlog, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase
endif else begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, /xlog, $
	xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=4.0, $
	xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
endelse


end










;==================================================================================


; --------------------------
;
; Actually do the work to
; determine the cumulative
; mass profile.
;
;
; 
; --------------------------
pro processandplot_onemassprofile, frun, snapnum, xmin, xmax, bins, $
					linecolor=linecolor, $
					totalmass=totalmass, $
					gas_to_tot=gas_to_tot


ok=fload_snapshot_bh(frun,snapnum)

; gas 
radius=fload_gas_xyz('r')
mass=fload_gas_mass(1)
print, "gas mass= ", total(mass)
mass= 1.0d+10 * mass

; everything
r_all= fload_all_xyz('r')
m_all= fload_all_mass(1)
print, "total mass= ", total(m_all)
m_all= 1.0d+10 * m_all

; bookkeeping
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
    rs[i]= currentr

    idx= where(radius le currentr)
    gasmass_within_r= total(mass(idx))

    idx= where(r_all le currentr)
    totalmass_within_r= total(m_all(idx))

    cumulative_mass[i]= gasmass_within_r
    if keyword_set(totalmass) then cumulative_mass[i]= totalmass_within_r
    if keyword_set(gas_to_tot) then cumulative_mass[i]= gasmass_within_r/totalmass_within_r
endfor

oplot, rs, cumulative_mass, psym=-3, linestyle=0, color= linecolor, thick=4.0


end




;==================================================================================




; --------------------------
;
; Actually do the work to
; determine the entropy
; profile.
;
;
; 
; --------------------------
pro processandplot_oneentropyprofile, frun, snapnum, xmin, xmax, bins, $
					linecolor=linecolor, $
					xrwt=xrwt


ok=fload_snapshot_bh(frun,snapnum)

; gas 
radius=fload_gas_xyz('r')
radius= radius/0.7
entropy=fload_gas_entropy(1) / (0.7^(4./3.))

; only x-ray gas
if keyword_set(xrwt) then begin
        xray=fload_gas_xray_luminosity(1,/diffuse_hotgas)
        idx= where(xray gt 0.0)
        if idx(0) ne -1 then begin
                radius= radius(idx)
                entropy= entropy(idx)
        endif else begin
                print, "No HOT particles!"
                return
        endelse
endif


; bookkeeping
rs = fltarr(bins)
ent_pro= fltarr(bins)

tempxmin= alog10(xmin)
tempxmax= alog10(xmax)

;--------------------
; process the cumulative profile

binsize = float((tempxmax-tempxmin))/bins

for i=1, bins do begin
    sm_r= tempxmin + (i-1)*binsize
    lg_r= tempxmin + (i)*binsize
    currentr= 10^(0.5*(sm_r+lg_r))
    rs[i-1]= currentr

    idx= where((radius ge 10^sm_r) and (radius le 10^lg_r))

    if idx(0) ne -1 then ent_pro[i-1]= mean(entropy(idx))
endfor

idx= where(ent_pro gt 0.0)
if idx(0) ne -1 then begin
	ent_pro= ent_pro(idx)
	rs= rs(idx)
endif

if linecolor eq 150 then begin
        symsize= 1.0
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        oplot, rs, ent_pro, psym=-8, linestyle=0, color= linecolor, thick=3.0
endif

if linecolor eq 50 then begin
        symsize= 1.2
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=3.0
        oplot, rs, ent_pro, psym=-8, linestyle=2, color= linecolor, thick=3.0
endif


end







;=================================================================================



; --------------------------
;
; Actually do the work to
; determine the cooling time.
;
;
; 
; --------------------------
pro processandplot_onetcoolprofile, frun, snapnum, xmin, xmax, bins, $
                                        linecolor=linecolor, $
					xrwt=xrwt, Metallicity=Metallicity


ok=fload_snapshot_bh(frun,snapnum)

; gas 
radius=fload_gas_xyz('r')
radius=radius/0.7
cooltime= fload_gas_coolingtime(Metallicity) / 0.7 / 0.7

; only x-ray gas
if keyword_set(xrwt) then begin
        xray=fload_gas_xray_luminosity(1,/diffuse_hotgas)
        idx= where(xray gt 0.0)
        if idx(0) ne -1 then begin
                radius= radius(idx)
                cooltime= cooltime(idx)
        endif else begin
                print, "No HOT particles!"
                return
        endelse
endif

; bookkeeping
rs = fltarr(bins)
ct_pro= fltarr(bins)

tempxmin= alog10(xmin)
tempxmax= alog10(xmax)

;--------------------
; process the cumulative profile

binsize = float((tempxmax-tempxmin))/bins

for i=1, bins do begin
    sm_r= tempxmin + (i-1)*binsize
    lg_r= tempxmin + (i)*binsize
    currentr= 10^(0.5*(sm_r+lg_r))
    rs[i-1]= currentr

    idx= where((radius ge 10^sm_r) and (radius le 10^lg_r))

    if idx(0) ne -1 then ct_pro[i-1]= mean(cooltime(idx))
endfor

idx= where(ct_pro gt 0.0)
if idx(0) ne -1 then begin
        ct_pro= ct_pro(idx)
        rs= rs(idx)
endif


if linecolor eq 150 then begin
        symsize= 1.0
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        oplot, rs, ct_pro, psym=-8, linestyle=0, color= linecolor, thick=3.0
endif

if linecolor eq 50 then begin
        symsize= 1.2
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=3.0
        oplot, rs, ct_pro, psym=-8, linestyle=2, color= linecolor, thick=3.0
endif



end






;=================================================================================



; --------------------------
;
; Actually do the work to
; determine the radial
; mass-to-light ratio.
;
; 
; --------------------------
pro processandplot_onemlprofile, frun, snapnum, xmin, xmax, bins, $
                                        linecolor=linecolor


ok=fload_snapshot_bh(frun,snapnum)

; stars - light
TTime= float(fload_time(1))

radius= fload_allstars_xyz('r')

N= fload_npart(2)+fload_npart(3)+fload_npart(4)
m= 1.0e+10*fload_allstars_mass(1)
age=fload_allstars_age(1)
age=float(TTime-age)
zmets=fload_allstars_z(1)


; get the luminosities
print, "load luminosities"
load_all_stellar_luminosities, N, TTime, m, age, zmets, $
        Lum_Bol, Lum_U, Lum_B, Lum_V, Lum_R, Lum_I, Lum_J, Lum_H, Lum_K, $
        Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z
        ;/notsolar  - yes, we'll take it in solar now

; trap for any NaN's
idx=where(finite(Lum_B) eq 0)
if idx(0) ne -1 then begin
        Lum_Bol(idx)= 100.0 & Lum_U(idx)= 100.0 & Lum_B(idx)= 100.0 & Lum_V(idx)= 100.0 & Lum_R(idx)= 100.0
        Lum_I(idx)= 100.0 & Lum_J(idx)= 100.0 & Lum_H(idx)= 100.0 & Lum_K(idx)= 100.0
        Lum_sdss_u(idx)= 100.0 & Lum_sdss_g(idx)= 100.0 & Lum_sdss_r(idx)= 100.0
        Lum_sdss_i(idx)= 100.0 & Lum_sdss_z(idx)= 100.0
        print, "Fixed ",n_elements(idx), " luminosities when trapping for NaN values."
endif
 

; everything - mass
r_all= fload_all_xyz('r')
m_all= fload_all_mass(1)
print, "total mass= ", total(m_all)
m_all= 1.0d+10 * m_all


; bookkeeping
rs = fltarr(bins)
ml_pro= fltarr(bins)

tempxmin= alog10(xmin)
tempxmax= alog10(xmax)

;--------------------
; process the cumulative profile

binsize = float((tempxmax-tempxmin))/bins

for i=1, bins do begin
    lg_r= tempxmin + (i)*binsize
    currentr= 10^(lg_r)
    rs[i-1]= currentr

    idx= where(radius le currentr)
    if idx(0) ne -1 then luminosity= total(Lum_B(idx)) else luminosity= 1.0

    idx= where(r_all le currentr)
    if idx(0) ne -1 then mass = total(m_all(idx)) else mass= 0.0

    ml_pro[i-1]= mass/luminosity
endfor

idx= where(ml_pro gt 0.0)
if idx(0) ne -1 then begin
        ml_pro= ml_pro(idx)
        rs= rs(idx)
endif


oplot, rs, ml_pro, psym=-3, linestyle=0, color= linecolor, thick=4.0


end






; --------------------------
;
; Actually do the work to
; determine the gas density
; profile.
;
;
; 
; --------------------------
pro processandplot_onedensityprofile, frun, snapnum, xmin, xmax, bins, $
					linecolor=linecolor, $
					xrwt=xrwt


ok=fload_snapshot_bh(frun,snapnum)

; gas 
radius=fload_gas_xyz('r')
radius= radius/0.7
rho=0.7*0.7*fload_gas_rho(1) * 674.59
;rho=fload_gas_rho(1) * 404.59

; only x-ray gas
if keyword_set(xrwt) then begin
        xray=fload_gas_xray_luminosity(1,/diffuse_hotgas)
        idx= where(xray gt 0.0)
        if idx(0) ne -1 then begin
                radius= radius(idx)
                rho= rho(idx)
        endif else begin
                print, "No HOT particles!"
                return
        endelse
endif


; bookkeeping
rs = fltarr(bins)
rho_pro= fltarr(bins)

tempxmin= alog10(xmin)
tempxmax= alog10(xmax)

;--------------------
; process the cumulative profile

binsize = float((tempxmax-tempxmin))/bins

for i=1, bins do begin
    sm_r= tempxmin + (i-1)*binsize
    lg_r= tempxmin + (i)*binsize
    currentr= 10^(0.5*(sm_r+lg_r))
    rs[i-1]= currentr

    idx= where((radius ge 10^sm_r) and (radius le 10^lg_r))

    if idx(0) ne -1 then rho_pro[i-1]= mean(rho(idx))
endfor

idx= where(rho_pro gt 0.0)
if idx(0) ne -1 then begin
	rho_pro= rho_pro(idx)
	rs= rs(idx)
endif

if linecolor eq 150 then begin
        symsize= 1.0
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        oplot, rs, rho_pro, psym=-8, linestyle=0, color= linecolor, thick=3.0
endif

if linecolor eq 50 then begin
        symsize= 1.2
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=3.0
        oplot, rs, rho_pro, psym=-8, linestyle=2, color= linecolor, thick=3.0
endif


end









;=========================================================================









; -----------------------------
;  The spherical gas density
; -----------------------------
pro gas_density, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "gas_density, junk"
   print, "  "
   return                       
endif                           

;filename='gasprofile.eps'
filename='gasprofile2.eps'
                
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=14, newysize=14



                
; ----------------------------- 
; set up constants      
; -----------------------------


frun1="/raid4/tcox/vc3bvc3b"
frun2="/raid4/tcox/vc3bvc3b_no"
snapnum= 30
                        
bins = 25
                
;xmax = 3
;xmin = -1
xmax = 300
xmin = 1.0

xaxistitle= 'r (kpc)'




;---------------------------
;  Print it
;---------------------------
;

x0= 0.15
x1= 0.99
y0= 0.15
y1= 0.99



; plot it
; -------
; gas density
;
;yaxistitle= '!7q!3(10 M!D!9n!3!N kpc!E-3!N)'
;ymax = 1.0e-2
;ymin = 1.0e-9
yaxistitle= '!7q!3(cm!E-3!N)'
ymax = 3.0e-2
ymin = 1.0e-7


!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, /xlog, $
xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
;xtickformat='(a1)', ytickformat='exp_label', $
xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase


; use existing procedures in 
;   determine_profiles.pro
;processandplot_onedensityprofile, frun1, snapnum, xmin, xmax, bins, linecolor= 150
;processandplot_onedensityprofile, frun1, snapnum, xmin, xmax, bins, linecolor= 150, /hotgasonly
;processandplot_onedensityprofile, frun1, snapnum, xmin, xmax, bins, linecolor= 50, /density_in_cm3
processandplot_onedensityprofile, frun1, snapnum, xmin, xmax, bins, linecolor= 50, /density_in_cm3, /xrwt

;processandplot_onedensityprofile, frun2, snapnum, xmin, xmax, bins, linecolor= 50
;processandplot_onedensityprofile, frun2, snapnum, xmin, xmax, bins, linecolor= 50, /hotgasonly
;processandplot_onedensityprofile, frun2, snapnum, xmin, xmax, bins, linecolor= 150, /density_in_cm3
processandplot_onedensityprofile, frun2, snapnum, xmin, xmax, bins, linecolor= 150, /density_in_cm3, /xrwt



; -----------------
;  Plot Extras
; -----------------


; want these on the top graph, not the bottom
symsize= 1.2
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
oplot, [2.8, 3.7], [5.0e-6,5.0e-6], thick=3.0, psym=-8, color= 50, linestyle=2
;xyouts, 0.67, 0.93, 'black hole (RS)', /normal, charthick=3.0, size=1.33, color=50
xyouts, 0.37, 0.40, 'BH (RS)', /normal, charthick=3.0, size=1.33, color=50


symsize= 1.0
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;oplot, [18.0, 24.4], [0.405,0.405], thick=3.0, psym=-2, color= 150
oplot, [2.8, 3.7], [2.5e-6,2.5e-6], thick=3.0, psym=-8, color= 150
;xyouts, 0.67, 0.90, 'no black hole', /normal, charthick=3.0, size=1.33, color=150
xyouts, 0.37, 0.36, 'std (RS)', /normal, charthick=3.0, size=1.33, color=150



;--------------------------------------
;--------------------------------------

device, /close



end

