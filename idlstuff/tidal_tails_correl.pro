;-----------------------------------------------------------------------------------
;
;  Histograms of particle orbital properties
;
;
;
;
;-----------------------------------------------------------------------------------

pro toomrecorrel, junk


if not keyword_set(junk) then begin
        print, " "
        print, " toomrecorrel, junk"
        print, " "
        print, " "
        return 
endif




; =====================
; =====================


frun= "data1/tides/disk_pro1_wpot"
snapnum= 20



print, "==================================================="
;ok= fload_snapshot_bh(frun,snapnum,/nopot_in_snap)
ok= fload_snapshot_bh(frun,snapnum)


G= 43007.1
Mass= 10.0    ; for toomre disk

;v=fload_disk_v('mag')
v=fload_allstars_v('tot')
print, "velocity= ", min(v), max(v)
r=fload_disk_xyz('r')
m=fload_disk_mass(1)
print, "radius= ", min(r), max(r)
print, "mass= ", min(m), max(m)


; good approx (or exact) for toomre runs
;energy= 0.5 * v * v - G * Mass / r
energy= fload_allstars_energy(1)
print, "energy= ", min(energy), max(energy)
print, "fraction with energy > 0 = ", n_elements(where(energy gt 0.0))*1.0/n_elements(energy)

;stop


;vphi=fload_allstars_v('phi')
;print, min(vphi), max(vphi)
;vphi=fload_disk_v('phi')
;print, min(vphi), max(vphi)
vtan=fload_allstars_v('tan')
print, "vtan= ", min(vtan), max(vtan)
vr=fload_allstars_v('r')
print, "vr= ", min(vr), max(vr)

ll= r * vtan

;
; really only valid for toomre disk
ecc_2= 1 + 2.*energy*ll*ll/G/G/Mass/Mass
print, "eccentricity^2= ", min(ecc_2), max(ecc_2)


;
; est. ecc for real disk
;------------------------------------
mdm= fload_halo_mass(1)
rdm= fload_halo_xyz('r')

;for i=0L, n_elements(r)-1 do begin
;	idx= where(rdm le r[i]) & if idx(0) ne -1 then mass_in_dm= total(mdm(idx)) else mass_in_dm= 0.0
;	idx= where(r le r[i]) & if idx(0) ne -1 then mass_in_s= total(m(idx)) else mass_in_s= 0.0
;	mass_in= mass_in_dm + mass_in_s
;	v_c= sqrt(G * mass_in / r[i])

	;ecc= vtan[i]/v_c
	;ecc= vtan[i]/v_c

;	if (i mod 5000) eq 0 then print, "i, r, mass_in, v_c, ecc= ", i, r[i], mass_in, v_c, ecc, "    vtan/vtot= ", vtan[i]/v[i]

;	ecc_2[i]= ecc * ecc
;endfor

ecc= vtan/v



;=================================================


; -------------------
filename='thist.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


; -----------------------------
; set up constants
; -----------------------------

; quantity to histogram
yaxistitle= "!6 Energy (arb. units)"
ymax =  5.3
ymin = -5.0

xaxistitle= "!6 V!D!7h!6!N / V!Dtot!N"
xmax =  1.0
xmin = -0.0



; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.95
y1= 0.98



;contour_makegeneralpic, radius, flux, xmax, xmin, ymax, ymin, $
contour_makegeneralpic, ecc, energy, xmax, xmin, ymax, ymin, $
                                pixels= 96, $
                                NxNImage=NxNImage

tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal



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
        ;ytickformat='(a1)', $
        ytitle=yaxistitle, $
        /nodata, /noerase




; -------------
;  Done
; -------------

device, /close



end




