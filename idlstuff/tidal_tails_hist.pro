;-----------------------------------------------------------------------------------
;
;  Histograms of particle orbital properties
;
;
;
;
;-----------------------------------------------------------------------------------

pro toomrehist, junk


if not keyword_set(junk) then begin
        print, " "
        print, " toomrehist, junk"
        print, " "
        print, " "
        return 
endif



; -------------------
filename='thist.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


; -----------------------------
; set up constants
; -----------------------------

; quantity to histogram
;xaxistitle= "!6 Energy (arb. units)"
;xmax =  5.3
;xmin = -5.0
;xmax =  0.3
;xmin = -0.8
;xaxistitle= "!6 Eccentricity"
xaxistitle= "!6 V!D!7h!6!N / V!Dtot!N"
xmax =  1.0
xmin = -0.0
;xaxistitle= "!6 Radius (h!E-1!N kpc)"
;xmax = 25.0
;xmin = 0.0

; number (histogram)
yaxistitle= ' '
;ymax = 11.0
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






; =====================
; =====================


;load_and_print_hist, "energy", "data1/tides/disk_pro1", 0
;load_and_print_hist, "energy", "data1/tides/disk_pro1", 0

;load_and_print_hist, "energy", "data1/tides/disk_std_wpot", 0, xmin=xmin, xmax=xmax, oplotit=20
;load_and_print_hist, "energy", "data1/tides/disk_pro1_wpot", 10, xmin=xmin, xmax=xmax, oplotit= 50
;load_and_print_hist, "energy", "data1/tides/disk_pro1_wpot", 20, xmin=xmin, xmax=xmax, oplotit= 150
;load_and_print_hist, "energy", "data1/tides/disk_pro1_wpot", 30, xmin=xmin, xmax=xmax, oplotit= 100
;load_and_print_hist, "energy", "data1/tides/disk_pro1_wpot", 40, xmin=xmin, xmax=xmax, oplotit= 200

load_and_print_hist, "eccentricity", "data1/tides/disk_std_wpot", 0, xmin=xmin, xmax=xmax, oplotit=20
load_and_print_hist, "eccentricity", "data1/tides/disk_pro1_wpot", 20, xmin=xmin, xmax=xmax, oplotit= 150

;load_and_print_hist, "energy", "data1/tides/toomredisk", 10, xmin=xmin, xmax=xmax, oplotit=20
;load_and_print_hist, "energy", "data1/tides/toomredisk_dv_pro1", 10, xmin=xmin, xmax=xmax, oplotit= 150
;load_and_print_hist, "energy", "data1/tides/toomredisk_dv_pro2", 10, xmin=xmin, xmax=xmax, oplotit= 50
;load_and_print_hist, "energy", "data1/tides/toomredisk_dv_pro3", 10, xmin=xmin, xmax=xmax, oplotit= 100
;load_and_print_hist, "energy", "data1/tides/toomredisk_dv_pro4", 10, xmin=xmin, xmax=xmax, oplotit= 200


;
; toomre plots
;---------------
;load_and_print_hist, "energy", "data1/tides/toomredisk", 10, xmin=xmin, xmax=xmax, oplotit=20
;load_and_print_hist, "energy", "data1/tides/toomredisk_dv_pro4", 10, xmin=xmin, xmax=xmax, oplotit= 150
;load_and_print_hist, "eccentricity", "data1/tides/toomredisk", 10, xmin=xmin, xmax=xmax, oplotit=20
;load_and_print_hist, "eccentricity", "data1/tides/toomredisk_dv_pro4", 10, xmin=xmin, xmax=xmax, oplotit= 150
;load_and_print_hist, "radius", "data1/tides/toomredisk", 10, xmin=xmin, xmax=xmax, oplotit=20
;load_and_print_hist, "radius", "data1/tides/toomredisk_dv_pro4", 10, xmin=xmin, xmax=xmax, oplotit= 150
;load_and_print_hist, "radius", "data1/tides/toomredisk_dv_pro4", 50, xmin=xmin, xmax=xmax, oplotit= 200
				




; -------------
;  Done
; -------------

device, /close



end






;=======================================================================================






pro load_and_print_hist, quantitystring, frun, snapnum, oplotit=oplotit, $
			xmin=xmin, xmax=xmax

if not keyword_set(oplotit) then oplotit= 150

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
ecc_2 = ecc * ecc

;------------------------------------
;
; determine quantity to :
;
; energy
if quantitystring eq "energy" then begin
	quantity= energy
	mannorm= 550.0     ; good for std disk
	;quantity= energy / 2.0e5   ; good for toomredisk
	;mannorm= 87.0     ; good for toomredisk
	goto, dohist
endif


; j
if quantitystring eq "j" then begin
	quantity= fload_stars_j(1)
	goto, dohist
endif

if quantitystring eq "radius" then begin
	mannorm= 87.0
	quantity= r
	goto, dohist
endif

if quantitystring eq "eccentricity" then begin
	idx= where(ecc_2 lt 0.0)
	if idx(0) ne -1 then begin
		print, "number less than zero= ", n_elements(ecc_2)
		ecc_2(idx)= abs(ecc_2(idx))
	endif
	quantity= sqrt(ecc_2)
	mannorm= 231.0     ; good for toomredisk
	mannorm= 70426.0     ; good for std disk
	goto, dohist
endif





; 
; calc. histogram & print
;
dohist:
print, "quantity= ", quantitystring
print, "N, min, max= ", n_elements(quantity), min(quantity), max(quantity)
;mannorm= 90.0
print, "setting normalization to= ", mannorm
;temp= process_histogram(quantity, xmax=xmax, xmin=xmin, levels=11, oplotit=oplotit, mannorm=13000.0)
;temp= process_histogram(quantity, xmax=xmax, xmin=xmin, levels=11, oplotit=oplotit, normalization=normalization)
;temp= process_histogram(quantity, xmax=xmax, xmin=xmin, levels=31, oplotit=oplotit, normalization=normalization)
;temp= process_histogram(quantity, xmax=xmax, xmin=xmin, levels=1000, oplotit=oplotit, normalization=normalization)
temp= process_histogram(quantity, xmax=xmax, xmin=xmin, levels=51, oplotit=oplotit, mannorm=mannorm)
print, min(temp), max(temp)



;
; any text?
;
if keyword_set(msg) then xyouts, 0.12, 0.90, msg, size=1.2, color=oplotit, /normal, charthick=3.0


end




