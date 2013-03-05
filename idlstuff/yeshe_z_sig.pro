;------------------------------------------------------------------------
;
;    Random Procedures related to Metals
;
;
;
;
;------------------------------------------------------------------------




;--------------------------------------------------------------------------



; routine to select points
; ---------------------------
pro select_thispoint, pointselection, thispsym, thiscolor


;  open (or filled) square
; --------------------------
if pointselection eq 1 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 8
        symcolor= 150
endif

;  x's
; -----
if pointselection eq 2 then begin
        symsel= 7
        symcolor= 100
endif

;  open circle
; --------------
if pointselection eq 3 then begin
        symsize= 1.0
        ;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
        symsel= 8
        symcolor= 50
endif

;  open triangle
; ----------------
if pointselection eq 4 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
        usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
        symsel= 8
        symcolor= 0
endif

;  asterisk
; ----------
if pointselection eq 5 then begin
        symsel= 2
        symcolor= 200
endif

;  filled square
; --------------------------
if pointselection eq 6 then begin
        symsize= 1.0
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 8
        symcolor= 20
endif

;  diamond
; ----------
if pointselection eq 7 then begin
        symsel= 4
        symcolor= 120
endif

;  triangle (open)
; -----------------
if pointselection eq 8 then begin
        symsel= 5
        symcolor= 180
endif

;  plus sign
; ----------
if pointselection eq 9 then begin
        symsel= 1
        symcolor= 170
endif


thispsym=symsel
thiscolor=symcolor


end












;=====================================================================================





pro plot_met_sigma, junk, filename=filename


if not keyword_set(junk) then begin
   print, "  "
   print, " plot_met_sigma, junk, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='metsigma.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
;setup_plot_stuff, 'ps', filename=filename, colortable=2
setup_plot_stuff, 'ps', filename=filename, colortable=4


; -----------------------------
; set up constants
; -----------------------------

xmax= 600.0
xmin= 30.0
xaxistitle='!7r!6 (km s!E-1!N)'

ymax = 0.5
ymin = -1.5
yaxistitle= "!6Log <Z> (Z!D!9n!6!N)"


if not keyword_set(snapnum) then snapnum= 25

; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
        ;/ylog, $
	/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; -----------------------------------------------


; lower gas fractions
; --------------------
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3x2_e", 30, pointt= 2
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3y_e", 30, pointt= 2
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3z_e", 30, pointt= 2

;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3x2_h", 30, pointt= 2
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3y_h", 30, pointt= 2
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3z_h", 30, pointt= 2

;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3x2_k", 30, pointt= 2
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3y_k", 30, pointt= 2
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3z_k", 30, pointt= 2


; higher gas fractions
; --------------------
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3u_e", 30, pointt= 4
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3v_e", 30, pointt= 4
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3w_e", 30, pointt= 4

;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3u_h", 30, pointt= 4
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3v_h", 30, pointt= 4
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3w_h", 30, pointt= 4

;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3u_k", 30, pointt= 4
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3v_k", 30, pointt= 4
;do_one_galaxy_mets, "/raid4/tcox/gfs/vc3vc3w_k", 30, pointt= 4


; e, f=0.8, q=1.0
; ------------------
fruns= ["/raid4/tcox/bs/b0e", $
        "/raid4/tcox/bs/b1e", $
        "/raid4/tcox/bs/b2e", $
        "/raid4/tcox/bs/b3e", $
        "/raid4/tcox/bs/b4e", $
        "/raid4/tcox/bs/b5e", $
        "/raid4/tcox/bs/b6e"]

snapnums= [50, 30, 30, 30, 30, 30, 30]

xpt1= 120.0
xpt2= 200.0
xpt3= 230.0

ylbl= ymin+0.2
pointt= 2
plot_series_mets, fruns, snapnums, pointt= pointt
select_thispoint, pointt, thispsym, thiscolor
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80e_q1', size=1.2, color=thiscolor, /data, charthick=4.0




; e, f=0.4, q=1.0
; ------------------
fruns= ["/raid4/tcox/ds/d0e", $
	"/raid4/tcox/ds/d1e", $
	"/raid4/tcox/ds/d2e", $
	"/raid4/tcox/ds/d3e7", $        ; warning this is q=0.25
	"/raid4/tcox/ds/d4e", $
	"/raid4/tcox/ds/d5e", $
	"/raid4/tcox/ds/d6e"]

snapnums= [50, 30, 30, 30, 30, 30, 30]

ylbl= ymin+0.3
pointt= 1
plot_series_mets, fruns, snapnums, pointt= pointt
select_thispoint, pointt, thispsym, thiscolor
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40e_q1', size=1.2, color=thiscolor, /data, charthick=4.0




; e, f=0.4, q=0.25
; ------------------
fruns= ["/raid4/tcox/ds/d0e2_q", $
	"/raid4/tcox/ds/d1e2_q", $
	"/raid4/tcox/ds/d2e2_q", $
	"/raid4/tcox/ds/d3e7", $
	"/raid4/tcox/ds/d4e2_q", $
	"/raid4/tcox/ds/d5e2_q", $
	"/raid4/tcox/ds/d6e2_q"]

snapnums= [50, 30, 30, 30, 30, 30, 30]

ylbl= ymin+0.4
pointt= 7
plot_series_mets, fruns, snapnums, pointt= pointt
select_thispoint, pointt, thispsym, thiscolor
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40e_q25', size=1.2, color=thiscolor, /data, charthick=4.0




; h, f=0.4, q=1.0
; ----------------
;fruns= ["/raid4/tcox/ds/d0h", $
;	"/raid4/tcox/ds/d1h", $
;	"/raid4/tcox/ds/d2h", $
;	"/raid4/tcox/ds/d3h7", $
;	"/raid4/tcox/ds/d4h", $
;	"/raid4/tcox/ds/d5h", $
;	"/raid4/tcox/ds/d6h"]

;snapnums= [50, 30, 30, 30, 30, 30, 30]

;plot_series_mets, fruns, snapnums, pointt=3



; k, f=0.4, q=1.0
; -----------------
;fruns= ["/raid4/tcox/ds/d0k", $
;	"/raid4/tcox/ds/d1k", $
;	"/raid4/tcox/ds/d2k", $
;	"/raid4/tcox/ds/d3k7", $
;	"/raid4/tcox/ds/d4k", $
;	"/raid4/tcox/ds/d5k", $
;	"/raid4/tcox/ds/d6k"]

;snapnums= [50, 30, 30, 30, 30, 30, 30]

;plot_series_mets, fruns, snapnums, pointt=5



; e, f=0.2, q=0.25
; ------------------
fruns= ["/raid4/tcox/es/e0e", $
	"/raid4/tcox/es/e1e", $
	"/raid4/tcox/es/e2e", $
	"/raid4/tcox/es/e3e", $
	"/raid4/tcox/es/e4e", $
	"/raid4/tcox/es/e5e", $
	"/raid4/tcox/es/e6e"]

snapnums= [50, 30, 30, 30, 30, 30, 30]

ylbl= ymin+0.5
pointt= 3
plot_series_mets, fruns, snapnums, pointt= pointt
select_thispoint, pointt, thispsym, thiscolor
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '20e_q25', size=1.2, color=thiscolor, /data, charthick=4.0



; e, f=0.05, q=0.25
; ------------------
fruns= ["/raid4/tcox/zs/z0e", $
	"/raid4/tcox/zs/z1e", $
	"/raid4/tcox/zs/z2e", $
	"/raid4/tcox/zs/z3e", $
	"/raid4/tcox/zs/z4e", $
	"/raid4/tcox/zs/z5e", $
	"/raid4/tcox/zs/z6e"]

snapnums= [50, 30, 30, 30, 30, 30, 30]


ylbl= ymin+0.6
pointt= 5
plot_series_mets, fruns, snapnums, pointt= pointt
select_thispoint, pointt, thispsym, thiscolor
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '5e_q25', size=1.2, color=thiscolor, /data, charthick=4.0






; print extras
; -------------

;xyouts, 0.55, 0.90, 'no black hole', size=1.5, color=0, /normal, charthick=4.0
;xyouts, 0.55, 0.85, '!6time course', size=1.5, color=0, /normal, charthick=4.0
;xyouts, 0.55, 0.80, 'dark -> ', size=1.5, color=20, /normal, charthick=4.0
;xyouts, 0.75, 0.80, ' light', size=1.5, color=160, /normal, charthick=4.0

;xyouts, 0.35, 0.27, 'B-luminosity weighted', size=1.5, color=0, /normal, charthick=4.0

;xyouts, 0.45, 0.35, 'New Stars', size=1.5, color=0, /normal, charthick=4.0
;xyouts, 0.45, 0.31, 'mass weighted (same as)', size=1.5, color=0, /normal, charthick=4.0
;xyouts, 0.45, 0.27, 'unweighted average', size=1.5, color=0, /normal, charthick=4.0

;xyouts, 0.35, 0.35, 'New Stars', size=1.5, color=0, /normal, charthick=4.0
;xyouts, 0.35, 0.31, 'B-luminosity weighted', size=1.5, color=0, /normal, charthick=4.0



; done
; -----
device, /close


end
   




;--------------------------------------------------------------------------




;  Do a series of runs
; -----------------------------------------------
pro plot_series_mets, fruns, snapnums, pointt=pointt


ns= n_elements(fruns)
ssigma= fltarr(ns)
savg_z= fltarr(ns)

for i=0, ns-1 do begin

	savg_z[i]= do_one_galaxy_mets(fruns[i], snapnums[i])

	;---------------------------

	read_sigma_file, fruns[i], time, Asigxy, Asigxz, Asigyz, Asigavg, Asigerr, $
                                Bsigxy, Bsigxz, Bsigyz, Bsigavg, Bsigerr, $
                                gsigavg, gsigerr

	sigmaidx= n_elements(time)-1
	ssigma[i]= Asigavg[sigmaidx]

	;---------------------------

endfor
    
select_thispoint, pointt, thispsym, thiscolor
    
; make it log
idx= where(savg_z le 0.0)
if idx(0) ne -1 then savg_z(idx)= 1.0e-6
savg_z= alog10(savg_z)

oplot, ssigma, savg_z, psym=-thispsym, color=thiscolor, thick=3.0

print, "-----------------------------------------------"



end







;--------------------------------------------------------------------------





;  Do the grunt work
; -----------------------------------------------
function do_one_galaxy_mets, frun, snapnum


ok= fload_snapshot_bh(frun,snapnum)


;---------------------------
nsz= fload_newstars_z(1)
;nsm= fload_newstars_mass(1)
avg_z= mean(nsz)/0.02
print, "average Z= ",avg_z


;---------------------------
;do_lum_weighted= 1
do_lum_weighted= 0
if do_lum_weighted eq 1 then begin
    TTime= float(fload_time(1))
TTime= 10.0

    ndisk= fload_npart(2)                               ; disk particles
    nbulge= fload_npart(3)                              ; bulge particles
    nstars= fload_npart(4)
    npart= long(ndisk) + long(nbulge) + long(nstars)
    print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars
    N= npart
    m= 1.0e+10*fload_allstars_mass(1)
    age=fload_allstars_age(1)
    age=float(TTime-age)
    zmets=fload_allstars_z(1)


    ; get the luminosities
    ;  - in units of solar luminosities
    print, "load luminosities"
    load_all_stellar_luminosities, N, TTime, m, age, zmets, $
                                Lum_Bol, Lum_U, Lum_B, Lum_V, Lum_R, Lum_I, Lum_J, Lum_H, Lum_K, $
                                Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z   ;, $
                                ;/notsolar

    ; trap for any NaN's
    idx=where(finite(Lum_B) eq 0)
    ;idx=where(finite(Lum_K) eq 0) 
    if idx(0) ne -1 then begin
        Lum_Bol(idx)= 100.0 & Lum_U(idx)= 100.0 & Lum_B(idx)= 100.0 & Lum_V(idx)= 100.0 & Lum_R(idx)= 100.0
        Lum_I(idx)= 100.0 & Lum_J(idx)= 100.0 & Lum_H(idx)= 100.0 & Lum_K(idx)= 100.0
        Lum_sdss_u(idx)= 100.0 & Lum_sdss_g(idx)= 100.0 & Lum_sdss_r(idx)= 100.0
        Lum_sdss_i(idx)= 100.0 & Lum_sdss_z(idx)= 100.0
        print, "Fixed ",n_elements(idx), " luminosities when trapping for NaN values."
    endif
    
    
    ;print, "Total U-band Lum= ", total(Lum_U)
    ;print, "Total B-band Lum= ", total(Lum_B)
    ;print, "Total V-band Lum= ", total(Lum_V)
    ;print, "Total K-band Lum= ", total(Lum_K)
    ;print, "Total sdss_u-band Lum= ", total(Lum_sdss_u)
    
; new star's
;
avg_z= mean(nsz)/0.02
print, "New Star average Z= ",avg_z

;nsMass= fload_newstars_mass(1)

;avg_z= total(nsz*nsMass)/total(nsMass)/0.02
;print, "New Star mass-weighted average Z= ",avg_z

nsLum_B= Lum_B[ndisk+nbulge:ndisk+nbulge+nstars-1]
avg_z_lum= total(nsz*nsLum_B)/total(nsLum_B)/0.02
print, "New Star luminosity-weighted average Z= ",avg_z_lum
print, "New Star Lum_B min/max", min(nsLum_B), max(nsLum_B)
;avg_z= avg_z_lum


  ;do_allstars= 1
  do_allstars= 0
  if do_allstars eq 1 then begin
    ; luminosity weighted
    avg_z_lumwt= total((zmets/0.02)*Lum_B)/total(Lum_B)
    print, "Luminosity-weighted average Z=", avg_z_lumwt

    avg_z_lumwt2= total(zmets*Lum_B)/total(Lum_B)/0.02
    print, "Luminosity-weighted average Z #2=", avg_z_lumwt2

    Mass= fload_allstars_mass(1)
    avg_z_masswt= total(zmets*Mass)/total(Mass)/0.02
    print, "Mass-weighted average Z=", avg_z_masswt

    ; no weighting
    avg_z= mean(zmets)/0.02
    print, "Average Z=", avg_z
    print, "All Star Lum_B min/max", min(Lum_B), max(Lum_B)

    ;avg_z= avg_z_lumwt
    ;avg_z= avg_z_masswt

  endif

endif



return, avg_z


end




;--------------------------------------------------------------------------



