;------------------------------------------------------------------------
;
;    Random Procedures related to Metals
;
;
;
;
;------------------------------------------------------------------------










;------------------------------------------------------------------------
;
;      Metal Fractions versus Mass
;     ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro metal_v_size, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "metal_v_size, junk"
   print, "  "
   return
endif

filename='metals_stars.eps'
;filename='metals_ubnd.eps'
;filename='metals_xray.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4 


;--------------------------------------
;--------------------------------------

xaxistitle = "Total Mass (10!E10!N M!D!9n!3!N)"
xmax= 10000.0
xmin= 1.0

;yaxistitle = "Unbound Metal Fraction (Z!D!9n!3!N)"
;yaxistitle = "Unbound Metal Fraction"
yaxistitle = "New Star Metal Fraction"
;yaxistitle = "X-ray Gas Metal Fraction"
ymax = 1.0
ymin = 0.0
;ymin = 0.1
;ymin = 0.01



;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;/ylog, $
	/xlog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; snapnum is the LAST snapshot

;process_and_plot_onegalaxy, "/raid4/tcox/vc3vc3e", 30, pointt= 2
;process_and_plot_onegalaxy, "/raid4/tcox/vc3vc3e_no", 107, pointt= 2

;process_and_plot_onegalaxy, "/raid4/tcox/vc3vc3e_gf6", 107, pointt= 3
;process_and_plot_onegalaxy, "/raid4/tcox/vc3vc3e_gf8", 107, pointt= 3
;process_and_plot_onegalaxy, "/raid4/tcox/vc3vc3e_gf8a", 107, pointt= 3

;process_and_plot_onegalaxy, "/raid4/tcox/vc1vc1", 30, pointt= 4
;process_and_plot_onegalaxy, "/raid4/tcox/vc2vc2", 30, pointt= 4


;process_and_plot_onegalaxy, "/raid4/tcox/vc4vc4a", 30, pointt= 6
;process_and_plot_onegalaxy, "/raid4/tcox/vc5vc5a", 30, pointt= 6
;process_and_plot_onegalaxy, "/raid4/tcox/vc6vc6a", 30, pointt= 6

;process_and_plot_onegalaxy, "/raid4/tcox/vc4vc4_no", 30, pointt= 7
;process_and_plot_onegalaxy, "/raid4/tcox/vc5vc5_no", 30, pointt= 7
;process_and_plot_onegalaxy, "/raid4/tcox/vc6vc6_no", 30, pointt= 7


;process_and_plot_onegalaxy, "/raid4/tcox/As/A1", 40, pointt= 5
;process_and_plot_onegalaxy, "/raid4/tcox/As/A2", 40, pointt= 5
;process_and_plot_onegalaxy, "/raid4/tcox/As/A3", 154, pointt= 5
;process_and_plot_onegalaxy, "/raid4/tcox/As/A4", 40, pointt= 5
;process_and_plot_onegalaxy, "/raid4/tcox/As/A5", 40, pointt= 5

process_and_plot_onegalaxy, "/raid4/tcox/vc3vc3e", 22, pointt= 2
process_and_plot_onegalaxy, "/raid4/tcox/vc3vc3h", 22, pointt= 2

process_and_plot_onegalaxy, "/raid2/tcox/vc3HGe", 22, pointt= 5
process_and_plot_onegalaxy, "/raid2/tcox/vc3HGh", 22, pointt= 5

device, /close



end




;--------------------------------------------------------------------------




;  Do the grunt work
; -----------------------------------------------
pro process_and_plot_onegalaxy, frun, snapnum, pointt=pointt



; -------------------------------------------------
; need to load id_accounting
;   (and order by id number)
idlistname= frun+"/hotgas_idlist.txt"
idlist_fmfile= load_id_list(idlistname)
print, "id's in file= ",n_elements(idlist_fmfile)
sidx=sort(idlist_fmfile)
idlist_fmfile= idlist_fmfile(sidx)
; -------------------------------------------------
; need to load id_accounting
;   (and order by id number)
idlistname= frun+"/wind_idlist.txt"
idlist_fmfile_2= load_id_list(idlistname)
print, "id's in file= ",n_elements(idlist_fmfile_2)
sidx=sort(idlist_fmfile_2)
idlist_fmfile_2= idlist_fmfile_2(sidx)
; -------------------------------------------------


ok= fload_snapshot_bh(frun,snapnum)

dm_mass= fload_halo_mass(1)
d_mass= fload_disk_mass(1)
b_mass= fload_bulge_mass(1)

gasz= fload_gas_metallicity(1)
gasm= fload_gas_mass(1)
nsz= fload_newstars_z(1)
nsm= fload_newstars_mass(1)
gids= fload_gas_id(1)


; total metals in stars
metals_in_stars= total(nsz*nsm)


; total metals in gas
metals_in_gas= total(gasz*gasm)


; total metals
total_metals= metals_in_stars+metals_in_gas


        ; Hot X-ray Gas
        ; ---------------------
        ; 
        ; select desired id's 
        idx= intarr(n_elements(gids))
        for ii=0L,n_elements(idlist_fmfile)-1 do begin
                ;if idlist_fmfile[i] eq lstid then duplicates= duplicates+1
                inlst= where(gids eq idlist_fmfile[ii])
                if inlst(0) ne -1 then begin
                   idx(inlst)= 1
                endif else begin
                   print, "couldn't find id=",idlist_fmfile[ii]
                endelse
                ;lstid= idlist_fmfile[i]
        endfor
        midx= where(idx eq 1)
        print, "matching gid's= ",n_elements(midx)

metals_in_xraygas= total(gasz(midx)*gasm(midx))



        ; Wind (E>=0) Gas
        ; ---------------------
        ; 
        ; select desired id's 
        idx= intarr(n_elements(gids))
        for ii=0L,n_elements(idlist_fmfile_2)-1 do begin
                inlst= where(gids eq idlist_fmfile_2[ii])
                if inlst(0) ne -1 then begin
                   idx(inlst)= 1
                endif else begin
                   print, "couldn't find id=",idlist_fmfile_2[ii]
                endelse
        endfor
        midx= where(idx eq 1)
        print, "matching gid's= ",n_elements(midx)

metals_in_unboundgas= total(gasz(midx)*gasm(midx))


totalmass= total(dm_mass)+total(d_mass)+total(b_mass)+total(gasm)+total(nsm)


select_thispoint, pointt, thispsym, thiscolor


print, "---------------------------------"
print, "total mass = ", totalmass
print, "total metals = ", total_metals
print, "newstar f= ", metals_in_stars/total_metals
print, "unbound f= ", metals_in_unboundgas/total_metals
print, "xray    f= ", metals_in_xraygas/total_metals
print, "---------------------------------"



;oplot, [totalmass], [metals_in_unboundgas/total_metals], psym=thispsym, color=thiscolor, thick=3.0
;oplot, [totalmass], [metals_in_xraygas/total_metals], psym=thispsym, color=thiscolor, thick=3.0
oplot, [totalmass], [metals_in_stars/total_metals], psym=thispsym, color=thiscolor, thick=3.0


end




;--------------------------------------------------------------------------




;==================================
;  Load the ID list 
;   (works for either ID list)
;==================================
function load_id_list, idfilename


spawn, 'wc '+idfilename, result
lines= long(result)
idlist= lonarr(lines-5)

get_lun, unit
openr, unit, idfilename

textjunk= ''
readf, unit, textjunk
readf, unit, textjunk
readf, unit, textjunk
readf, unit, textjunk
readf, unit, textjunk
readf, unit, idlist
close, unit
free_lun, unit


return, idlist


end




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

bins = 50

xmax= 600.0
xmin= 30.0
xaxistitle='!7r!6 (km s!E-1!N)'

ymax = 1.2
ymin = 0.2
yaxistitle= "!6<Z> (Z!D!9n!6!N)"


if not keyword_set(snapnum) then snapnum= 25

; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
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

plot_series_mets, fruns, snapnums, pointt=2




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

plot_series_mets, fruns, snapnums, pointt=1



; e, f=0.4, q=0.25
; ------------------
;fruns= ["/raid4/tcox/ds/d0e2_q", $
;	"/raid4/tcox/ds/d1e2_q", $
;	"/raid4/tcox/ds/d2e2_q", $
;	"/raid4/tcox/ds/d3e7", $
;	"/raid4/tcox/ds/d4e2_q", $
;	"/raid4/tcox/ds/d5e2_q", $
;	"/raid4/tcox/ds/d6e2_q"]

;snapnums= [50, 30, 30, 30, 30, 30, 30]

;plot_series_mets, fruns, snapnums, pointt=2



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

plot_series_mets, fruns, snapnums, pointt=3


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

plot_series_mets, fruns, snapnums, pointt=5






; print extras
; -------------

;xyouts, 0.55, 0.90, 'no black hole', size=1.5, color=0, /normal, charthick=4.0
;xyouts, 0.55, 0.85, '!6time course', size=1.5, color=0, /normal, charthick=4.0
;xyouts, 0.55, 0.80, 'dark -> ', size=1.5, color=20, /normal, charthick=4.0
;xyouts, 0.75, 0.80, ' light', size=1.5, color=160, /normal, charthick=4.0

;xyouts, 0.35, 0.27, 'B-luminosity weighted', size=1.5, color=0, /normal, charthick=4.0

xyouts, 0.45, 0.35, 'New Stars', size=1.5, color=0, /normal, charthick=4.0
xyouts, 0.45, 0.31, 'mass weighted (same as)', size=1.5, color=0, /normal, charthick=4.0
xyouts, 0.45, 0.27, 'unweighted average', size=1.5, color=0, /normal, charthick=4.0

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
do_lum_weighted= 1
;do_lum_weighted= 0
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



