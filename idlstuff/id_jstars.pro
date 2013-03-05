;=========================================
;
;  Do the grunt work of figuring out
; what are the id's of certain particles,
; and then write these to a files.
;
;=========================================

;
; Write a file with list of
;  new star particle ids that
;  satisfy some criteria
;-------------------------------
pro generate_idlist_j, frun, snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "generate_idlist_j, frun, snapnum"
   return
endif

;ok= fload_snapshot(frun,snapnum)
ok= fload_snapshot_bh(frun,snapnum)


if fload_npart(4) gt 0 then begin

    ; load new star info
    nsids= fload_newstars_id(1)
    r= fload_newstars_xyz('r')
    jz= fload_newstars_j(23)

    ; write id list
    ; ----------------------------------------
    idx= where((r le 3.0) and (jz le 0.0))
    idlist= nsids(idx)
    id_filename= frun+"/id_list_ns_j_p_in.txt"
    print, ' writing r<3.0, jz>0 new star ids to '+id_filename
    write_id_list, idlist, id_filename, cmt=cmt

    ; ----------------------------------------
    idx= where((r gt 3.0) and (jz le 0.0))
    idlist= nsids(idx)
    id_filename= frun+"/id_list_ns_j_p_out.txt"
    print, ' writing r>3.0, jz>0 new star ids to '+id_filename
    write_id_list, idlist, id_filename, cmt=cmt

    ; ----------------------------------------
    idx= where((r le 3.0) and (jz le 0.0))
    idlist= nsids(idx)
    id_filename= frun+"/id_list_ns_j_m_in.txt"
    print, ' writing r<3.0, jz<0 new star ids to '+id_filename
    write_id_list, idlist, id_filename, cmt=cmt

    ; ----------------------------------------
    idx= where((r le 3.0) and (jz le 0.0))
    idlist= nsids(idx)
    id_filename= frun+"/id_list_ns_j_m_out.txt"
    print, ' writing r>3.0, jz<0 new star ids to '+id_filename
    write_id_list, idlist, id_filename, cmt=cmt

    ; ----------------------------------------
endif




if fload_npart(2)+fload_npart(3) gt 0 then begin

    ; load old star info
    osids= fload_oldstars_id(1)
    r= fload_oldstars_xyz('r')
    jz= fload_oldstars_j(23)

    ; write id list
    ; ----------------------------------------
    idx= where((r le 3.0) and (jz le 0.0))
    idlist= osids(idx)
    id_filename= frun+"/id_list_os_j_p_in.txt"
    print, ' writing r<3.0, jz>0 old star ids to '+id_filename
    write_id_list, idlist, id_filename, cmt=cmt

    ; ----------------------------------------
    idx= where((r gt 3.0) and (jz le 0.0))
    idlist= osids(idx)
    id_filename= frun+"/id_list_os_j_p_out.txt"
    print, ' writing r>3.0, jz>0 old star ids to '+id_filename
    write_id_list, idlist, id_filename, cmt=cmt

    ; ----------------------------------------
    idx= where((r le 3.0) and (jz le 0.0))
    idlist= osids(idx)
    id_filename= frun+"/id_list_os_j_m_in.txt"
    print, ' writing r<3.0, jz<0 old star ids to '+id_filename
    write_id_list, idlist, id_filename, cmt=cmt

    ; ----------------------------------------
    idx= where((r le 3.0) and (jz le 0.0))
    idlist= osids(idx)
    id_filename= frun+"/id_list_os_j_m_out.txt"
    print, ' writing r>3.0, jz<0 old star ids to '+id_filename
    write_id_list, idlist, id_filename, cmt=cmt

    ; ----------------------------------------
endif




end






;==================================
;  Write the ID list
;==================================
pro write_id_list, idlist, idfilename, cmt=cmt

if not keyword_set(cmt) then cmt=' '

; ----------------------------
;  Write id list to a file
;
get_lun, unit
openw, unit, idfilename

printf, unit, '#  file: '+idfilename
printf, unit, "#  list of new star ids "
printf, unit, "#  "
printf, unit, '#  '+cmt
printf, unit, "#  "
for i=0L,n_elements(idlist)-1 do printf, unit, idlist[i]

close, unit

end





;==================================
;  Load the ID list
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


return, idlist


end





;==============================================
;
;  Plots to display newstar properties as 
;  a function of binding energy.
;
;==============================================


;
; Generate plot showing the new star particle
; ages, as a function of binding energy
;----------------------------------------------
pro plot_newstar_j_age, frun, snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "plot_newstar_j_age, frun, snapnum"
   return
endif

ok= fload_snapshot(frun,snapnum)

;nsids= fload_newstars_id(1)
;pids= fload_newstars_parentid(1)
;r= fload_newstars_xyz('r')
;r= fload_newstar_xyz('rxy')

;pe= fload_newstar_energy(1,/potential)
;ke= fload_newstar_energy(1,/kinetic)
e_tot= fload_newstars_energy(1)

; this is actually, the formation time
age= fload_newstars_age(1)
; convert to age
Ttime= fload_time(1)
age= Ttime-age


; now order these
; -----------------
eindx= sort(e_tot)
e_tot= e_tot(eindx)
age= age(eindx)


; now bin data
; --------------

; we'll use equal particle number
; bins for right now

bins= 20

pts_in_bin= n_elements(age)/bins

Aage= fltarr(bins)
Aageerr= fltarr(bins)
binde= fltarr(bins)

for i=0,bins-1 do begin

	si= i*pts_in_bin
	ei= si+pts_in_bin-1
	if i eq bins-1 then ei= n_elements(age)-1

	be= e_tot[si:ei]
	be_moment= moment(be)
	binde[i]= be_moment[0]

	ag= age[si:ei]
	ag_moment= moment(ag)
	Aage[i]= ag_moment[0]
	Aageerr[i]= sqrt(ag_moment[1])
endfor



;----------------
; Print it up
;----------------

;filename=frun+'/newstarage.eps'
filename='newstarage.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

yaxistitle = "Average Age (Gyr)"
xaxistitle = "Binding Energy (GU)"
xmax = 10.0
xmin = -50.0
ymax = 4.0
ymin = 0.0

;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; center 1
thispsym= 6
thiscolor= 70
oplot, binde, Aage, thick=4.0, psym=-thispsym, color=thiscolor
oplot, binde, Aage+Aageerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
oplot, binde, Aage-Aageerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1

xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0

device, /close


end





;
; Generate plot showing the new star particle
; radius, as a function of binding energy
;----------------------------------------------
pro plot_newstar_radius_e, frun, snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "plot_newstar_radius_e, frun, snapnum"
   return
endif

ok= fload_snapshot(frun,snapnum)

;do_newstars= 1
do_newstars= 0
if do_newstars eq 1 then begin
    ;nsids= fload_newstars_id(1)
    ;pids= fload_newstars_parentid(1)
    ;r= fload_newstars_xyz('r')
    ;r= fload_newstar_xyz('rxy')

    ;pe= fload_newstars_energy(1,/potential)
    ;ke= fload_newstars_energy(1,/kinetic)
    e_tot= fload_newstars_energy(1)

    ; this is actually, the formation time
    radius= fload_newstars_xyz('r')
endif

do_oldstars= 1
;do_oldstars= 0
if do_oldstars eq 1 then begin
    e_tot= fload_oldstars_energy(1)
    radius= fload_oldstars_xyz('r')
endif


; now order these
; -----------------
eindx= sort(e_tot)
e_tot= e_tot(eindx)
radius= radius(eindx)


; now bin data
; --------------

; we'll use equal particle number
; bins for right now

bins= 20

pts_in_bin= n_elements(radius)/bins

Aradius= fltarr(bins)
Aradiuserr= fltarr(bins)
binde= fltarr(bins)

for i=0,bins-1 do begin

	si= i*pts_in_bin
	ei= si+pts_in_bin-1
	if i eq bins-1 then ei= n_elements(radius)-1

	be= e_tot[si:ei]
	be_moment= moment(be)
	binde[i]= be_moment[0]

	ag= radius[si:ei]
	ag_moment= moment(ag)
	Aradius[i]= ag_moment[0]
	Aradiuserr[i]= sqrt(ag_moment[1])
endfor



;----------------
; Print it up
;----------------

;filename=frun+'/newstarr.eps'
filename='newstarr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

yaxistitle = "Average Radius (kpc)"
xaxistitle = "Binding Energy (GU)"
xmax = 10.0
xmin = -50.0
ymax = 100.0
ymin = 0.01

;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; center 1
thispsym= 6
thiscolor= 70
oplot, binde, Aradius, thick=4.0, psym=-thispsym, color=thiscolor
oplot, binde, Aradius+Aradiuserr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
oplot, binde, Aradius-Aradiuserr, thick=1.0, psym=-3, color=thiscolor, linestyle=1

xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0

device, /close


end


















;
;  Find Intersection of two arrays
; ------------------------------------
FUNCTION SetIntersection, a, b
minab = Min(a, Max=maxa) > Min(b, Max=maxb) ;Only need intersection of ranges
maxab = maxa < maxb

   ; If either set is empty, or their ranges don't intersect: result = NULL.

IF maxab LT minab OR maxab LT 0 THEN RETURN, -1
r = Where((Histogram(a, Min=minab, Max=maxab) NE 0) AND  $
          (Histogram(b, Min=minab, Max=maxab) NE 0), count)

IF count EQ 0 THEN RETURN, -1 ELSE RETURN, r + minab
END



