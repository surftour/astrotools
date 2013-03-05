pro id_image, junk


   frun="/data/tcox/Sbc201a-u4"

   idsnap=60
   ;idsnap_idfilename= 'id_list.txt'
   ;generate_idlist_newstar_r, frun, idsnap, idlistname=idsnap_idfilename
   generate_idlist_newstar_r, frun, idsnap

   imagesnap=56
   ;imagesnap=0
   ;xlen= 120.0
   xlen= 80.0
   filename='contour.eps'
   listofids=load_id_list(idsnap_idfilename)

   ; or, should I just do this in procedure
   contour_allstars, frun, imagesnap, xlen, 'ps', filename=filename, /nolabels

end






;=========================================
;
;  Do the grunt work of figuring out
; what the id's are of certain particles
;
;=========================================

;
; Write a file with list of
;  new star particle ids that
;  satisfy some criteria
;-------------------------------
pro generate_idlist_newstar_r, frun, snapnum, idlistname=idlistname 

if not keyword_set(frun) then begin
   print, "  "
   print, "generate_idlist_newstar_r, frun, snapnum, idlistname=idlistname"
   return
endif

ok= fload_snapshot(frun,snapnum)
;ok= fload_snapshot_bh(frun,snapnum)


;do_newstars= 1
do_newstars= 0
if do_newstars eq 1 then begin
    nsids= fload_newstars_id(1)
    pids= fload_newstars_parentid(1)
    ;r= fload_newstars_xyz('r')
    ;r= fload_newstar_xyz('rxy')
    e_tot= fload_newstars_energy(1)
endif

;do_oldstars= 1
do_oldstars= 0
if do_oldstars eq 1 then begin
    nsids= fload_oldstars_id(1)
    pids= nsids*0.0
    e_tot= fload_oldstars_energy(1)
endif


do_allstars= 1
;do_allstars= 0
if do_allstars eq 1 then begin
    nsids= [fload_newstars_id(1),fload_oldstars_id(1)]
    nadaids= 0.0*fload_oldstars_id(1)
    pids= [fload_newstars_parentid(1),nadaids]
    r= [fload_newstars_xyz('r'),fload_oldstars_xyz('r')]
    ;e_tot= [fload_newstars_energy(1),fload_oldstars_energy(1)]
endif



; now order these
; -----------------
;eindx= sort(e_tot)
;e_tot= e_tot(eindx)
eindx= sort(r)
r= r(eindx)
nsids= nsids(eindx)
pids= pids(eindx)


; id selection
; --------------

; trial 1, at a specific radius
; ------------------------------
do_radius_decomp = 1
;do_radius_decomp = 0
if do_radius_decomp eq 1 then begin

	; ----------------------------
	;  Make sure we're getting 
	;  original ID's
	for idx=0L,n_elements(pids)-1 do begin
	   if ((pids[idx] lt nsids[idx]) and (pids[idx] gt 0)) then nsids[idx]= pids[idx]
	endfor

	R_e= 5.0

	; ----------------------------------------
	idx= where((r ge 0.0) and (r le R_e))
	idlist= nsids(idx)

	id_filename= "id_list_0to1.txt"
	print, ' writing 0-1 R_e ids to '+id_filename
	write_id_list, idlist, id_filename, cmt=cmt

	; ----------------------------------------
	idx= where((r gt R_e) and (r le 2.0*R_e))
	idlist= nsids(idx)

	id_filename= "id_list_1to2.txt"
	print, ' writing 1-2 R_e ids to '+id_filename
	write_id_list, idlist, id_filename, cmt=cmt

	; ----------------------------------------
	idx= where((r gt 2.0*R_e) and (r le 3.0*R_e))
	idlist= nsids(idx)

	id_filename= "id_list_2to3.txt"
	print, ' writing 2-3 R_e ids to '+id_filename
	write_id_list, idlist, id_filename, cmt=cmt

	; ----------------------------------------
	idx= where((r gt 3.0*R_e) and (r le 4.0*R_e))
	idlist= nsids(idx)

	id_filename= "id_list_3to4.txt"
	print, ' writing 3-4 R_e ids to '+id_filename
	write_id_list, idlist, id_filename, cmt=cmt

	; ----------------------------------------
	idx= where((r gt 4.0*R_e) and (r le 5.0*R_e))
	idlist= nsids(idx)

	id_filename= "id_list_4to5.txt"
	print, ' writing 4-5 R_e ids to '+id_filename
	write_id_list, idlist, id_filename, cmt=cmt

	; ----------------------------------------
	idx= where((r gt 5.0*R_e) and (r le 6.0*R_e))
	idlist= nsids(idx)

	id_filename= "id_list_5to6.txt"
	print, ' writing 5-6 R_e ids to '+id_filename
	write_id_list, idlist, id_filename, cmt=cmt

	; ----------------------------------------
	idx= where(r gt 6.0*R_e)
	idlist= nsids(idx)

	id_filename= "id_list_gt6.txt"
	print, ' writing 6+ R_e ids to '+id_filename
	write_id_list, idlist, id_filename, cmt=cmt

endif



; trial 2, in energy bins
; --------------------------
; we'll use equal particle number
; bins for right now
do_energy_decomp= 0
;do_energy_decomp= 1
if do_energy_decomp eq 1 then begin
  bins= 20

  pts_in_bin= n_elements(nsids)/bins

  for i=0,bins-1 do begin

	si= i*pts_in_bin
	ei= si+pts_in_bin-1
	if i eq bins-1 then ei= n_elements(nsids)-1

	be= e_tot[si:ei]
	be_moment= moment(be)
	cmt= 'average energy= '+strcompress(string(be_moment[0]),/remove_all)

	idlist= nsids[si:ei]
	pidlist= pids[si:ei]

	; ----------------------------
	;  Make sure we're getting 
	;  original ID's
	for idx=0,n_elements(idlist)-1 do begin
	   if ((pidlist[idx] lt idlist[idx]) and (pidlist[idx] gt 0)) then idlist[idx]= pidlist[idx]
	endfor


	; ----------------------------
	;  Write id list to a file
	;
	ifnum= '000'+strcompress(string(i),/remove_all)
	ifnum= strmid(ifnum,strlen(ifnum)-3,3)
	id_filename= "id_list_"+ifnum+".txt"
	print, "i= ",i,' writing '+cmt+' to '+id_filename
	write_id_list, idlist, id_filename, cmt=cmt

  endfor
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
pro plot_newstar_age_e, frun, snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "plot_newstar_age_e, frun, snapnum"
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



