pro map, initialdiskrun, frun, snapnum, filename=filename, xmax=xmax, ymax=ymax, $
					newstars=newstars, oldstars=oldstars, $
					diffaxes=diffaxes

if not keyword_set(filename) then filename="gasmapping.eps"
if not keyword_set(snapnum) then snapnum=0
if not keyword_set(frun) or not keyword_set(initialdiskrun) then begin
   print, "  "
   print, "map, initialdiskrun, frun, snapnum, filename=filename"
   print, "  "
   return
endif





;-------------------------------------
;  Load 0 Snapshot
;-------------------------------------

if fload_snapshot_bh(initialdiskrun, 0, /nopot_in_snap) lt 0 then begin
	print, "PROBLEM: opening file"
	return
endif


ids_init= fload_gas_id(1)
r_init= fload_gas_xyz('r')

if keyword_set(oldstars) then begin
	ids_init= fload_oldstars_id(1)
	r_init= fload_oldstars_xyz('r')
endif

sortidx= sort(ids_init)
ids_init= ids_init(sortidx)
r_init= r_init(sortidx)
print, "n(ids_init)= ", n_elements(ids_init)

print, "fload_npart(0)= ", fload_npart(0)
print, "fload_npart(1)= ", fload_npart(1)
print, "fload_npart(2)= ", fload_npart(2)
print, "fload_npart(3)= ", fload_npart(3)
print, "fload_npart(4)= ", fload_npart(4)
print, "fload_npart(5)= ", fload_npart(5)
progenitor_n_gal= total(fload_npart(99))
print, "isolated galaxy progenitor_n_gal= ", progenitor_n_gal

;progenitor_highest_id= ids_init[n_elements(ids_init)-1]
;print, "isolated galaxy progenitor_highest_id= ", progenitor_highest_id


;-------------------------------------
;  Load Snapshot x
;-------------------------------------

if fload_snapshot_bh(frun, snapnum, /nopot_in_snap) lt 0 then begin
	print, "PROBLEM: opening file"
	return
endif


yaxistitle='!6Final Radius of Gas (kpc)'
ids_final= fload_gas_id(1)
r_final= fload_gas_xyz('r')

if keyword_set(newstars) then begin
	yaxistitle= '!6Final Radius of New Stars (kpc)'
	ids_final= fload_newstars_id(1)
	r_final= fload_newstars_xyz('r')
endif

if keyword_set(oldstars) then begin
	yaxistitle= '!6Final Radius of Old Stars (kpc)'
	ids_final= fload_oldstars_id(1)
	r_final= fload_oldstars_xyz('r')
endif

sortidx= sort(ids_final)
ids_final= ids_final(sortidx)
r_final= r_final(sortidx)

print, "fload_npart(0)= ", fload_npart(0)
print, "fload_npart(1)= ", fload_npart(1)
print, "fload_npart(2)= ", fload_npart(2)
print, "fload_npart(3)= ", fload_npart(3)
print, "fload_npart(4)= ", fload_npart(4)
print, "fload_npart(5)= ", fload_npart(5)
print, "total(fload_npart(99))= ", total(fload_npart(99))
print, "frun n_particle_to_search_for_in_final_remnant= ", fload_npart(0)


;
;
; Now, do the mapping between initial and final
;
print, " "
print, " loop N= ", n_elements(ids_final)
print, " "
initialradius= fltarr(n_elements(ids_final))
finalradius= fltarr(n_elements(ids_final))
initi= 0L
for i=0L, n_elements(ids_final)-1 do begin

	if (i mod 50000) eq 0 then print, "i= ",i

	thisid= ids_final[i]
	if thisid gt progenitor_n_gal then thisid= thisid - progenitor_n_gal

	finalradius[i]= r_final[i]

	loopi= 1
	;idx= where(ids_init eq thisid)
	while thisid ne ids_init[initi] do begin

		initi = initi + 1

		if initi ge n_elements(ids_init) then begin
			initi= 0L
			loopi = loopi + 1
		endif

		if loopi gt 3 then begin
			print, " "
			print, " STOPPING: can't find indx after running through list twice"
			print, " "
			stop
		endif
	endwhile

	initialradius[i]= r_init[initi]

	if i lt 10 then begin
	print, "---"
	print, "init  ", initi, ids_init[initi], r_init[initi], initialradius[i]
	print, "final ", i, thisid, ids_final[i], r_final[i], finalradius[i]
	endif

endfor


; manually do log
;initialradius= alog10(initialradius)
;finalradius= alog10(finalradius)
;xmin= -1.5  & xmax= 2.0
;ymin= -1.5  & ymax= 2.5


;----------------------------------------
; Send it to plotting routines
;----------------------------------------

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename=filename, colortable=4

if keyword_set(oldstars) then fload_newcolortable, 10
if keyword_set(newstars) then fload_newcolortable, 25


xaxistitle='!6Initial Radius of Gas (kpc)'
if not keyword_set(xmax) then xmax = 100.0
if not keyword_set(xmin) then xmin = 0.0
if not keyword_set(ymax) then ymax = 100.0
if not keyword_set(ymin) then ymin = 0.0

if keyword_set(diffaxes) then begin
	; log
	;xmin= 0.2
	;xmax= 3.0
	;ymin= 0.2
	;ymax= 3.0
	;xaxistitle='!6 Log Initial Radius (kpc)'
	;yaxistitle='!6 Log Final Radius (kpc)'
	;initialradius= alog10(initialradius)
	;finalradius= alog10(finalradius)

	; 1/r
	xmin= -3.0
	xmax= 1.3
	ymin= -3.0
	ymax= 1.3
	xaxistitle='!6 1 / Initial Radius (kpc)'
	yaxistitle='!6 1 /  Final Radius (kpc)'
	initialradius=  1.0 / initialradius
	finalradius= 1.0 / finalradius
	initialradius= alog10(initialradius)
	finalradius= alog10(finalradius)

endif

contour_makegeneralplot, initialradius, finalradius, xmax, xmin, ymax, ymin, 'ps', filename=filename, $
                        ;msg=msg, $
                        ;xaxistitle='!6Initial Radius (kpc)', $
                        ;xaxistitle='!6Initial Radius of Gas (kpc)', $
			xaxistitle=xaxistitle, $
			;yaxistitle='!6Final Radius of Gas (kpc)', $
			yaxistitle=yaxistitle, $
			;yaxistitle='!6Final Radius of New Stars (kpc)', $
                        ;xaxistitle='!6Log Initial Radius (kpc)', $
			;yaxistitle='!6Log Final Radius (kpc)', $
                        colortable=colortable

xyouts, 0.20, 0.88, frun, size=1.33, color=0, /normal, charthick= 3.0


device, /close



end


