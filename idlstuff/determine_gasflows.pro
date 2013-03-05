


pro doit, junk

	frun= "data/ds/vc3vc3f"
	mass_v_radius, frun
	plot_mass_within_r, frun

end




;===================================================================================

pro mass_v_radius, frun


if not keyword_set(frun) then begin
        print, " "
        print, " mass_v_radius, frun"
        print, " "
        print, " "
        print, " "
        return
endif


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snapshot* | wc ",result
nsnaps=long(result[0])


rarray= [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, $
	 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, $
	 20.0, 25.0, 30.0, 50.0]
Nrarray= n_elements(rarray)


;
; use BH information to determine
; the center
;
ok= fload_snapshot_bh(frun,0,/nopot_in_snap,/skip_center)
init_bhids= fload_blackhole_id(1)
init_bhids= long(init_bhids(sort(init_bhids)))
n_bh= n_elements(init_bhids)
print, "Blackhole ID's: ", init_bhids 



; ----------------------------------------
; This part loops through the snapshots
; and compiles mass v. radius information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"

        ; open snapshot
        ;ok=fload_snapshot(frun,snapnum)
        ok=fload_snapshot_bh(frun,i,/nopot_in_snap,/skip_center)

	exts= '0000'+strcompress(string(i), /remove_all)
	istring= strmid(exts,strlen(exts)-3,3)

        ; what time is it?
        time= fload_time(1)


	; now, write a file - one per snapshot per blackhole, that
	; contains the mass within specific radii for each component

        ;--------------------------------------------
	;  disk 1

	sid= long(1)
	npart= init_bhids[0]
	bhlbl= strcompress(string(init_bhids[0]),/remove_all)

	process_snap_mass_profile, frun, istring, bhlbl, sid, npart, Nrarray, rarray, time



	;--------------------------------------------
	; disk 2
	;

	sid= init_bhids[0]+long(1)
	npart= init_bhids[1]-init_bhids[0]
	bhlbl= strcompress(string(init_bhids[1]),/remove_all)

	process_snap_mass_profile, frun, istring, bhlbl, sid, npart, Nrarray, rarray, time


endfor


; ----------------------------------------
        
print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"
        




end







;===================================================================================










;---------------------------------------
; Process the amount of mass within r
;---------------------------------------
pro process_snap_mass_profile, frun, istring, bhlbl, sid, npart, $
			Nrarray, rarray, time

        center1= fload_1gal_center(sid,npart)

        rgas= fload_1gal_gas_xyz('r',sid,npart,center=center1)
        mgas= fload_1gal_gas_mass(1,sid,npart)
        rdisk= fload_1gal_disk_xyz('r',sid,npart,center=center1)
        mdisk= fload_1gal_disk_mass(1,sid,npart)
        rbulge= fload_1gal_bulge_xyz('r',sid,npart,center=center1)
        mbulge= fload_1gal_bulge_mass(1,sid,npart)
        rnewstars= fload_1gal_newstars_xyz('r',sid,npart,center=center1)
        mnewstars= fload_1gal_newstars_mass(1,sid,npart)
        rdm= fload_1gal_halo_xyz('r',sid,npart,center=center1)
        mdm= fload_1gal_halo_mass(sid,npart)


	centerlbl= string(center1[0])+', '+string(center1[1])+', '+string(center1[2])

        openw, 1, frun+'/massprofile_'+bhlbl+'_'+istring+'.txt', ERROR=err

        printf, 1, "#   massprofile.txt     snap="+istring+ "   bh= "+bhlbl
        printf, 1, "#      ---> mass, in gadget units (10^10 msolar/h), within each radius "
        printf, 1, "#     center=  "+centerlbl
        printf, 1, "#  radius   "
        printf, 1, "# (kpc/h)      gas          disk        bulge     newstars           dm  "

        ; now, do the radii
        for i=0,Nrarray-1 do begin

                idx= where(rgas le rarray[i])
                if idx(0) ne -1 then mass_gas= total(mgas(idx)) else mass_gas= 0.0
                idx= where(rdisk le rarray[i])
                if idx(0) ne -1 then mass_disk= total(mdisk(idx)) else mass_disk= 0.0
                idx= where(rbulge le rarray[i])
                if idx(0) ne -1 then mass_bulge= total(mbulge(idx)) else mass_bulge= 0.0
                idx= where(rnewstars le rarray[i])
                if idx(0) ne -1 then mass_newstars= total(mnewstars(idx)) else mass_newstars= 0.0
                idx= where(rdm le rarray[i])
                if idx(0) ne -1 then mass_darkmatter= total(mdm(idx)) else mass_darkmatter= 0.0

                printf, 1, FORMAT= '("  ",F6.3, 5(F11.6,"  "))', $
                        rarray[i], mass_gas, mass_disk, mass_bulge, mass_newstars, mass_darkmatter
        endfor
        close, 1


end







;===================================================================================








; ----------------------------
;  Read hotgas.txt file
; ----------------------------
pro read_massprofile_file, filename, radius, mass_g, mass_d, mass_b, mass_ns, mass_dm

thisfile= filename

spawn, "wc "+thisfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then hgas_data= fltarr(6,lines)

openr, 1, thisfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, hgas_data
close, 1


radius= hgas_data[0,*]
mass_g= hgas_data[1,*]
mass_d= hgas_data[2,*]
mass_b= hgas_data[3,*]
mass_ns= hgas_data[4,*]
mass_dm= hgas_data[5,*]


end






;===================================================================================




pro plot_mass_within_r, frun


if not keyword_set(frun) then begin
	print, " "
	print, " plot_mass_within_r, frun"
	print, " "
	print, " "
	return
endif


filename='massinr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


yaxistitle = "Mass (h!E-1!N10!E10!N)"
xaxistitle = "Time (h!E-1!NGyr)"


radius= 2.0

xmax = 2.2
xmin = 0

;ymax = 5.0
ymax = 3.2
;ymax = 2.0
;ymin = 0.0
ymin= 0.0
;ymin= 0.001


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


;
; ---------------


;plot_one_massprofile, frun, "gas", radius, /dototal, zcolor= 150
;plot_one_massprofile, frun, "gas+newstars", radius, /dototal, zcolor= 150


; plot: massprofile.eps (uncomment the following)
;----------------------------
xyouts, 0.65, 0.52, "Mass w/i r=", /normal, size=1.4, color=0, charthick= 3.0
plot_one_massprofile, frun, "gas+newstars", 0.5, /dototal, zcolor= 200
xyouts, 0.65, 0.48, '0.5 kpc/h', /normal, size=1.2, color= 200, charthick= 3.0
plot_one_massprofile, frun, "gas+newstars", 1.0, /dototal, zcolor= 175
xyouts, 0.65, 0.44, '1.0 kpc/h', /normal, size=1.2, color= 175, charthick= 3.0
plot_one_massprofile, frun, "gas+newstars", 1.5, /dototal, zcolor= 150
xyouts, 0.65, 0.40, '1.5 kpc/h', /normal, size=1.2, color= 150, charthick= 3.0
plot_one_massprofile, frun, "gas+newstars", 2.0, /dototal, zcolor= 125
xyouts, 0.65, 0.36, '2.0 kpc/h', /normal, size=1.2, color= 125, charthick= 3.0
plot_one_massprofile, frun, "gas+newstars", 2.5, /dototal, zcolor= 100
xyouts, 0.65, 0.32, '2.5 kpc/h', /normal, size=1.2, color= 100, charthick= 3.0
plot_one_massprofile, frun, "gas+newstars", 3.0, /dototal, zcolor=  75
xyouts, 0.65, 0.28, '3.0 kpc/h', /normal, size=1.2, color=  75, charthick= 3.0
plot_one_massprofile, frun, "gas+newstars", 5.0, /dototal, zcolor= 50
xyouts, 0.65, 0.24, '5.0 kpc/h', /normal, size=1.2, color=  50, charthick= 3.0
plot_one_massprofile, frun, "gas+newstars", 10.0, /dototal, zcolor= 0
xyouts, 0.65, 0.20, '10.0 kpc/h', /normal, size=1.2, color= 0, charthick= 3.0



; plot: massprofile_pergal.eps (uncomment the following)
;----------------------------
;plot_one_massprofile, frun, "gas+newstars", 1.5, zcolor= 150
;xyouts, 0.65, 0.30, '1.5 kpc/h', /normal, size=1.2, color= 150, charthick= 3.0
;xyouts, 0.65, 0.26, 'gal 1 - dotted', /normal, size=1.2, color= 150, charthick= 3.0
;xyouts, 0.65, 0.22, 'gal 2 - dashed', /normal, size=1.2, color= 150, charthick= 3.0



;  plot extras
; -------------
xyouts, 0.25, 0.88, fload_getid(frun), /normal, charthick=4.0, size=1.4, color=0

;rlbl= 'R<='+strcompress(string(radius),/remove_all)
;xyouts, 0.6, 0.80, rlbl, /normal, size=1.2, color= 150, charthick= 3.0



device, /close



end










;===================================================================================


pro plot_one_massprofile, frun, component, radius, dototal=dototal, zcolor=zcolor

; ----------------------------

spawn, "/bin/ls "+frun+"/massprofile_*_000.txt",result
if strlen(result[0]) le 0 then begin
	print, " "
	print, " PROBLEM: can't find any 'massprofile' files, try running the 'mass_v_radius' script"
	print, " "
	return
endif

; ----------------------------

; determine the number of
; snapshots in frun directory
spawn, "/bin/ls "+frun+"/snapshot* | wc ",result
nsnaps=long(result[0])

; ----------------------------

; grab BH info
;
ok= fload_snapshot_bh(frun,0,/nopot_in_snap,/skip_center)
init_bhids= fload_blackhole_id(1)
init_bhids= long(init_bhids(sort(init_bhids)))
n_bh= n_elements(init_bhids)
print, "Blackhole ID's: ", init_bhids

; ----------------------------


time= fltarr(nsnaps)
mass1= fltarr(nsnaps)
mass2= fltarr(nsnaps)
mass= fltarr(nsnaps)

for i=0, nsnaps-1 do begin

        print, "--------------------------------------"

        ; open snapshot
        ok=fload_snapshot_bh(frun,i,/nopot_in_snap,/skip_center,/header_only)

        exts= '0000'+strcompress(string(i), /remove_all)
        istring= strmid(exts,strlen(exts)-3,3)

        ; what time is it?
        time[i]= fload_time(1)


        ;--------------------------------------------
        ;  disk 1
        bhlbl= strcompress(string(init_bhids[0]),/remove_all)
	filename= frun+"/massprofile_"+bhlbl+"_"+istring+".txt"
	mass1[i]= grab_one_massprofile(filename,component,radius)

        ;--------------------------------------------
        ;  disk 2
        bhlbl= strcompress(string(init_bhids[1]),/remove_all)
	filename= frun+"/massprofile_"+bhlbl+"_"+istring+".txt"
	mass2[i]= grab_one_massprofile(filename,component,radius)

endfor


if keyword_set(dototal) then begin
	mass= mass1 + mass2
	oplot, time, mass, psym=-3, linestyle= 0, color=zcolor, thick= 4.0
endif else begin
	oplot, time, mass1, psym=-3, linestyle= 2, color=zcolor
	oplot, time, mass2, psym=-3, linestyle= 1, color=zcolor
endelse


end




;===================================================================================



function grab_one_massprofile, filename, component, xradius

read_massprofile_file, filename, radius, mass_g, mass_d, mass_b, mass_ns, mass_dm

idx= where(radius gt (xradius-0.001))
if idx(0) ne -1 then idxx= idx(0) else return, -1

if component eq 'gas+newstars' then return, mass_g(idxx)+mass_ns(idxx)
if component eq 'gas' then return, mass_g(idxx)
if component eq 'disk' then return, mass_d(idxx)
if component eq 'bulge' then return, mass_b(idxx)
if component eq 'newstars' then return, mass_ns(idxx)
if component eq 'dm' then return, mass_dm(idxx)

end



;===================================================================================


