
pro plot_one_disk, x, y, clr, $
			big= big, $
			small=small

	;
        ; plot ring
	;
        symsize= 0.6
        if keyword_set(small) then symsize= 0.2
        if keyword_set(big) then symsize= 1.2
        ;
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
        ;
        oplot, x, y, psym=-8, color=clr
        ;oplot, x, y, psym=-8, color= 150
        ;if i le 3 then oplot, x, y, psym=-8, color= 150
        ;if i gt 3 then oplot, x, y, psym=1, color= 150
        ;if i eq 4 then oplot, x, y, psym=-1, color= 0, linestyle= 1
        ;if i eq 5 then oplot, x, y, psym=-2, color= 0, linestyle= 2


        ;
        ; plot theta= 0 point differently
        ;
        symsize= 1.2
        if keyword_set(small) then symsize= 0.2
        if keyword_set(big) then symsize= 1.8
	;
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        ;
	;fload_newcolortable, 4
        oplot, [x[0]], [y[0]], psym=8, color= clr
        ;oplot, [x[0]], [y[0]], psym=8, color= 50
        ;if i le 3 then oplot, [x[0]], [y[0]], psym=8, color= 50
        ;if i gt 3 then oplot, [x[0]], [y[0]], psym=1, color= 50
        ;if i eq 4 then oplot, [x[0]], [y[0]], psym=1, color= 50, symsize= 2.0, thick= 4.0
        ;if i eq 5 then oplot, [x[0]], [y[0]], psym=2, color= 50, symsize= 2.0, thick= 4.0
    
end









; does the work for each panel
; -------------------------------
pro toomre_do_one_panel, frun, snapnum, xlen, $
	x0, y0, xs, ys, $
	draw_rings=draw_rings, $
	center=center, $
	pt1center=pt1center, $
	pt2center=pt2center, $
	remove_disk1=remove_disk1, $
	remove_disk2=remove_disk2, $
	remove_pt1=remove_pt1, $
	remove_pt2=remove_pt2, $
	remove_traj1=remove_traj1, $
	remove_traj2=remove_traj2



    if not keyword_set(center) then center=[0,0,0]


    ;  get data
    ; -------------
    if not keyword_set(loadedsnap) then begin
      ok= -1
      ;if ok lt 0 then ok= fload_snapshot_bh("/n/home/tcox/Tools/MakeToomreDisk/toomredisk.dat", 0, /ics, /nopot_in_snap)
      if ok lt 0 then ok= fload_snapshot_bh(frun, snapnum, /nopot_in_snap, /skip_center)
      if ok eq 1 then begin
        print, "PROBLEM: opening file"
        return
      endif
    endif



    ; now, plot the sucker
    ; -------------------------
    x1= x0 + xs
    y1= y0 + ys
    oplot_one_axes, x0, y0, x1, y1, "-1", xlen    ;, /add_j_vector

    ; time 
    xyouts, x0+0.025, y1-0.06, fload_timelbl(0.7,3,/noteq), /normal, size= 1.5, charthick=3.0, color= 0




    ; central mass
    ; ----------------
    cmx= fload_halo_xyz('x',center=[0,0,0])
    cmy= fload_halo_xyz('y',center=[0,0,0])
    cmz= fload_halo_xyz('z',center=[0,0,0])
    cmm= fload_halo_mass(1)
    cmids= fload_halo_id(1)
    n_ptmass= fload_npart(1)
    id_ptmass1= min(cmids)
    if n_ptmass gt 1 then id_ptmass2= max(cmids)

    sortid= sort(cmids)
    cmx= cmx(sortid)
    cmy= cmy(sortid)
    cmz= cmz(sortid)
    cmm= cmm(sortid)
    cmids= cmids(sortid)

    imgcenter= [0,0,0]
    if keyword_set(pt1center) then imgcenter= [cmx(0), cmy(0), cmz(0)]
    if keyword_set(pt2center) and n_ptmass gt 1 then imgcenter= [cmx(1), cmy(1), cmz(1)]

    cmx= cmx - imgcenter(0)
    cmy= cmy - imgcenter(1)
    cmz= cmz - imgcenter(2)



    ; draw path of each center
    ; --------------------------
    ; center 1
    if not keyword_set(remove_traj1) then begin
	read_center, time, center, filename=frun+"/center_"+strcompress(string(cmids(0)),/remove_all)+".txt"
	;oplot, center[0,*]-imgcenter(0), center[1,*]-imgcenter(1), psym=-3, color= 0, thick= 0.5, linestyle= 2
	oplot, center[0,*]-imgcenter(0), center[1,*]-imgcenter(1), psym=3, color= 0, thick= 0.5
    endif

    ; center 2
    if not keyword_set(remove_traj2) and n_ptmass gt 1 then begin
    	read_center, time, center, filename=frun+"/center_"+strcompress(string(cmids(1)),/remove_all)+".txt"
    	;oplot, center[0,*]-imgcenter(0), center[1,*]-imgcenter(1), psym=-3, color= 0, thick= 0.5, linestyle= 2
    	oplot, center[0,*]-imgcenter(0), center[1,*]-imgcenter(1), psym=3, color= 0, thick= 0.5
    endif


    ;
    ; current point in time
    ;
    symsize= 2.0
    usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
    if not keyword_set(remove_pt1) then oplot, [cmx(0)], [cmy(0)], psym=8, color= 0
    if not keyword_set(remove_pt2) and n_ptmass gt 1 then oplot, [cmx(1)], [cmy(1)], psym=8, color= 0



    ; disk test particles
    ; ---------------------
    xdisk= fload_allstars_xyz('x',center=[0,0,0])
    ydisk= fload_allstars_xyz('y',center=[0,0,0])
    zdisk= fload_allstars_xyz('z',center=[0,0,0])
    mdisk= fload_allstars_mass(1)
    idsdisk= fload_allstars_ids(1)
    n_disk= fload_npart(2)
    print, "n_disk, max/min ID: ", n_disk, max(idsdisk), min(idsdisk)

    sortid= sort(idsdisk)
    xdisk= xdisk(sortid)
    ydisk= ydisk(sortid)
    zdisk= zdisk(sortid)
    mdisk= mdisk(sortid)
    idsdisk= idsdisk(sortid)

    xdisk= xdisk - imgcenter(0)
    ydisk= ydisk - imgcenter(1)
    zdisk= zdisk - imgcenter(2)

    ;symsize= 0.6
    ;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
    ;oplot, xdisk, ydisk, psym=8, color= 150


    nrings= 5
    ringn= fltarr(nrings)
    ringn[0]= 12
    ringn[1]= 26
    ringn[2]= 43
    ringn[3]= 63
    ringn[4]= 87


    fload_newcolortable, 2
    zeroclr= 50
    dclr= 40


    ;
    ;   disk 1
    ;-------------------------------
    if not keyword_set(remove_disk1) then begin
      for i= 1, nrings do begin

	if i eq 1 then lowid= id_ptmass1+1 else lowid= total(ringn[0:i-2])+id_ptmass1+1
	;if i eq nrings then hiid= n_disk else hiid= lowid + ringn[i-1] - 1
	hiid= lowid + ringn[i-1]
	idx= where((idsdisk ge lowid) and (idsdisk) lt hiid)

	ixdisk= xdisk(idx)
	iydisk= ydisk(idx)

	ixdisk= [ixdisk,ixdisk[0]]
	iydisk= [iydisk,iydisk[0]]

	irdisk= sqrt(ixdisk*ixdisk + iydisk*iydisk)

	print, "low/hi= ", lowid, hiid, n_elements(idx), "   min/max r= ", min(irdisk), max(irdisk)

	plot_one_disk, ixdisk, iydisk, zeroclr + (i-1)*dclr

      endfor
    endif



    ;
    ;   disk 2
    ;-------------------------------
    if (n_ptmass gt 1) and (not keyword_set(remove_disk2)) then begin
      for i= 1, nrings do begin

        if i eq 1 then lowid= id_ptmass2+1 else lowid= total(ringn[0:i-2])+id_ptmass2+1
        hiid= lowid + ringn[i-1]
        idx= where((idsdisk ge lowid) and (idsdisk) lt hiid)

        print, "low/hi= ", lowid, hiid, n_elements(idx)

        ixdisk= xdisk(idx)
        iydisk= ydisk(idx)

        ixdisk= [ixdisk,ixdisk[0]]
        iydisk= [iydisk,iydisk[0]]

	plot_one_disk, ixdisk, iydisk, zeroclr + (i-1)*dclr

      endfor
    endif






; -------------
;  Done
; -------------



end




