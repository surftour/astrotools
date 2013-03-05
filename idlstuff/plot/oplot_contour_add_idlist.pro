;
;
;===============================================
pro oplot_contour_add_idlist, idlistfilename, center=center, $
				rotate_theta=rotate_theta, rotate_phi=rotate_phi
;
;
;
;
;
;


if not keyword_set(idlistfilename) then begin
   print, "  "
   print, "PROBLEM: idlistfilename"
   print, "  "
   print, "  "
   print, "  "
   return
endif

orig_center= center


		; -----------------------
		;
		;   ID's
		; --------------
		; 
		;   gas  only
		;
		; -----------------------

                ;overplot_id_list= 1
                overplot_id_list= 0
                if overplot_id_list eq 1 then begin

			fload_newcolortable, 4

                        ; need to load id_accounting
                        ;idlist_fmfile= load_id_list('id_list.txt')
                        ;idlist_fmfile= load_id_list('id_list_005.txt')
                        ;idlist_fmfile= load_id_list('id_list_018.txt')
			frun= fload_frun(1)
                        ;idlist_fmfile= fload_id_list(frun+'/idlist_suv_v1_029.txt')
                        ;idlist_fmfile= fload_id_list(frun+'/idlist_suv_v2_029.txt')
                        idlist_fmfile= fload_id_list_suv(frun+'/idlist_suv_jage1_030.txt')
                        ;idlist_fmfile= fload_id_list_suv(frun+'/idlist_suv_jage2_030.txt')
                        ;idlist_fmfile= fload_id_list(frun+'/id_list_test.txt')
                        print, "id's in file= ",n_elements(idlist_fmfile)
                        sidx=sort(idlist_fmfile)
                        idlist_fmfile= idlist_fmfile(sidx)

                        ; gas
                        ; -----
			if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin
                        	x= fload_gas_xyz('x',center=orig_center)
                        	y= fload_gas_xyz('y',center=orig_center)
                        	z= fload_gas_xyz('z',center=orig_center)
				; new rotation procedure
				process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new
				x= x_new
				y= y_new
				z= z_new
			endif else begin
                        	;x= fload_gas_xyz('x',center=[0,0,0])
                        	;y= fload_gas_xyz('y',center=[0,0,0])
                        	;z= fload_gas_xyz('z',center=[0,0,0])
                        	x= fload_gas_xyz('x',center=orig_center)
                        	y= fload_gas_xyz('y',center=orig_center)
                        	z= fload_gas_xyz('z',center=orig_center)
			endelse
                        gid= fload_gas_id(1)
                        ;idx= where(idlist_fmfile eq gid)
                        idx= intarr(n_elements(gid))
                        lstid= -1
                        duplicates= 0
                        for i=0L,n_elements(idlist_fmfile)-1 do begin
                            if idlist_fmfile[i] eq lstid then duplicates= duplicates+1
                            inlst= where(gid eq idlist_fmfile[i])
                            if inlst(0) ne -1 then begin
                                idx(inlst)= 1
                            endif else begin
                                ;print, "couldn't find id=",idlist_fmfile[i]
                            endelse
                            lstid= idlist_fmfile[i]
                        endfor
                        midx= where(idx eq 1)
                        print, "matching gid's= ",n_elements(midx)
                        print, "duplicate id's= ",duplicates
                        if midx(0) ne -1 then oplot, x(midx), y(midx), psym=7, color= 150, symsize= 0.2
                        ;if idx(0) ne -1 then oplot, x(midx), y(midx), psym=3, color= 150

                endif


		; -----------------------
		;
		;   ID's
		; --------------
		; 
		; gas + new stars
		;
		; -----------------------

                ;overplot_id_list= 1
                overplot_id_list= 0
                if overplot_id_list eq 1 then begin

			fload_newcolortable, 4

                        ; need to load id_accounting
                        ;idlist_fmfile= load_id_list('id_list.txt')
                        ;idlist_fmfile= load_id_list('id_list_005.txt')
                        ;idlist_fmfile= load_id_list('id_list_018.txt')
			frun= fload_frun(1)
                        ;idlist_fmfile= fload_id_list(frun+'/idlist_suv_v1_029.txt')
                        ;idlist_fmfile= fload_id_list(frun+'/idlist_suv_v2_029.txt')
                        ;idlist_fmfile= fload_id_list_suv(frun+'/idlist_suv_jage1_030.txt')
                        idlist_fmfile= fload_id_list_suv(frun+'/idlist_suv_jage2_030.txt')
                        ;idlist_fmfile= fload_id_list(frun+'/id_list_test.txt')
                        print, "id's in file= ",n_elements(idlist_fmfile)
                        sidx=sort(idlist_fmfile)
                        idlist_fmfile= idlist_fmfile(sidx)

                        ; gas
                        ; -----
			if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin
                        	x= fload_gas_xyz('x',center=orig_center)
                        	y= fload_gas_xyz('y',center=orig_center)
                        	z= fload_gas_xyz('z',center=orig_center)
				; new rotation procedure
				process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new
				x= x_new
				y= y_new
				z= z_new
			endif else begin
                        	;x= fload_gas_xyz('x',center=[0,0,0])
                        	;y= fload_gas_xyz('y',center=[0,0,0])
                        	;z= fload_gas_xyz('z',center=[0,0,0])
                        	x= fload_gas_xyz('x',center=orig_center)
                        	y= fload_gas_xyz('y',center=orig_center)
                        	z= fload_gas_xyz('z',center=orig_center)
			endelse
                        gid= fload_gas_id(1)
                        ;idx= where(idlist_fmfile eq gid)
                        idx= intarr(n_elements(gid))
                        lstid= -1
                        duplicates= 0
                        for i=0L,n_elements(idlist_fmfile)-1 do begin
                            if idlist_fmfile[i] eq lstid then duplicates= duplicates+1
                            inlst= where(gid eq idlist_fmfile[i])
                            if inlst(0) ne -1 then begin
                                idx(inlst)= 1
                            endif else begin
                                ;print, "couldn't find id=",idlist_fmfile[i]
                            endelse
                            lstid= idlist_fmfile[i]
                        endfor
                        midx= where(idx eq 1)
                        print, "matching gid's= ",n_elements(midx)
                        print, "duplicate id's= ",duplicates
                        if idx(0) ne -1 then oplot, x(midx), y(midx), psym=7, color= 150, symsize= 0.2
                        ;if idx(0) ne -1 then oplot, x(midx), y(midx), psym=3, color= 150

                        ; new stars
                        ; ---------
			if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin
                        	x= fload_newstars_xyz('x',center=orig_center)
                        	y= fload_newstars_xyz('y',center=orig_center)
                        	z= fload_newstars_xyz('z',center=orig_center)
				; new rotation procedure
				process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new
				x= x_new
				y= y_new
				z= z_new
			endif else begin
                        	x= fload_newstars_xyz('x',center=[0,0,0])
                        	y= fload_newstars_xyz('y',center=[0,0,0])
                        	z= fload_newstars_xyz('z',center=[0,0,0])
			endelse
                        nsid= fload_newstars_id(1)
                        ;pid= fload_newstars_parentid(1)
                        ;idx= where(idlist_fmfile eq nsid)
                        idx= intarr(n_elements(nsid))
                        for i=0L,n_elements(idlist_fmfile)-1 do begin
                            inlst= where(nsid eq idlist_fmfile[i])
                            if inlst(0) ne -1 then idx(inlst)= 1
                            ;inlst= where(pid eq idlist_fmfile[i])
                            ;if inlst(0) ne -1 then idx(inlst)= 1
                        endfor
                        midx=where(idx eq 1)
                        if midx(0) ne -1 then begin
                           print, "matching nsid's= ",n_elements(midx)
                           ;if idx(0) ne -1 then oplot, x(midx), y(midx), psym=7, color= 220
                           if idx(0) ne -1 then oplot, x(midx), y(midx), psym=3, color= 220
                        endif else begin
                           print, "no matching nsid's"
                        endelse
                endif


                ; --------------------------------------------------------------------
		;
		;    ID's
		;  --------
		;
		;   Stars:  Old 
		;
                ; --------------------------------------------------------------------

                ;overplot_oldid_list= 1
                overplot_oldid_list= 0
                if overplot_oldid_list eq 1 then begin

			fload_newcolortable, 4

                        ; need to load id_accounting
                        ;idlist_fmfile= load_id_list('id_list.txt')
                        ;idlist_fmfile= load_id_list('id_list_001.txt')
                        idlist_fmfile= load_id_list('id_list_005.txt')
                        ;idlist_fmfile= load_id_list('id_list_018.txt')
                        print, "id's in file= ",n_elements(idlist_fmfile)
                        sidx=sort(idlist_fmfile)
                        idlist_fmfile= idlist_fmfile(sidx)

                        ; oldstars
                        ; ----------
                        x= fload_oldstars_xyz('x',center=orig_center)
                        y= fload_oldstars_xyz('y',center=orig_center)
                        gid= fload_oldstars_id(1)
                        ;idx= where(idlist_fmfile eq gid)
                        idx= intarr(n_elements(gid))
                        lstid= -1
                        duplicates= 0
                        for i=0,n_elements(idlist_fmfile)-1 do begin
                            if idlist_fmfile[i] eq lstid then duplicates= duplicates+1
                            inlst= where(gid eq idlist_fmfile[i])
                            if inlst(0) ne -1 then begin
                                idx(inlst)= 1
                            endif else begin
                                ;print, "couldn't find id=",idlist_fmfile[i]
                            endelse
                            lstid= idlist_fmfile[i]
                        endfor
                        midx= where(idx eq 1)
                        print, "matching gid's= ",n_elements(midx)
                        print, "duplicate id's= ",duplicates
                        if idx(0) ne -1 then oplot, x(midx), y(midx), psym=7, color= 150, symsize= 0.2
                        ;if idx(0) ne -1 then oplot, x(midx), y(midx), psym=3, color= 150

                endif


                ; --------------------------------------------------------------------
		;
		;    ID's
		;  --------
		;
		;   Stars:  All of 'em
		;
                ; --------------------------------------------------------------------

                overplot_oldid_list= 1
                ;overplot_oldid_list= 0
                if overplot_oldid_list eq 1 then begin

			fload_newcolortable, 4

                        ; need to load id_accounting
                        ;idlist_fmfile= fload_id_list('/raid4/tcox/localgroup/v4/earth_idlist.txt')
                        ;idlist_fmfile= fload_id_list('/raid4/tcox/localgroup/bhires/earth_idlist.txt')
                        ;idlist_fmfile= fload_id_list('/n/scratch/hernquist_lab/tcox/tests/Sb10xSb10xhs_e_256/tidal_tail_idlist.txt')
                        idlist_fmfile= fload_id_list(idlistfilename)
			;frun= fload_frun(1)
                        ;idlist_fmfile= fload_id_list(frun+'/id_list_test.txt')
                        print, "id's in file= ",n_elements(idlist_fmfile)
                        sidx=sort(idlist_fmfile)
                        idlist_fmfile= idlist_fmfile(sidx)

                        ; allstars
                        ; ----------
                        if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin
                                x= fload_allstars_xyz('x',center=orig_center)
                                y= fload_allstars_xyz('y',center=orig_center)
                                z= fload_allstars_xyz('z',center=orig_center)
                                ; new rotation procedure
                                process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new
                                x= x_new
                                y= y_new
                                z= z_new
                        endif else begin
                                ;x= fload_allstars_xyz('x',center=[0,0,0])
                                ;y= fload_allstars_xyz('y',center=[0,0,0])
                                ;z= fload_allstars_xyz('z',center=[0,0,0])
                                x= fload_allstars_xyz('x',center=orig_center)
                                y= fload_allstars_xyz('y',center=orig_center)
                                z= fload_allstars_xyz('z',center=orig_center)
                        endelse
                        gid= fload_allstars_ids(1)
                        ;idx= where(idlist_fmfile eq gid)
                        idx= intarr(n_elements(gid))
                        lstid= -1
                        duplicates= 0
                        for i=0L,n_elements(idlist_fmfile)-1 do begin
                            if idlist_fmfile[i] eq lstid then begin
				duplicates= duplicates+1
				print, lstid, idlist_fmfile[i]
			    endif
                            inlst= where(gid eq idlist_fmfile[i])
                            if inlst(0) ne -1 then begin
                                idx(inlst)= 1 
                            endif else begin
                                ;print, "couldn't find id=",idlist_fmfile[i]
                            endelse
                            lstid= idlist_fmfile[i]
                        endfor
                        midx= where(idx eq 1)
                        print, "matching id's= ",n_elements(midx)
                        print, "duplicate id's= ",duplicates
                        ;if midx(0) ne -1 then oplot, x(midx), y(midx), psym=7, color= 150, symsize= 0.4
                        if midx(0) ne -1 then oplot, x(midx), y(midx), psym=3, color= 150 , symsize= 1.0

                endif



		; --------------------------------------------------------------------


		overplot_hotmetalgas_list= 0
		;overplot_hotmetalgas_list= 1
		if overplot_hotmetalgas_list eq 1 then begin

			fload_newcolortable, 4

                        ; need to load id_accounting
			;frun="/raid4/tcox/vc3vc3e"
			;frun="/raid4/tcox/vc3vc3e_no"
			;frun="/raid4/tcox/vc3vc3e_2"
			;frun="/raid4/tcox/vc3vc3k"
			idlistfilename= frun+"/hotgas_idlist.txt"
			;idlistname= frun+"/wind_idlist.txt"
			;idlistname= frun+"/wind_idlist_atfirstpassage.txt"
			;idlistname= frun+"/id_list_ns_j_p_in.txt"
			;idlistname= frun+"/id_list_ns_j_p_out.txt"
			;idlistname= frun+"/id_list_ns_j_m_in.txt"
			;idlistname= frun+"/id_list_ns_j_m_out.txt"
			;idlistname= frun+"/id_list_os_j_p_in.txt"
			;idlistname= frun+"/id_list_os_j_p_out.txt"
			;idlistname= frun+"/id_list_os_j_m_in.txt"
			;idlistname= frun+"/id_list_os_j_m_out.txt"
			;idlistfilename= "/home/tcox/Sbc10x/idlists/group_3_idlist.txt"
			;idlist_fmfile= load_id_list(idlistname)
			; --------------------------------------
			spawn, 'wc '+idlistfilename, result
			lines= long(result)
			idlist= lonarr(lines-5)

			get_lun, unit
			openr, unit, idlistfilename

			textjunk= ''
			readf, unit, textjunk
			readf, unit, textjunk
			readf, unit, textjunk
			readf, unit, textjunk
			readf, unit, textjunk
			readf, unit, idlist
			close, unit
			idlist_fmfile= idlist
			; --------------------------------------
			print, "id's in file= ",n_elements(idlist_fmfile)
			sidx=sort(idlist_fmfile)
			idlist_fmfile= idlist_fmfile(sidx)

			; gas particles
			; -------------
                        x= fload_gas_xyz('x',center=center)
                        y= fload_gas_xyz('y',center=center)
                        gid= fload_gas_id(1)

                        idx= intarr(n_elements(gid))
                        for i=0L,n_elements(idlist_fmfile)-1 do begin
                            inlst= where(gid eq idlist_fmfile[i])
                            if inlst(0) ne -1 then begin
                                idx(inlst)= 1
                            endif else begin
                                print, "couldn't find id=",idlist_fmfile[i]
                            endelse
                        endfor
                        midx= where(idx eq 1)
                        print, "matching gid's= ",n_elements(midx)
                        if midx(0) ne -1 then oplot, x(midx), y(midx), psym=7, color= 150, symsize= 0.2
                        ;if midx(0) ne -1 then oplot, x(midx), z(midx), psym=7, color= 150, symsize= 0.2


		endif


		; --------------------------------------------------------------------

		overplot_tdwarf_list= 0
		;overplot_tdwarf_list= 1
		if overplot_tdwarf_list eq 1 then begin

			fload_newcolortable, 4

			tdwarfs_overplot_idlists, "/home/tcox/Sbc10x/idlists/group_3_idlist.txt", ptcolor=150
			tdwarfs_overplot_idlists, "/home/tcox/Sbc10x/idlists/group_5_idlist.txt", ptcolor=200
			;tdwarfs_overplot_idlists, "/home/tcox/Sbc10x/idlists/group_10_idlist.txt", ptcolor=100
			tdwarfs_overplot_idlists, "/home/tcox/Sbc10x/idlists/group_6_idlist.txt", ptcolor=100
			tdwarfs_overplot_idlists, "/home/tcox/Sbc10x/idlists/group_15_idlist.txt", ptcolor=50

		endif


		; --------------------------------------------------------------------



; -------------
;  Done
; -------------



end


