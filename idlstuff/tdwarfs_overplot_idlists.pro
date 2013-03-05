pro tdwarfs_overplot_idlists, idfname, ptcolor=ptcolor



                        idlistfilename= idfname
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
                        if midx(0) ne -1 then oplot, x(midx), y(midx), psym=7, color= ptcolor, symsize= 0.2
                        ;if midx(0) ne -1 then oplot, x(midx), z(midx), psym=7, color= 150, symsize= 0.2


end
