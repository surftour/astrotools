pro oplot_bh_separation, x0, y0, xz=xz, projonly=projonly

                        bhid= fload_blackhole_id(1)
			if n_elements(bhid) lt 2 then return

                        ; bh 1
                        bhidlbl1= strcompress(string(bhid[0]),/remove_all)
                        read_center, ttime1, center1, filename=fload_frun(1)+"/centers_bh_"+bhidlbl1+".txt"
                        idx= where(ttime1 ge fload_time(1)-1.0e-5)
                        bhcen1= center1(*,idx(0))

                        ; bh 2
                        bhidlbl2= strcompress(string(bhid[1]),/remove_all)
                        read_center, ttime2, center2, filename=fload_frun(1)+"/centers_bh_"+bhidlbl2+".txt"
                        idx= where(ttime2 ge fload_time(1)-1.0e-5)
                        bhcen2= center2(*,idx(0))


                        rdiff= bhcen1-bhcen2
                        rsep= sqrt(rdiff(0)*rdiff(0) + rdiff(1)*rdiff(1) + rdiff(2)*rdiff(2))
                        rprojsep= sqrt(rdiff(0)*rdiff(0) + rdiff(1)*rdiff(1))
                        if keyword_set(xz) then rprojsep= sqrt(rdiff(0)*rdiff(0) + rdiff(2)*rdiff(2))

                        rsep=rsep/0.7
                        rprojsep=rprojsep/0.7
                        print, "separations: ", rsep, rprojsep

	if not keyword_set(projonly) then begin
              rseplbl= strcompress(string(rsep),/remove_all)
	      digs= 1
              if (rsep) ge 1.0 then digs= 1
              if (rsep) ge 10.0 then digs= 2
              rseplbl = strmid(rseplbl,0,digs+2)        ; T=0.x (4+digits after decimal)
              rseplbl= '!7D!6R= '+rseplbl+' kpc'
              ;xyouts, 0.07, 0.90, rseplbl, /normal, size= 1.4, charthick=3.0, color= 0
              xyouts, x0, y0, rseplbl, /normal, size= 1.4, charthick=3.0, color= 0
	endif

                        rsepprojlbl= strcompress(string(rprojsep),/remove_all)
			digs= 1
                        if (rsep) ge 1.0 then digs= 1
                        if (rsep) ge 10.0 then digs= 2
                        rsepprojlbl = strmid(rsepprojlbl,0,digs+2)        ; T=0.x (4+digits after decimal)
                        rsepprojlbl= '!7D!6R!Dproj!N= '+rsepprojlbl+' kpc'
                        xyouts, x0, y0-0.025, rsepprojlbl, /normal, size= 1.4, charthick=3.0, color= 0


end


