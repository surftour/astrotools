
pro testit, junk

;frun= "data1/tides/toomrepro1"
;frun= "data1/tides/toomrepro2"
;frun= "data1/tides/toomrepro3"
;frun= "data1/tides/toomrepro4"

frun= "data1/tides/toomredisk_dv_pro2"

;snapnum=100
snapnum=300
;snapnum=375
;snapnum=500

thistempfile= "temp.eps"

thisxlen= 50.0

;toomre_image, frun, snapnum, filename=thistempfile, xlen= thisxlen
;toomre_image, frun, snapnum, filename=thistempfile, xlen= thisxlen, /pt1center
;toomre_image, frun, snapnum, filename=thistempfile, xlen= thisxlen, /pt1center, /remove_disk2
;toomre_image, frun, snapnum, filename=thistempfile, xlen= thisxlen, /pt1center, /remove_disk2, /remove_traj
;toomre_image, frun, snapnum, filename=thistempfile, xlen= thisxlen, /pt2center, /remove_disk1, /remove_traj

toomre_image, frun, snapnum, filename=thistempfile, xlen= thisxlen, /pt1center

end

;=================================================================



pro do_seq, frun

   starti= 0
   ;starti= 1
   ;starti= 6
   ;starti= 51
   ;starti= 101

   ;spawn, "/bin/ls "+frun+"/snap* | wc ",result
   ;endi=long(result[0])-1
   ;endi= 1
   ;endi= 5
   ;endi= 50
   ;endi= 100
   endi= 200
   ;endi= 201
   ;endi= 1500
   ;endi= 5000

   thisxlen= 21.0
   ;thisxlen= 52.5     ; sidelen= 150 kpc
   ;thisxlen= 60.0
   ;thisxlen= 84.0

   print, "frun= ", frun
   print, "starti= ", starti
   print, "endi= ", endi

   spawn, "mkdir "+frun+"/movie_toomre"

   for i=starti,endi do begin


        exts='0000'+strcompress(string(i),/remove_all)


        thisi= i
	;thistempfile= "temp.eps"
	;thistempfile= "temp1.eps"
	;thistempfile= "temp2.eps"
	;thistempfile= "temp3.eps"
	thistempfile= "temp4.eps"

	;toomre_image, frun, thisi, filename=thistempfile, xlen= thisxlen
	;toomre_image, frun, thisi, filename=thistempfile, xlen= thisxlen, /pt1center
	toomre_image, frun, thisi, filename=thistempfile, xlen= thisxlen, /pt1center, /remove_traj1, /remove_traj2

        thisfile= frun+'/movie_toomre/img_'+strmid(exts,strlen(exts)-4,4)+'.jpg'
        cmd= "convert -antialias -quality 95 "+thistempfile+" "+thisfile
        spawn, cmd

   endfor

end




;=================================================================




