

;=================================================================



pro do_seq, frun

   starti= 0
   ;starti= 1
   ;starti= 6
   ;starti= 51
   starti= 101

   spawn, "/bin/ls "+frun+"/snap* | wc ",result
   endi=long(result[0])-1
   ;endi= 5
   ;endi= 50
   ;endi= 100
   endi= 201

   print, "frun= ", frun
   print, "starti= ", starti
   print, "endi= ", endi

   spawn, "mkdir "+frun+"/movie_kinematics"

   for i=starti,endi do begin


        exts='0000'+strcompress(string(i),/remove_all)


        thisi= i
	;thistempfile= "temp.eps"
	thistempfile= "temp1.eps"
	do_one_angle, frun, 0.0, 0.0, snapnum= thisi, filename=thistempfile, xlen= 10.0, reff=3.0
        thisfile= frun+'/movie_kinematics/img_'+strmid(exts,strlen(exts)-3,3)+'.jpg'
        cmd= "convert -quality 95 "+thistempfile+" "+thisfile
        spawn, cmd

        ;thisi= i
	;do_one_angle, frun, 90.0, 0.0, snapnum= thisi, filename="temp_90.eps", xlen= 10.0, reff=3.0
        ;thisfile= frun+'/movie_kinematics/img_90_'+strmid(exts,strlen(exts)-3,3)+'.jpg'
        ;cmd= "convert -quality 95 temp_90.eps "+thisfile
        ;spawn, cmd

        ;thisi= i
	;do_one_angle, frun, 90.0, 90.0, snapnum= thisi, filename="temp_9090.eps", xlen= 10.0, reff=3.0
        ;thisfile= frun+'/movie_kinematics/img_9090_'+strmid(exts,strlen(exts)-3,3)+'.jpg'
        ;cmd= "convert -quality 95 temp_9090.eps "+thisfile
        ;spawn, cmd

   endfor

end




;=================================================================




