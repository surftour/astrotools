;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;  
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro do_seq, frun

   starti= 0
   ;starti= 1
   ;starti= 35
   ;starti= 94

   ;endi= 3
   ;endi= 11
   ;endi= 97
   endi= 999
   ;spawn, "/bin/ls "+frun+"/snap* | wc ",result
   ;endi=long(result[0])-1


   ; controlled kick test
   frun= "/n/data1/ITC_lab/lblecha/control_kick_test/v0_t0/"
   epsfile= "temp_laura2.eps"
   ; entire system just sits there
   ;frun= "/n/data1/ITC_lab/lblecha/kick_1gal/debugged_output/resx2_softdv4_dmonly_vk0"
   ;epsfile= "temp_laura.eps"


   print, "frun= ", frun
   print, "starti= ", starti
   print, "endi= ", endi

   ;spawn, "mkdir "+frun+"/movie_kin"
   spawn, "mkdir "+frun+"/movie_laura"

   for i=starti,endi do begin

	thisi= i

	exts='0000'+strcompress(string(thisi),/remove_all)

	;thisfile= frun+'/movie_kin/img_gal1_0_0_'+strmid(exts,strlen(exts)-3,3)+'.jpg'
	;thisfile= '/n/home/tcox/movie_laura/img_'+strmid(exts,strlen(exts)-3,3)+'.jpg'
	thisfile= '/n/home/tcox/movie_laura/img_v0t0_'+strmid(exts,strlen(exts)-3,3)+'.jpg'

	;
	;
	;
        contour_dm, frun, thisi, 50.0, 'ps', filename=epsfile, /showbhs, center=[0, 0,0]

	cmd= "convert -quality 95 "+epsfile+" "+thisfile
	spawn, cmd

   endfor

end












;==================================================================================
;
;
;
;

