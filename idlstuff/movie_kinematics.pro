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

   ;starti= 0
   ;starti= 1
   ;starti= 35
   starti= 94

   ;endi= 3
   ;endi= 11
   endi= 97
   ;spawn, "/bin/ls "+frun+"/snap* | wc ",result
   ;endi=long(result[0])-1

   print, "frun= ", frun
   print, "starti= ", starti
   print, "endi= ", endi

   spawn, "mkdir "+frun+"/movie_kin"

   for i=starti,endi do begin

	thisi= i

	exts='0000'+strcompress(string(thisi),/remove_all)

	thisfile= frun+'/movie_kin/img_gal1_0_0_'+strmid(exts,strlen(exts)-3,3)+'.jpg'

	;
	;  need to .run kinematics_info
	;
	;
	pa= 0.0  ; needed to reset position angle
	do_one_angle, frun, 0.0, 0.0, snapnum=thisi, filename="temp_kin.eps", xlen=7.0, $
                                        R_e= R_e, $
                                        Vrot_maj=Vrot_maj, avgVrot_maj=avgVrot_maj, $
                                        Vrot_min=Vrot_min, avgVrot_min=avgVrot_min, $
                                        pa=pa, Sig=sig, $
                                        fit_a=fit_a, fit_b=fit_b, $
                                        fit_ellip=fit_ellip, a4diva=a4diva




	cmd= "convert -quality 95 temp_kin.eps "+thisfile
	spawn, cmd

   endfor

end












;==================================================================================
;
;
;
;

