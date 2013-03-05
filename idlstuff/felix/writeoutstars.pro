PRO writeoutstars,namein,numin,name,num,center,$
                  outstarspos,outstarsvel,outnewstarspos,outnewstarsvel,$
		  meanreff

     IF ( (name eq namein) and (num eq numin)) THEN BEGIN
        openw,3,'vp.'+strcompress(name,/remove_all)+'.dat'
	FOR dir=0,2 DO BEGIN     
	   getnewcenter,center,dir,num,xxx,yyy,zzz,xcenter,ycenter,zcenter	             

	   distance2D = sqrt((outstarspos(xxx,*)-xcenter)^2 + (outstarspos(yyy,*)-ycenter)^2)
	   vel2D      = outstarsvel(zzz,*)
	   
	   FOR i=0L,n_elements(outstarspos(0,*))-1 DO BEGIN
	      printf,3,distance2D(i)/meanreff,vel2D(i)
	   ENDFOR
	ENDFOR
        close,3	    	

        openw,3,'vp.new.'+strcompress(name,/remove_all)+'.dat'
	FOR dir=0,2 DO BEGIN     
	   getnewcenter,center,dir,num,xxx,yyy,zzz,xcenter,ycenter,zcenter	             

	   distance2D = sqrt((outnewstarspos(xxx,*)-xcenter)^2 + (outnewstarspos(yyy,*)-ycenter)^2)
	   vel2D      = outnewstarsvel(zzz,*)
	   
	   FOR i=0L,n_elements(outnewstarspos(0,*))-1 DO BEGIN
	      printf,3,distance2D(i)/meanreff,vel2D(i)
	   ENDFOR
	ENDFOR
        close,3	    	
	
	
     ENDIF

END
