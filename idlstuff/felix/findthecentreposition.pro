PRO findthecentreposition,pos,mass,num,N,center,centerpositionfile

  if (file_test(centerpositionfile) eq 1) then begin    
     xc = 0.0
     yc = 0.0
     zc = 0.0
     openr,1,centerpositionfile
     readf,1,xc,yc,zc
     center    = fltarr(3);     
     center(0) = xc
     center(1) = yc
     center(2) = zc     
     close,1
  endif else begin   


     print,"Computing the centre position ..."
     minnumber = 100
     rmax      = 100.0

     center    = fltarr(3);

     FOR i=0,2 DO BEGIN
	center(i) = total(pos(i,*)*mass(*),/double)/total(mass(*),/double);
     ENDFOR   

     distance = (pos(0,*)-center(0)) * (pos(0,*)-center(0)) + $
        	(pos(1,*)-center(1)) * (pos(1,*)-center(1)) + $
        	(pos(2,*)-center(2)) * (pos(2,*)-center(2))

     count = N

     WHILE (count gt minnumber) DO BEGIN
	; print,count,center(0),center(1),center(2),sqrt(rmax)
	rmax = 0.8 * rmax

	select = where(distance lt rmax,count)

	oldcenter = center
	FOR i=0,2 DO BEGIN
	   center(i) = total(pos(i,select)*mass(select),/double)/total(mass(select),/double)
	ENDFOR   

	distance = (pos(0,*)-center(0)) * (pos(0,*)-center(0)) + $
        	   (pos(1,*)-center(1)) * (pos(1,*)-center(1)) + $
        	   (pos(2,*)-center(2)) * (pos(2,*)-center(2))

	;; VERIFICATION: OK
	;oplot,[oldcenter(0),center(0)],[oldcenter(1),center(1)],color=250,thick=3;
	;xcirc = sqrt(rmax) * cos(lindgen(101)/100.0*2.0*!Pi)
	;ycirc = sqrt(rmax) * sin(lindgen(101)/100.0*2.0*!Pi)

       ; oplot,[oldcenter(0),center(0)],[oldcenter(1),center(1)],color=50,thick=0;     
       ; oplot,center(0)+xcirc,center(1)+ycirc,color=50,thick=0;          
       ENDWHILE
    openw,1,centerpositionfile
    printf,1,center(0),center(1),center(2)
    close,1
       
  endelse
  
END
