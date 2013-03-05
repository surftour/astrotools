PRO fitromanowsky,ysigma,xsigma,xstars,ystars,thick=thick

   interpolationpoint=0.5

   symbolsize = 1.5
   symsize = 1.5
   colorstouse=[6,7,8,9]
   minxtouse  = [0.2,0.2,0.2,0.2]

   scale = 0.0
   for i=3,0,-1 do begin
   
      if (i eq 0) then begin
	 ; read mendez data
	 men = read_ascii("men.dat")

	 ; scale to the effective radius of 73.8 arcsec
	 men.field1(0,*) = men.field1(0,*) / 73.8
      
         scale = 153.0/interpol(ysigma,      xsigma,[interpolationpoint,interpolationpoint])	 
	 scale = 1.0/scale(0)

         n = n_elements(men.field1(0,*))
	 
	 romstars = {field1 : fltarr(6,n-8)}	 
	 rompn    = {field1 : fltarr(6,8)}
	 
	 romstars.field1(0:1,0:n-1-8) = men.field1(0:1,0:n-1-8)
	 rompn.field1(0:1,0:7)        = men.field1(0:1,n-8 + indgen(8))

	 for k=2,5,2 do begin 
	    romstars.field1(k,*) = men.field1(2,0:n_elements(men.field1(0,*))-1-8)
	    rompn.field1(k,*)    = men.field1(2,n_elements(men.field1(0,*))-8 + indgen(8))
         endfor	 
	 for k=3,5,2 do begin 
	    romstars.field1(k,*) = -men.field1(2,0:n_elements(men.field1(0,*))-1-8)
	    rompn.field1(k,*)    = -men.field1(2,n_elements(men.field1(0,*))-8 + indgen(8))
         endfor	 	 
	 
	 romdm    = {field1 : fltarr(2,100)}
	 romnodm  = {field1 : fltarr(2,100)}
      endif else begin 

	 ; read romanowsky data
	 romdm    = read_ascii("rom."+strcompress(string(i),/remove_all)+".dm.dat")
	 romnodm  = read_ascii("rom."+strcompress(string(i),/remove_all)+".nodm.dat")
	 rompn    = read_ascii("rom."+strcompress(string(i),/remove_all)+".pn.dat")      
	 romstars = read_ascii("rom."+strcompress(string(i),/remove_all)+".stars.dat")      	 

	 ; get the old scaling for nostalgy
	 scale = interpol(romdm.field1(1,*),romdm.field1(0,*),[interpolationpoint,interpolationpoint])/$
		 interpol(ysigma,      xsigma,[interpolationpoint,interpolationpoint]) + $
        	 interpol(romnodm.field1(1,*),romnodm.field1(0,*),[interpolationpoint,interpolationpoint])/$
		 interpol(ysigma,      xsigma,[interpolationpoint,interpolationpoint])
	 scale = 2.0/scale(0)
	 
      endelse

      if (i eq 2) then begin
         scalingfactor        = 35.0/54.9
	 rompn.field1(0,*)    = rompn.field1(0,*)    * scalingfactor
	 romstars.field1(0,*) = romstars.field1(0,*) * scalingfactor	 
	 rompn.field1(2,*)    = rompn.field1(2,*)    * scalingfactor
	 rompn.field1(3,*)    = rompn.field1(3,*)    * scalingfactor	 
      endif

      ; making the best fits
      n1 = n_elements(romstars.field1(0,*))
      n2 = n_elements(rompn.field1(0,*))
      n  = n1 + n2
      xxx = fltarr(n)
      yyy = fltarr(n)
      zzz = fltarr(n)
      xxx(0:n1-1) = romstars.field1(0,*)
      xxx(n1:n-1) = rompn.field1(0,*)
      
      yyy(0:n1-1) = romstars.field1(1,*)
      yyy(n1:n-1) = rompn.field1(1,*)
      
      zzz(0:n1-1) = 0.5*(abs(romstars.field1(2,*)) + abs(romstars.field1(3,*)))
      zzz(n1:n-1) = 0.5*(abs(rompn.field1(4,*))    + abs(rompn.field1(5,*)))

      select = where(xxx gt minxtouse(i))
      xxx = xxx(select)
      yyy = yyy(select)
      zzz = zzz(select)
          
      ;; print,"before:", scale
      sss    = interpol(ystars,xstars,xxx)
      
      scale1 = 0
      scale2 = 1
      scale  = 0.5
      value1 = total( ((yyy*scale1-sss)/zzz)^2)
      value2 = total( ((yyy*scale2-sss)/zzz)^2)
      value  =  1e10
      oldval = -1e10
      WHILE (abs(scale1-scale2) gt 0.00001) DO BEGIN
         value  = total( ((yyy*scale-sss)/zzz)^2)	       
         IF (value2 gt value1) THEN BEGIN
	    scale2 = scale
	    value2 = value
	 ENDIF ELSE BEGIN
	    scale1 = scale
	    value1 = value
         ENDELSE
         scale  = 0.5 * (scale1 + scale2)	 
      ENDWHILE      
      print,"after: ", scale

      if (i eq 3) then begin
         oplot,romdm.field1(0,*),romdm.field1(1,*)*scale,color=3,thick=thick
         oplot,romnodm.field1(0,*),romnodm.field1(1,*)*scale,color=3,thick=thick
      endif	 
                  
      ;plotsym,0,/fill;
      usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
      oplot,rompn.field1(0,*),rompn.field1(1,*)*scale,symsize=symbolsize,psym=8,color=colorstouse(i)
      error_bars,rompn.field1(0,*),rompn.field1(1,*)*scale+rompn.field1(4,*)*scale,rompn.field1(1,*)*scale+rompn.field1(5,*)*scale,color=colorstouse(i),thick=thick   

      if (i gt 0) then begin
         ; there are no x-axis error bars for the mendez data
         error_bars,rompn.field1(1,*)*scale,rompn.field1(0,*)+rompn.field1(2,*),rompn.field1(0,*)+rompn.field1(3,*),/rotate,color=colorstouse(i),thick=thick   
      endif

      ;plotsym,3,/fill;
      oplot,romstars.field1(0,*),romstars.field1(1,*)*scale,symsize=symbolsize*0.75,psym=7,color=colorstouse(i)
      error_bars,romstars.field1(0,*),romstars.field1(1,*)*scale+romstars.field1(2,*)*scale,romstars.field1(1,*)*scale+romstars.field1(3,*)*scale,color=colorstouse(i),thick=thick   

      ; CHECK : OK 
      ;oplot,xxx,yyy*scale,psym=4,color=0
      ;error_bars,xxx,(yyy+zzz)*scale,(yyy-zzz)*scale,color=0,thick=thick
     

      if (i eq 2) then begin
         ; draw the arrow
	 shiftit = 0 
	 xshift1  = 0.8
	 xshift2  = -0.225
	 n = 0;n_elements(rompn.field1(0,*))-1
	 ;; arrow,rompn.field1(0,n)+xshift1,rompn.field1(1,n)*scale+shiftit,rompn.field1(0,n)/scalingfactor+xshift2,rompn.field1(1,n)*scale+shiftit,/data,color=colorstouse(i)
         ;plotsym,0,thick=thick
         oplot,[rompn.field1(0,n)/scalingfactor,rompn.field1(0,n)/scalingfactor],$
	       [rompn.field1(1,n)*scale,rompn.field1(1,n)*scale],symsize=symbolsize*0.9,psym=8,color=colorstouse(i)
      endif

      ; print the data for Gary
      openw,1,"stackeddata.stars."+strcompress(string(i),/remove_all)+".dat"
         for j=0,n_elements(romstars.field1(0,*))-1 do begin
	    printf,1,romstars.field1(0,j),romstars.field1(1,j)*scale,romstars.field1(2,j)*scale,romstars.field1(3,j)*scale
	 endfor      
      close,1

      openw,1,"stackeddata.pn."+strcompress(string(i),/remove_all)+".dat"
        for j=0,n_elements(rompn.field1(0,*))-1 do begin
	    if (i eq 0) then begin
	       printf,1,rompn.field1(0,j),rompn.field1(1,j)*scale,rompn.field1(4,j)*scale,rompn.field1(5,j)*scale,0.0,0.0
	    endif else begin
	       printf,1,rompn.field1(0,j),rompn.field1(1,j)*scale,rompn.field1(4,j)*scale,rompn.field1(5,j)*scale,rompn.field1(2,j),rompn.field1(3,j)
	    endelse   
	 endfor      
      close,1
	 	 
   endfor
  
END
