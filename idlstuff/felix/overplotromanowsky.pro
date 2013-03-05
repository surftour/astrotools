PRO overplotromanowsky,ysigma,xsigma,interpolationpoint,thick=thick

   symbolsize = 1.5
   symsize = 1.5
   colorstouse=[6,7,8,9]

   for i=1,3 do begin
      ; read romanowsky data
      romdm    = read_ascii("rom."+strcompress(string(i),/remove_all)+".dm.dat")
      romnodm  = read_ascii("rom."+strcompress(string(i),/remove_all)+".nodm.dat")
      rompn    = read_ascii("rom."+strcompress(string(i),/remove_all)+".pn.dat")      
      romstars = read_ascii("rom."+strcompress(string(i),/remove_all)+".stars.dat")      	 

      if (i eq 2) then begin
         scalingfactor = 35.0/54.9
      endif

      ; overplot the romanowsky data    
      scale = interpol(romdm.field1(1,*),romdm.field1(0,*),[interpolationpoint,interpolationpoint])/$
	      interpol(ysigma,      xsigma,[interpolationpoint,interpolationpoint]) + $
              interpol(romnodm.field1(1,*),romnodm.field1(0,*),[interpolationpoint,interpolationpoint])/$
	      interpol(ysigma,      xsigma,[interpolationpoint,interpolationpoint])
      scale = 2.0/scale(0)

      ;plotsym,0,/fill;
      usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
      oplot,rompn.field1(0,*),rompn.field1(1,*)*scale,symsize=symbolsize,psym=8,color=colorstouse(i)
      error_bars,rompn.field1(0,*),rompn.field1(1,*)*scale+rompn.field1(4,*)*scale,rompn.field1(1,*)*scale+rompn.field1(5,*)*scale,color=colorstouse(i),thick=thick   

      error_bars,rompn.field1(1,*)*scale,rompn.field1(0,*)+rompn.field1(2,*),rompn.field1(0,*)+rompn.field1(3,*),/rotate,color=colorstouse(i),thick=thick   

      ;plotsym,3,/fill;
      oplot,romstars.field1(0,*),romstars.field1(1,*)*scale,symsize=symbolsize*0.75,psym=7,color=colorstouse(i)
      error_bars,romstars.field1(0,*),romstars.field1(1,*)*scale+romstars.field1(2,*)*scale,romstars.field1(1,*)*scale+romstars.field1(3,*)*scale,color=colorstouse(i),thick=thick   

      if (i eq 2) then begin
         ; draw the arrow
	 shiftit = 0 
	 xshift1  = -0.92
	 xshift2  = 0.2
	 n = 0;n_elements(rompn.field1(0,*))-1
	 arrow,rompn.field1(0,n)+xshift1,rompn.field1(1,n)*scale+shiftit,rompn.field1(0,n)*scalingfactor+xshift2,rompn.field1(1,n)*scale+shiftit,/data,color=colorstouse(i)
         ;plotsym,0,thick=thick
	 usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
         oplot,[rompn.field1(0,n)*scalingfactor,rompn.field1(0,n)*scalingfactor],$
	       [rompn.field1(1,n)*scale,rompn.field1(1,n)*scale],symsize=symbolsize*0.9,psym=8,color=colorstouse(i)
      endif

      if (i eq 1) then begin
         oplot,romdm.field1(0,*),romdm.field1(1,*)*scale,color=3,thick=thick
         oplot,romnodm.field1(0,*),romnodm.field1(1,*)*scale,color=3,thick=thick
      endif	 
	 	 
   endfor

   ; read mendez data
   men = read_ascii("men.old.dat")
   
   ; scale to the effective radius of 91.0 arcsec
   men.field1(0,*) = men.field1(0,*) / 91.0

   scale = 153.0/$
	  interpol(ysigma,      xsigma,[interpolationpoint,interpolationpoint])
   scale = 1.0/scale(0);;

   ; mendez PN
   selectvalues = where(men.field1(2,*) gt 0)
   ;plotsym,0,/fill;
   oplot,men.field1(0,selectvalues),men.field1(1,selectvalues)*scale,symsize=symbolsize,psym=8,color=colorstouse(0)   
   error_bars,men.field1(0,selectvalues),men.field1(1,selectvalues)*scale+men.field1(2,selectvalues)*scale,$ 
                                        men.field1(1,selectvalues)*scale-men.field1(2,selectvalues)*scale,color=colorstouse(0),thick=thick   
   ; mendez stars
   selectvalues = where(men.field1(2,*) lt 0)
   ;plotsym,3,/fill;
   oplot,men.field1(0,selectvalues),men.field1(1,selectvalues)*scale,symsize=symbolsize,psym=7,color=colorstouse(0)   
END
