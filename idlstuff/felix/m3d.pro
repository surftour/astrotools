pro m3d,ppp,vvv,mass,index,xxx,yyy,zzz,$
                    xcenter,ycenter,zcenter,bin,mindist,maxdist,x,y,count
		    
	
   n          = long(maxdist-mindist)/bin
   x          = fltarr(n)
   y          = fltarr(n)
   count      = fltarr(n)
   distance3D = sqrt((ppp(xxx,index)-xcenter)^2 + (ppp(yyy,index)-ycenter)^2 + (ppp(zzz,index)-zcenter)^2)      
   ind        = long((distance3D(*)-mindist)/bin)

   cummass = 0.0
   for i=0L,n-1 do begin
     iii      = where(ind eq i)
     if (min(iii) ge 0) then begin
        x(i)     = max(distance3D(iii))
        cummass  = cummass + total(mass(iii))
        y(i)     = cummass
        count(i) = n_elements(iii)
     endif	
   endfor	    

   cleanup,x,y,count

end
