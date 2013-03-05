PRO density3d,ppp,vvv,mass,index,xxx,yyy,zzz,$
                    xcenter,ycenter,zcenter,bin,mindist,maxdist,x,y,count

   n          = long(maxdist-mindist)/bin
   x          = fltarr(n)
   y          = fltarr(n)
   count      = fltarr(n)
   distance3D = sqrt((ppp(xxx,index)-xcenter)^2 + (ppp(yyy,index)-ycenter)^2 + (ppp(zzz,index)-zcenter)^2)   
   vel3d      = sqrt(vvv(0,*)^2 + vvv(1,*)^2 + vvv(2,*)^2)      
   ind        = long((distance3D(*)-mindist)/bin)

   for i=0L,n-1 do begin
     iii      = where(ind eq i)
     if (min(iii) ge 0) then begin
        x(i)     = mean(distance3D(iii))
        y(i)     = total(mass(iii)) / (4.*!Pi/3.) / ((mindist+(i+1.0)*bin)^3-(mindist+i*bin)^3) / bin^3
        count(i) = n_elements(iii)
      endif	
   endfor
   
   cleanup,x,y,count
   
END
