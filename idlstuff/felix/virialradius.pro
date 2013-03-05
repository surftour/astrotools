PRO virialradius,ppp,mmm,xxx,yyy,zzz,xcenter,ycenter,zcenter,radius

   ; we do as if all the haloes were at z=0
   ; this is the value of virial density (w/o 4pi/3) in units of Msun kpc which are the units of the simulations
   densityvalue = 323.0 * 0.3 * 2.775e11 * 0.7^2 * 1e-9 
   
   distance3D   = sqrt((ppp(xxx,*)-xcenter)^2 + (ppp(yyy,*)-ycenter)^2 + (ppp(zzz,*)-zcenter)^2)
   mass3D       = double(mmm * (1.0e10 / (4.0*!Pi/3.0))); changing into units of Msun and including the density factor
   result       = sort(distance3D)
   
   distance3D   = distance3D(result)
   mass3D       = mass3D(result)
   
   ; convert mass to cumulative mass
   FOR i=1L,n_elements(distance3D)-1L DO BEGIN
      mass3D(i) = mass3D(i) + mass3D(i-1L)
   ENDFOR
   
   ; get the radius of 200 times the mean density
   i = 10L
   WHILE (mass3D(i)/distance3D(i)/distance3D(i)/distance3D(i) gt densityvalue) DO BEGIN
      ;;print,i,mass3D(i)/distance3D(i)/distance3D(i)/distance3D(i),densityvalue
      i = i + 1
   ENDWHILE
   
   radius = distance3D(i-1)
   print,'Virialradius',radius

END
