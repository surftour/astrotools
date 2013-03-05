pro readcenters, filename,distance,time

 openr, 1, filename
 junk= ''
 lines= 0
 readf, 1, junk     ; galaxy 1
 readf, 1, lines
 data= fltarr(4,lines)
 readf, 1, data
 time= data[0,*]
 gal1xyz= data[1:3,*]
 readf, 1, junk     ; galaxy 1
 readf, 1, lines
 data= fltarr(4,lines)
 readf, 1, data
 gal2xyz= data[1:3,*]
 
 distance = sqrt( (gal1xyz(0,*)-gal2xyz(0,*))^2 + $
                  (gal1xyz(1,*)-gal2xyz(1,*))^2 + $ 
                  (gal1xyz(2,*)-gal2xyz(2,*))^2 )
 
 close,1
 
 ;window,3
 ;plot,gal1xyz(0,*),gal1xyz(1,*)
 ;oplot,gal2xyz(0,*),gal2xyz(1,*),color=250
 ;window,0
end
