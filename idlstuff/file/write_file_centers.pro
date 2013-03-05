;===============================================================================
;
;  not currently used
;
;
;pro write_centerpositions, time, center_1, center_2, filename=filename
;
;
;
;openw, 1, filename
;
;; galaxy 1
;printf, 1, 'galaxy 1'
;
;lines= n_elements(time)
;printf, 1, lines
;
;data= fltarr(4,lines)
;data[0,*]= time
;c1= transpose(center_1)
;data[1:3,*]= c1
;printf, 1, data
;
;; galaxy 2
;printf, 1, 'galaxy 2'
;printf, 1, lines
;c2= transpose(center_2)
;data[1:3,*]= c2
;printf, 1, data
;
;
;close, 1
;
;
;end
;
;
;
;
;
;===============================================================================
;
;
;
;pro write_center, time, center, filename=filename, bhid=bhid
;
;
;
;openw, 1, filename
;
;; galaxy 1
;printf, 1, bhid
;
;lines= n_elements(time)
;printf, 1, lines
;
;data= fltarr(4,lines)
;data[0,*]= time
;c1= transpose(center)
;data[1:3,*]= c1
;printf, 1, data
;
;close, 1 
;
;
;end
;











