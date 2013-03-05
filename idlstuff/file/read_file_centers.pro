;===============================================================================
;
;  not currently used
;
;
;pro read_centerpositions, time, center_1, center_2, filename=filename
;        
;        
;        
;; read this with the following
;openr, 1, filename, ERROR=err
;junk= ''
;lines= 0
;
;if err NE 0 then begin
;	print, " "
;	print, " ERROR found when opening file. "
;	print, " file: ", filename
;	print, " err: ", err
;	print, " "
;	close, 1
;	return
;endif
;
; galaxy 1
;readf, 1, junk
;readf, 1, lines
;data= fltarr(4,lines)
;readf, 1, data  
;time= data[0,*] 
;gal1_xyz= data[1:3,*]
;center_1= gal1_xyz
;                
;; galaxy 2
;readf, 1, junk
;readf, 1, lines
;data= fltarr(4,lines)
;readf, 1, data
;gal2_xyz= data[1:3,*]
;center_2= gal2_xyz
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
;===============================================================================



pro read_file_centers, time, center, filename=filename
        
        
cmd= "wc "+filename
spawn, cmd, result
wclines= long(result[0])
        
; read this with the following
print, "opening: ", filename
openr, 1, filename, ERROR=err
junk= ''
lines= 0

if err NE 0 then begin
        print, " "
        print, " ERROR found when opening file. "
        print, " file: ", filename
        print, " err: ", err
        print, " "
        close, 1
        return
endif

; galaxy 1
;readf, 1, bhid
;readf, 1, junk
;readf, 1, lines
;if lines lt 0 then lines= wclines-2
lines= wclines
data= fltarr(4,lines)
readf, 1, data  
time= data[0,*]
gal1_xyz= data[1:3,*]
center= gal1_xyz
                
           
close, 1   


end        





