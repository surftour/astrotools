; --------------------------------
;  Read desika FWHM file
; ----------------------------------
pro read_file_desika_FWHM, dfile, time, FWHMxy, FWHMxz, FWHMyz


spawn, "wc "+dfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then d_data= fltarr(4,lines)

openr, 1, dfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, d_data
close, 1


time= d_data[0,*]
FWHMxy= d_data[1,*]
FWHMxz= d_data[2,*]
FWHMyz= d_data[3,*]


end




