
;==================================================================================
pro read_file_colors, frun, time, mags, $
			funky=funky

cfile= frun+'/colors.txt'
if keyword_set(funky) then cfile= frun+'/colors.txt.disk=1.0'

spawn, "wc "+cfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then c_data= fltarr(15,lines)

openr, 1, cfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, c_data
close, 1


time= c_data[0,*]
mags= c_data[1:14,*]

end




;==================================================================================

