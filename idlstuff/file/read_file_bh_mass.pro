; --------------------------------
;  Read bh_mass file
; ----------------------------------
pro read_file_bh_mass, frun, time, bhmass

bhfile= frun+'/bh_mass.txt'

spawn, "wc "+bhfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then bh_data= fltarr(2,lines)

openr, 1, bhfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, bh_data
close, 1


time= bh_data[0,*]
bhmass= bh_data[1,*]

end


