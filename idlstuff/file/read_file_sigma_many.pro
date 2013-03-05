; --------------------------------
;  Read many sigma file
; ----------------------------------
pro read_many_sigma_file, sigfile, time, sig1, sig2, sig3, sig4, sig5, sig6

;sigfile= frun+'/sigma.txt'

spawn, "wc "+sigfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then sig_data= fltarr(7,lines)

openr, 1, sigfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, sig_data
close, 1


time= sig_data[0,*]
sig1= sig_data[1,*]
sig2= sig_data[2,*]
sig3= sig_data[3,*]
sig4= sig_data[4,*]
sig5= sig_data[5,*]
sig6= sig_data[6,*]


end

