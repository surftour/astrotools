; --------------------------------
;  Read sigma file
; ----------------------------------
pro read_file_sigma, frun, snapext, time, sigxy, sigxz, sigyz, sigavg, sigerr

;if not keyword_set(sigfile) then sigfile= frun+'/sigma.txt' else sigfile=frun+'/'+sigfile
sigfile= frun

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


snapext= sig_data[0,*]
time= sig_data[1,*]
sigxy= sig_data[2,*]
sigxz= sig_data[3,*]
sigyz= sig_data[4,*]
sigavg= sig_data[5,*]
sigerr= sig_data[6,*]


end



