; --------------------------------
;  Read R_e file
; ----------------------------------
;pro read_file_R_e, frun, time, re, reerr, re2, re2err, gre, greerr
pro read_file_r_e, frun, snapext, time, re, reerr

;refile= frun+'/R_e.txt'
refile= frun

spawn, "wc "+refile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then re_data= fltarr(4,lines)

openr, 1, refile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, re_data
close, 1


snapext= re_data[0,*]
time= re_data[1,*]
re= re_data[2,*]
reerr= re_data[3,*]

;time= re_data[0,*]
;re= re_data[1,*]
;reerr= re_data[2,*]
;re2= re_data[3,*]
;re2err= re_data[4,*]
;gre= re_data[5,*]
;greerr= re_data[6,*]


end

