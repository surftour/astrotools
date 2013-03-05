


; ----------------------------
;  Read hotgas.txt file
; ----------------------------
pro read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, $
				gas_tot, gas_hot, gas_cold, gas_sf

hgasfile= frun+'/hotgas.txt'

spawn, "wc "+hgasfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then hgas_data= fltarr(9,lines)

openr, 1, hgasfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, hgas_data
close, 1


time= hgas_data[0,*]
temp_keV_X= hgas_data[1,*]
temp_keV= hgas_data[2,*]
temp_K= hgas_data[3,*]
entropy= hgas_data[4,*]
gas_tot= hgas_data[5,*]
gas_hot= hgas_data[6,*]
gas_cold= hgas_data[7,*]
gas_sf= hgas_data[8,*]


end



