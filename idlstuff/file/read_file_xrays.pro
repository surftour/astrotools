
; ----------------------------
;  Read xrays.txt file
; ----------------------------
pro read_file_xrays, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, $
				xray, xray_sf, $
				xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

hgasfile= frun+'/xrays.txt'

spawn, "wc "+hgasfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then hgas_data= fltarr(11,lines)

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

; with new xrays.txt file (i.e. it has X-ray emission weighted Z)
z_X= hgas_data[2,*]
mass_xraygas= hgas_data[3,*]
entropy= hgas_data[4,*]
xray= hgas_data[5,*]
xray_sf= hgas_data[6,*]
xray_rs_s= hgas_data[7,*]
xray_rs_h= hgas_data[8,*]
xray_rs0_s= hgas_data[9,*]
xray_rs0_h= hgas_data[10,*]

end





; -------------------------------------------------------------------------------------------


