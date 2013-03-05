
;========================================
;========================================
;
;
; procedure reads phil's list of 
; directories, and his processed data
;
;
function load_list_from_phil_file, philfile, isstring=isstring


phildir= '/n/home/tcox/Documents/Data_PhilsSimSersicFits/'

; interesting files in here:
;
;       burstmassfraction
;       n_s_pure
;       r_eff_projected_all
;       snapdirs (make sure to set /isstring)
;
;

phil_file_todo= phildir+philfile

spawn, "wc "+phil_file_todo,result
lines=long(result)
datalines=lines(0)-5

openr, 1, phil_file_todo
junk=''

phillist= fltarr(datalines)
if keyword_set(isstring) then phillist= strarr(datalines)

; read the data
for i=0,datalines-1 do begin

        readf, 1, junk
	;print, junk

	; if it's multiple items per line
        tempjunk= strsplit(junk,/extract,count=count)
	;tempjunk= strcompress(junk,/remove_all)

	if tempjunk(0) ne "#" then begin
        	if keyword_set(isstring) then phillist(i)= tempjunk(0) else phillist(i)= float(tempjunk(0))
	endif else begin
		i=i-1
	endelse

endfor

close, 1


return, phillist


end




;========================================
;========================================

pro doit, junk

	print, "--------------------------------------------------------------"
	print, " "
	print, "opening: snapdirs.txt"
	fruns= load_list_from_phil_file("snapdirs.txt", /isstring)
	print, " "
	print, "opening: r_eff_projected_all.txt"
	r_eff= load_list_from_phil_file("r_eff_projected_all.txt")
	print, " "
	print, "opening: n_s_pure.txt"
	n_s= load_list_from_phil_file("n_s_pure.txt")
	n_s_m= load_list_from_phil_file("n_s_pure_m.txt")
	n_s_p= load_list_from_phil_file("n_s_pure_p.txt")
	print, " "
	print, "opening: burstmassfraction.txt"
	fb= load_list_from_phil_file("burstmassfraction.txt", /isstring)
	print, " "
	print, "--------------------------------------------------------------"

stop
end




;========================================
;========================================




