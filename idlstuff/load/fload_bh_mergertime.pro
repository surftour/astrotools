



;function fload_bh_mergertime, frun

;forward_function get_blackhole_data

;get_blackhole_data, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass, bh_mdot_edd

; get bh merger time
;time_merge= 0
;idx= where(bh_num eq 1)
;if idx(0) ne -1 then time_merge= bhtime(idx(0))

;return, time_merge

;end






;=====================================================


function get_blackhole_data, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass, bh_mdot_edd


; get bh data
;bhfile= '/home/tcox/data/bh/'+fload_getid(frun)+'.bh'
;bhfile= frun+'/'+fload_getid(frun)+'.bh'
spawn, "/bin/ls "+frun+"/*.bh", result
if strlen(result) eq 0 then begin
	print, "-----------------------------------------"
	print, "  Can't file: "+frun+"/*.bh"
	print, "-----------------------------------------"
	return, -1
endif

bhfile=result
bhfile=strcompress(bhfile,/remove_all)


spawn, "wc "+bhfile,result
lines=long(result)
lines=lines(0)
if lines gt 0 then bhdata= fltarr(5,lines)

; open file
openr, 1, bhfile, ERROR=err

if (err NE 0) then begin
        print, "  "
        print, "Problem: ",!ERR_STRING
        print, "  "
	close, 1
	ERR= 0
	bhtime= [0]
	bhmsg= "No Black Hole"
endif else begin
	;close, unit
	print, "opening: ",bhfile
	bhdata= read_ascii(bhfile)
	bhtime= bhdata.field1[0,*]
	bh_num= bhdata.field1[1,*]
	bh_mass= bhdata.field1[2,*]
	bh_mdot_gu= bhdata.field1[3,*]
	bh_mdot_sunyr= bhdata.field1[4,*]
	bh_totalmass= bhdata.field1[5,*]
	bh_mdot_edd= bhdata.field1[6,*]
	;readf,1,bhdata
	;bhtime= bhdata[0,*]
	;sfrsfr= bhdata[1,*]
	close, 1
endelse

return, 0

end




;=====================================================









function fload_bh_mergertime, frun

ok= get_blackhole_data(frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass, bh_mdot_edd)

if ok lt 0 then begin
	return, ok
endif else begin
	; get bh merger time
	time_merge= 0
	idx= where(bh_num eq 1)
	if idx(0) ne -1 then time_merge= bhtime(idx(0))

	return, time_merge
endelse

end










