;==================================
;  Load the ID list
;==================================
function fload_id_list_suv, idfilename


spawn, 'wc '+idfilename, result
lines= long(result)
data= fltarr(lines-5,4)
idlist= lonarr(lines-5)

openr, 2, idfilename

textjunk= ''
readf, 2, textjunk
readf, 2, textjunk
readf, 2, textjunk
readf, 2, textjunk
readf, 2, textjunk
readf, 2, data
close, 2

idlist(*)= long(data(*,0))

return, idlist


end


