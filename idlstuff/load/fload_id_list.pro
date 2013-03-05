;==================================
;  Load the ID list
;==================================
function fload_id_list, idfilename


spawn, 'wc '+idfilename, result
lines= long(result)
idlist= lonarr(lines-5)

openr, 2, idfilename

textjunk= ''
readf, 2, textjunk
readf, 2, textjunk
readf, 2, textjunk
readf, 2, textjunk
readf, 2, textjunk
readf, 2, idlist
close, 2


return, idlist


end


