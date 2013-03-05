;===============================================================================
pro write_line_to_file, filename, x, $
			col1txt=col1txt, $
			comment=comment


if not file_test(filename) then begin
        openw, 1, filename
endif else begin
        openu, 1, filename, /append
endelse



if keyword_set(comment) then goto, addcomment


;---------------
;
numstring= strcompress(string(n_elements(x)-1),/remove_all)
fmtstring= '(F8.4,"  ",'+numstring+'(F10.4,"  "))'
if keyword_set(col1txt) then fmtstring='(" '+col1txt+' ",F8.4,"  ",'+numstring+'(F10.4,"  "))'

printf, 1, FORMAT=fmtstring, x
close, 1
return


;---------------
;
addcomment:
comment= string(comment)
printf, 1, comment
close, 1
return


end



;===============================================================================




