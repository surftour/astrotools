function fload_getid, frun



  if strlen(frun) GT 0 then begin

    fulldir = strcompress(string(frun),/remove_all)
    tpos=strpos(fulldir,"tjcox")
    if tpos LT 0 then begin
	fpos= 0 
    endif else begin
 	fpos=tpos+5
    endelse
    fid=strmid(fulldir,fpos)

    if strmid(fulldir,0,4) eq 'Data' then fid=strmid(fulldir,5)

    if strmid(fulldir,0,4) eq 'pool' then fid=strmid(fulldir,5)

    return, fid

  endif

end


