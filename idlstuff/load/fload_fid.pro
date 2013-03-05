function fload_fid, dummy


  COMMON FileInfo

  ;fid= 'nothing yet'
  fid= ' '

  fulldir = strcompress(string(rundir),/remove_all)

  ; if it's a /data/tcox/ directory
  if strmid(fulldir,1,4) eq 'data' then begin
    tpos=strpos(fulldir,"tcox")
    if tpos LT 0 then fpos= 0 else fpos=tpos+5
    fid=strmid(fulldir,fpos)
  endif

  ; if it's a data/ directory
  if strmid(fulldir,0,5) eq 'data/' then fid=strmid(fulldir,5)

  ; if it's a data1/ directory
  if strmid(fulldir,0,6) eq 'data1/' then fid=strmid(fulldir,6)

  ; if it's /home/tcox/Data/ directory
  if strmid(fulldir,1,4) eq 'home' then fid=strmid(fulldir,16)

  ; if it's just Data/ directory
  if strmid(fulldir,0,4) eq 'Data' then fid=strmid(fulldir,5)

  ; if it's pool/ directory
  if strmid(fulldir,0,4) eq 'pool' then fid=strmid(fulldir,5)

  ; if it's execute
  if strmid(fulldir,0,7) EQ 'execute' then fid=strmid(fulldir,8)

  ; if it's /raid#/tcox/
  if strmid(fulldir,1,4) EQ 'raid' then fid=strmid(fulldir,12)

  ; if it's /odyssey/tcox/
  if strmid(fulldir,1,7) EQ 'odyssey' then fid=strmid(fulldir,14)


  ; just return the first 2 characters
  if string(dummy) eq 'firsttwo' then fid=strmid(fid,0,2)


  return, fid

end


