function fload_frun, dummy, runidonly=runidonly


  COMMON FileInfo

  fid= 'nothing yet'

  fulldir = strcompress(string(rundir),/remove_all)
  
  if keyword_set(runidonly) then fulldir=strmid(rundir,strpos(rundir,'/',/reverse_search)+1)

  return, fulldir

end


