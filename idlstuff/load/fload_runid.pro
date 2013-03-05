function fload_runid, dummy


  COMMON FileInfo

  fid= 'nothing yet'

  runid=strmid(rundir,strpos(rundir,'/',/reverse_search)+1)

  return, runid

end


