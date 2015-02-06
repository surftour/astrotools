pro myname,name,fullname=fullname,fullpath=fullpath,snapnr=snapnr,$
                model=model,cluste=cluster,quiet=quiet
   CD,'.',CURRENT=old_dir

   stripart=STRMID(name,0,1)
   IF stripart EQ '/' then fullname=name ELSE fullname=old_dir+'/'+name
   IF NOT(IS_DEF(quiet)) THEN print,'Analyzing filename: <',fullname,'>'

   len=STRLEN(fullname)

   IF (len GT 3) THEN snapnr=STRMID(fullname,len-3,3) else snapnr='xxx'
   IF NOT(IS_DEF(quiet)) THEN print,'Snapshot NR: ',snapnr

   pos1=STRPOS(fullname,'/',/REVERSE_SEARCH)
   IF (pos1-1 GT 0) THEN fullpath=STRMID(fullname,0,pos1)
   IF NOT(IS_DEF(quiet)) THEN print,'FullPath: <',fullpath,'>'

   IF (pos1-1 GT 0) THEN pos2=STRPOS(fullname,'/',pos1-1,/REVERSE_SEARCH)
   IF ((pos1 GT 0) AND (pos2 GT 0)) THEN $
      model=STRMID(fullname,pos2+1,pos1-pos2-1) $
   ELSE $
      model='xxx'
   IF NOT(IS_DEF(quiet)) THEN print,'Name of model: <',model,'>'

   IF (pos2-1 GT 0) THEN pos3=STRPOS(fullname,'/',pos2-1,/REVERSE_SEARCH) ELSE pos3=-1
   IF ((pos2-1 GT 0) AND (pos3 GT 0)) THEN $
      cluster=STRMID(fullname,pos3+1,pos2-pos3-1) $
   ELSE $
      cluster='xxx'
   IF NOT(IS_DEF(quiet)) THEN print,'Name of cluster: <',cluster,'>'


end
