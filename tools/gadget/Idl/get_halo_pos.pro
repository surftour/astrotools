; general syntax:
; get_halo_pos,name,id=id,x0=x0,y0=y0,z0=z0,m200=m200,r200=r200,hid=hid,sigv=sigv,c=c,redshift=redshift
;     name     : filename
;     id       : Halo id (default='a')
;     x0,y0,z0 : position in kpc/h
;     m200,r200: Vir. properties in Msol/h and Mpc/h
;     sigv     : velocity dispersion
;     c        : nfw conzentration (if nfw fit is present)
;     redshift : redshift of the output
;
; Example:
; get_halo_pos,'snap_092',x0=x0,y0=y0,z0=z0
; print,x0,y0,z0

pro get_halo_pos,name,id=id,x0=x0,y0=y0,z0=z0,m200=m200,r200=r200,hid=hid,sigv=sigv,c=c,redshift=redshift
   IF NOT(IS_DEF(id)) THEN id='a'

   myname,name,fullpath=fullpath,snapnr=snapnr,/quiet,model=model

   readmatrix,aa,fullpath+'/Post/main_prog.'+id+'_dm.gv',/quiet

   jj=WHERE(aa(*,0) EQ double(snapnr))

   IF jj(0) GT -1 THEN BEGIN
      m200=aa(jj(0),4)
      r200=aa(jj(0),6)
      x0=aa(jj(0),9)*1000
      y0=aa(jj(0),10)*1000
      z0=aa(jj(0),11)*1000
      hid=aa(jj(0),1)
   END ELSE BEGIN
      print,'WARNING: No halo found for given file'
      m200=-1
      r200=-1
      x0=-1
      y0=-1
      z0=-1
      hid=-1
      sigv=-1
      c=-1
   END
   IF hid GT 0 THEN BEGIN
;      readmatrix,aa,fullpath+'/Post/SO.'+snapnr+'.gv',/quiet,LastLines=-1
      readmatrix,aa,fullpath+'/Post/SO.'+snapnr+'.gv',/quiet,ReadLines=hid
      IF model EQ 'dm' THEN sigv=aa(hid-1,12) else sigv=aa(hid-1,14)
      ssf=SIZE(FINDFILE(fullpath+'/nfw_fit_'+id+'_logr.dat'))
      readnew,name,h,'HEAD',/quiet
      redshift=h.redshift
      IF ssf(0) NE 0 THEN BEGIN
         readmatrix,aa,fullpath+'/nfw_fit_'+id+'_logr.dat',/quiet,o=1
         dst=ABS(aa(*,0)-redshift)
         ii=SORT(dst)
         c=aa(ii(0),3)
      END
   END
end
