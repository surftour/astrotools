function tdwarfs_read_fofgroup,fdata_base,ngroups,nsubgroups,npart,fofcat,fofprop,subcat,subprop,ids,mass,pos,vel,types,rho,u,nel,nHI,hsml,sfr,age,z

	;fof cat
	fdata = fdata_base+'.fofcat'
	openr,1,fdata
	ngroups = 0L
	readu,1,ngroups
	print,'ngroups = ',ngroups

	fofcat = {length:lonarr(ngroups),offset:lonarr(ngroups)}
	length = lonarr(ngroups)
	offset = lonarr(ngroups)
	readu,1,length
	readu,1,offset
	fofcat.length(*) = length
	fofcat.offset(*) = offset
	length = 0L
	offset = 0L
	close,1
	;print,fofcat.length(0),fofcat.offset(0)

	;sub cat
	fdata = fdata_base+'.subcat'
	openr,1,fdata
	nsubgroups = 0L
	readu,1,nsubgroups
	print,'nsubgroups = ',nsubgroups

	if nsubgroups gt 0 then begin
	subcat = {length:lonarr(nsubgroups),offset:lonarr(nsubgroups)}
	length = lonarr(nsubgroups)
	offset = lonarr(nsubgroups)
	readu,1,length
	readu,1,offset
	subcat.length(*) = length
	subcat.offset(*) = offset
	length = 0L
	offset = 0L
	endif 
	close,1
	;print,subcat.length(0),fofcat.offset(0)

	;fof props
	fdata = fdata_base+'.fofprop'
	openr,1,fdata
	ngb = 0L
	readu,1,ngb
	;print,'ngroups.prop = ',ngb
	if(ngb ne ngroups) then begin
		print,'Catalogue lengths do not agree!'
		close,1
	       	return,1
	endif
	fofprop = {cm:fltarr(ngroups,3),m_tot:fltarr(ngroups),m:fltarr(5,ngroups),sfr:fltarr(ngroups)}
	cm = fltarr(ngroups,3)
	readu,1,cm
	fofprop.cm = cm
	cm = 0.0
	m = fltarr(ngroups)
	readu,1,m
	fofprop.m_tot= m*10.0^(10.0)
	;for i=0,5 do begin
	;	readu,1,m
	;	fofprop.m(i,*) = m
	;	print,i
	;endfor
	;m = 0.0
	;sfr = fltarr(ngroups)
	;readu,1,sfr
	;fofprop.sfr = sfr
	;sfr =  0.0
	close,1

	;sub props
	fdata = fdata_base+'.subprop'
	openr,1,fdata
	ngb = 0L
	readu,1,ngb
	;print,'ngroups.prop = ',ngb
	if ngb gt 0 then begin
	if(ngb ne nsubgroups) then begin
		print,'Subhalo catalogue lengths do not agree!'
		close,1
	       	return,1
	endif
	subprop = {cm:fltarr(nsubgroups,3),m_tot:fltarr(nsubgroups),m:fltarr(5,nsubgroups),sfr:fltarr(nsubgroups)}
	cm = fltarr(nsubgroups,3)
	readu,1,cm
	subprop.cm = cm
	cm = 0.0
	m = fltarr(nsubgroups)
	readu,1,m
	subprop.m_tot= m*10.0^(10.0)
	;for i=0,5 do begin
	;	readu,1,m
	;	subprop.m(i,*) = m
	;	print,i
	;endfor
	;m = 0.0
	;sfr = fltarr(ngroups)
	;readu,1,sfr
	;subprop.sfr = sfr
	;sfr =  0.0
	endif
	close,1

	;ids
	fdata = fdata_base+'.ids'
	openr,1,fdata
	npart = 0L
	readu,1,npart
	print,'npart = ',npart
	ids = lonarr(npart)
	readu,1,ids
	close,1

	;mass
	fdata = fdata_base+'.mass'
	openr,1,fdata
	npb = 0L
	readu,1,npb
	;print,'npart.mass = ',npb
	if(npb ne npart) then begin
		close,1
		print,'Mass npart ne ids npart!'
		return,1
	endif
	mass = fltarr(npart)
	readu,1,mass
	close,1
	mass(*) = mass(*)*10.0^(10.0)

	;pos
	fdata = fdata_base+'.pos'
	openr,1,fdata
	npb = 0L
	readu,1,npb
	if(npb ne npart) then begin
		close,1
		print,'Pos npart ne ids npart!'
		return,1
	endif
	pos = fltarr(3,npart)
	readu,1,pos
	close,1

	;vel
	fdata = fdata_base+'.vel'
	openr,1,fdata
	npb = 0L
	readu,1,npb
	if(npb ne npart) then begin
		close,1
		print,'Vel npart ne ids npart!'
		return,1
	endif
	vel = fltarr(3,npart)
	readu,1,vel
	close,1

	;types
	fdata = fdata_base+'.types'
	openr,1,fdata
	npb = 0L
	readu,1,npb
	if(npb ne npart) then begin
		close,1
		print,'Types npart ne ids npart!'
		return,1
	endif
	types = lonarr(npart)
	readu,1,types
	close,1

	;u
	fdata = fdata_base+'.u'
	openr,1,fdata
	npb = 0L
	readu,1,npb
	if(npb ne npart) then begin
		close,1
		print,'u npart ne ids npart!'
		return,1
	endif
	u = fltarr(npart)
	readu,1,u
	close,1

	;rho
	fdata = fdata_base+'.rho'
	openr,1,fdata
	npb = 0L
	readu,1,npb
	if(npb ne npart) then begin
		close,1
		print,'rho npart ne ids npart!'
		return,1
	endif
	rho= fltarr(npart)
	readu,1,rho
	close,1

	;nHI
	fdata = fdata_base+'.nh'
	openr,1,fdata
	npb = 0L
	readu,1,npb
	if(npb ne npart) then begin
		close,1
		print,'NH npart ne ids npart!'
		return,1
	endif
	nHI = fltarr(npart)
	readu,1,nHI
	close,1

	;nel
	fdata = fdata_base+'.nel'
	openr,1,fdata
	npb = 0L
	readu,1,npb
	if(npb ne npart) then begin
		close,1
		print,'nel npart ne ids npart!'
		return,1
	endif
	nel= fltarr(npart)
	readu,1,nel
	close,1

	;hsml
	fdata = fdata_base+'.hsml'
	openr,1,fdata
	npb = 0L
	readu,1,npb
	if(npb ne npart) then begin
		close,1
		print,'hsml npart ne ids npart!'
		return,1
	endif
	hsml= fltarr(npart)
	readu,1,hsml
	close,1


	;sfr
	fdata = fdata_base+'.sfr'
	openr,1,fdata
	npb = 0L
	readu,1,npb
	if(npb ne npart) then begin
		close,1
		print,'sfr npart ne ids npart!'
		return,1
	endif
	sfr= fltarr(npart)
	readu,1,sfr
	close,1


	;age
	fdata = fdata_base+'.age'
	openr,1,fdata
	npb = 0L
	readu,1,npb
	if(npb ne npart) then begin
		close,1
		print,'age npart ne ids npart!'
		return,1
	endif
	sfr= fltarr(npart)
	readu,1,age
	close,1

	;z
	fdata = fdata_base+'.zm'
	openr,1,fdata
	npb = 0L
	readu,1,npb
	if(npb ne npart) then begin
		close,1
		print,'Z npart ne ids npart!'
		return,1
	endif
	z= fltarr(npart)
	readu,1,z
	close,1


	return,0

end
