pro make_uvfile,fname,lamuv
	structin = mrdfits(fname,1)
	dims = size(structin.LAMBDA)
	nmet = dims[1]
	nage = dims[2]
	nlam = dims[3]
	uvlum = dblarr(nmet,nage)
	FOR ii=0,nmet-1 DO BEGIN
		FOR jj=0,nage-1 DO BEGIN
			tabinv,alog(structin.LAMBDA[ii,jj,*]),[alog(lamuv)],ix
			uvlum[ii,jj] = exp(interpolate(alog(structin.LAMBDA_LLAMBDA[ii,jj,*]),ix))/lamuv
		ENDFOR
	ENDFOR
	outstruct = {BC03IN:fname,LAMBDA:lamuv,FLAMBDA_UV:uvlum,$
			AGE:structin.AGE,METALLICITY:structin.METALLICITY}
	mwrfits,outstruct,'uvtable.fits'
END

function jy_get_uv,uvfn,stars_age,stars_mets,stars_mass
	lsol = double(1.99e33)
	structin = mrdfits(uvfn,1)
	help,structin,/str
	tabinv,alog(structin.metallicity),alog(stars_mets),im
	tabinv,alog(structin.age),alog(stars_age*1.0e9),ja
	;uvout = bilinear(alog(structin.FLAMBDA_UV),im,ja)
        uvout = interpolate(alog(structin.FLAMBDA_UV[4,*]),ja)
	return,exp(uvout)*stars_mass*1.0e10/lsol
END
