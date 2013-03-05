pro fload_blackhole_details, frun, $
		ids= ids, $
		time= time, $
		bhmass= bhmass, $
		mdot= mdot, $
		local_rho= local_rho, $
		local_c= local_c


if not keyword_set(frun) then begin
	print, " "
	print, "fload_blackhole_details, frun"
	print, " "
	return
endif


spawn, "/bin/ls "+frun+"/blackhole_details_* | wc ", result
n_details= long(result[0])
;print, n_details

if n_details le 0 then begin
	time= -1
	return
endif

found_a_detail= 0

for i=0,n_details-1 do begin

	detailnum= strcompress(string(i),/remove_all)
	bhfile= frun+'/blackhole_details_'+detailnum+'.txt'

	spawn, "wc "+bhfile,result
	lines= long(result)

	if lines gt 0 then begin
	  lid= fltarr(lines)
	  ltime= fltarr(lines)
	  lbhmass= fltarr(lines)
	  lmdot= fltarr(lines)
	  lrho= fltarr(lines)
	  lc= fltarr(lines)

	  spawn, "grep ThisTask "+frun+"/blackhole_details_"+detailnum+".txt | wc", result
	  badlines= long(result[0])

	  openr, 1, bhfile
	  print, "opening: ", bhfile, "   lines= ", lines, " (",badlines,")"
	  junk=''
	  for ii=0L,(lines(0)-1) do begin
		readf, 1, junk
		tempjunk= strsplit(junk,/extract,count=count)
		;if count ne 6 then begin
		;	lid(ii)= -1
		;	ltime(ii)= 0.0
		;	lbhmass(ii)= 0.0
		;	lmdot(ii)= 0.0
		;	lrho(ii)= 0.0
		;	lc(ii)= 0.0
		;endif else begin
		if strmid(tempjunk(0),0,8) ne "ThisTask" and count eq 14 then begin       ; this is for Laura's run
			idplus= strcompress(tempjunk(0))             ; BH=id
			lid(ii)= float(strmid(idplus,3,strlen(idplus)-1))
			ltime(ii)= float(tempjunk(1))
			lbhmass(ii)= float(tempjunk(2))
			lmdot(ii)= float(tempjunk(3))
			lrho(ii)= float(tempjunk(4))
			lc(ii)= float(tempjunk(5))
		endif
		;endelse
	  endfor
	  close, 1

	  idx=where(lid ge 0)
	  lid= lid(idx)
	  ltime= ltime(idx)
	  lbhmass= lbhmass(idx)
	  lmdot= lmdot(idx)
	  lrho= lrho(idx)
	  lc= lc(idx)

	  if found_a_detail eq 0 then begin
		ids= lid
		time= ltime
		bhmass= lbhmass
		mdot= lmdot
		local_rho= lrho
		local_c= lc

		found_a_detail= 1
	  endif else begin
		ids= [ids, lid]
		time= [time, ltime]
		bhmass= [bhmass, lbhmass]
		mdot= [mdot, lmdot]
		local_rho= [local_rho, lrho]
		local_c= [local_c, lc]
	  endelse

	endif

endfor



end






