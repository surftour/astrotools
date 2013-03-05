function fload_mergertime, frun

; this needs the centers.txt file to be in the
; output directory

  read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=frun+"/centers.txt"

  if n_elements(cp_time) lt 2 then begin
	read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=frun+"/centers_bh.txt"
  endif

  if n_elements(cp_time) lt 2 then begin
	read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=frun+"/centers_den.txt"
  endif

  mergertime= 0

  for i=0,n_elements(cp_time)-1 do begin
	dx= (cp_cen1[0,i]-cp_cen2[0,i])
	dy= (cp_cen1[1,i]-cp_cen2[1,i])
	dz= (cp_cen1[2,i]-cp_cen2[2,i])
	rdiff= sqrt(dx*dx + dy*dy + dz*dz)

	if ((rdiff lt 0.1) and (mergertime le 0)) then mergertime= cp_time[i]

  endfor


  return, mergertime

end


