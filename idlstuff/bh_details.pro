

;====================================================================




pro save_bhdetails_at_snaptimes, frun


if not keyword_set(frun) then begin
	print, " "
	print, " save_bhdetails_at_snaptimes, frun "
	print, " "
	return
endif



;
; constants
;
GRAVITY= 6.672e-8
CC= 2.9979d+10
PROTONMASS= 1.6726e-24
THOMPSON= 6.6524e-25
UnitTime_in_s= 3.08568e+16



; this assumes they are ordered
;  0 through


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snapshot_* | wc ",result
nsnaps=long(result[0])

if not keyword_set(startsnap) then startsnap= 0
if not keyword_set(endsnap) then endsnap= nsnaps


time= fltarr(nsnaps)




; what happens if we don't assume it's a major
; merger, but instead we check for BH's and do it from
; there

ok= fload_snapshot_bh(frun,0,/nopot_in_snap)
init_bhids= fload_blackhole_id(1)
init_bhids= long(init_bhids(sort(init_bhids)))
n_bh= n_elements(init_bhids)
print, "Blackhole ID's: ", init_bhids


;
; write files
;
for bhi= 0, n_bh-1 do begin

        bhlbl= strcompress(string(init_bhids[bhi]),/remove_all)
	bhfilename= frun+"/bh_"+bhlbl+".txt"
	print, "opening: ", bhfilename
	openw, 1, bhfilename
	printf, 1, "# frun= "+frun
	printf, 1, "# bhid= "+bhlbl
	printf, 1, "#  "
	printf, 1, "#    Time       BH Mass   Part. Mass        M_dot        M_dot  "
	printf, 1, "#  (Gyr/h)     (Gad. U)     (Gad. U)    (M_sun/Yr)       (Edd.) "
	close, 1
endfor




; get bolometric luminosity from the accretion rate
; --------------------------------------------------
fload_blackhole_details, frun, ids=bhdids, time=bhdtime, bhmass=bhdmass, mdot=bhdmdot, local_rho=lr, local_c=lc


sorttime= sort(bhdtime)
sorted_time= bhdtime(sorttime)
sorted_bhid= bhdids(sorttime)
sorted_bhm= bhdmass(sorttime)
sorted_mdot= bhdmdot(sorttime)





; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------


for i=startsnap,endsnap-1 do begin

        print, "--------------------------------------"


        ; open snapshot
        if i ne 0 then begin
          if i ge 1000 then begin
                 ok=fload_snapshot_bh(frun,i,/skip_center,/nopot_in_snap,/do_four)
          endif else begin
                 ok=fload_snapshot_bh(frun,i,/skip_center,/nopot_in_snap)
          endelse
        endif

        ; what time is it?
        time= fload_time(1)


        ;
        ; ---------------------------------------------------
        for bhi= 0, n_bh-1 do begin

		bhid= init_bhids[bhi]

                print, "bhi= ", bhi
		print, "bhid= ", bhid

		idx= where((sorted_bhid eq bhid) and (sorted_time ge time))
		if idx(0) ne -1 then begin

			idxrange= 3
			if idx(0) le idxrange then lowidx= 0 else lowidx= idx(0)-idxrange
			if idx(0) ge n_elements(sorted_time)-idxrange then hiidx= n_elements(sorted_time)-1 else hiidx= idx(0)+idxrange

			pm= fload_blackhole_mass(2,idtofollow=bhid)

			; instantaneous values
			;bhm= sorted_bhm(idx(0))
			;mdot= sorted_mdot(idx(0))
			;meddington = 4.0 * !PI * GRAVITY * CC * PROTONMASS / (0.1 * CC * CC * THOMPSON) * bhm * UnitTime_in_s
			;mdotedd= mdot/meddington

			; time averaged
			bhm= sorted_bhm(lowidx:hiidx)
			mdot= sorted_mdot(lowidx:hiidx)
			ids= sorted_bhid(lowidx:hiidx)
			meddington = 4.0 * !PI * GRAVITY * CC * PROTONMASS / (0.1 * CC * CC * THOMPSON) * bhm * UnitTime_in_s
			mdotedd= mdot/meddington
			ididx= where(ids eq bhid)
			bhm= mean(bhm(ididx))
			mdot= mean(mdot(ididx))
			mdotedd= mean(mdotedd(ididx))
		endif else begin
			bhm= 0.0
			pm= 0.0	
			mdot= 0.0
			mdotedd= 0.0
		endelse

		; add info to bh file
		bhlbl= strcompress(string(init_bhids[bhi]),/remove_all)
		bhfilename= frun+"/bh_"+bhlbl+".txt"
		openu, 1, bhfilename, /append
		printf, 1, FORMAT= '("   ", F6.3, "    ",4(F11.8,"  "))', time, bhm, pm, mdot, mdotedd
		close, 1


        endfor

endfor


end







;======================================================




pro compile_bhdetails, frun, dt


if not keyword_set(frun) then begin
	print, " "
	print, " compile_bhdetails, frun, dt "
	print, " "
	return
endif


if not keyword_set(dt) then dt= 0.01   ; display every 10 million years.

;
; constants
;
GRAVITY= 6.672e-8
CC= 2.9979d+10
PROTONMASS= 1.6726e-24
THOMPSON= 6.6524e-25
UnitTime_in_s= 3.08568e+16





; what happens if we don't assume it's a major
; merger, but instead we check for BH's and do it from
; there

ok= fload_snapshot_bh(frun,0,/nopot_in_snap)
init_bhids= fload_blackhole_id(1)
init_bhids= long(init_bhids(sort(init_bhids)))
n_bh= n_elements(init_bhids)
print, "Blackhole ID's: ", init_bhids


;
; write files
;
for bhi= 0, n_bh-1 do begin

        bhlbl= strcompress(string(init_bhids[bhi]),/remove_all)
	bhfilename= frun+"/bh_"+bhlbl+".txt"
	print, "opening: ", bhfilename
	openw, 1, bhfilename
	printf, 1, "# frun= "+frun
	printf, 1, "# bhid= "+bhlbl
	printf, 1, "#  "
	printf, 1, "#    Time       BH Mass   Part. Mass        M_dot        M_dot  "
	printf, 1, "#  (Gyr/h)     (Gad. U)     (Gad. U)    (M_sun/Yr)       (Edd.) "
	close, 1
endfor




; get bolometric luminosity from the accretion rate
; --------------------------------------------------
fload_blackhole_details, frun, ids=bhdids, time=bhdtime, bhmass=bhdmass, mdot=bhdmdot, local_rho=lr, local_c=lc


sorttime= sort(bhdtime)
sorted_time= bhdtime(sorttime)
sorted_bhid= bhdids(sorttime)
sorted_bhm= bhdmass(sorttime)
sorted_mdot= bhdmdot(sorttime)


bhdetails_maxtime= max(bhdtime)



; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

mastertime= 0.0

repeat begin

        print, "--------------------------------------"
	print, " "
        print, "mastertime= ", mastertime
	print, " "


        ;
        ; ---------------------------------------------------
        for bhi= 0, n_bh-1 do begin

		bhid= init_bhids[bhi]

                print, "bhi= ", bhi
		print, "bhid= ", bhid

		;idx= where(sorted_time gt mastertime)
		;if idx(0) eq -1 then break

		idx= where((sorted_bhid eq bhid) and (sorted_time ge mastertime))
		if idx(0) ne -1 then begin

			idxrange= 3
			if idx(0) le idxrange then lowidx= 0 else lowidx= idx(0)-idxrange
			if idx(0) ge n_elements(sorted_time)-idxrange then hiidx= n_elements(sorted_time)-1 else hiidx= idx(0)+idxrange

			;pm= fload_blackhole_mass(2,idtofollow=bhid)
			pm= 0.0

			; instantaneous values
			;bhm= sorted_bhm(idx(0))
			;mdot= sorted_mdot(idx(0))
			;meddington = 4.0 * !PI * GRAVITY * CC * PROTONMASS / (0.1 * CC * CC * THOMPSON) * bhm * UnitTime_in_s
			;mdotedd= mdot/meddington

			; time averaged
			bhm= sorted_bhm(lowidx:hiidx)
			mdot= sorted_mdot(lowidx:hiidx)
			ids= sorted_bhid(lowidx:hiidx)
			meddington = 4.0 * !PI * GRAVITY * CC * PROTONMASS / (0.1 * CC * CC * THOMPSON) * bhm * UnitTime_in_s
			mdotedd= mdot/meddington
			ididx= where(ids eq bhid)
			bhm= mean(bhm(ididx))
			mdot= mean(mdot(ididx))
			mdotedd= mean(mdotedd(ididx))
		endif else begin
			bhm= 0.0
			pm= 0.0	
			mdot= 0.0
			mdotedd= 0.0
		endelse

		; add info to bh file
		bhlbl= strcompress(string(init_bhids[bhi]),/remove_all)
		bhfilename= frun+"/bh_"+bhlbl+".txt"
		openu, 1, bhfilename, /append
		printf, 1, FORMAT= '("   ", F6.3, "    ",4(F11.8,"  "))', mastertime, bhm, pm, mdot, mdotedd
		close, 1


        endfor

	mastertime= mastertime + dt

endrep until (mastertime gt bhdetails_maxtime)


end












;======================================================





; --------------------------------
;  Read bh_mass file
; ----------------------------------
pro read_bh_file, frun, bhid, time, bhm, pm, mdot, mdotedd, $
			 bhfile=bhfile, mergertime=mergertime

bhlbl= strcompress(string(bhid),/remove_all)
if not keyword_set(bhfile) then bhfile= frun+'/bh_'+bhlbl+'.txt'

spawn, "wc "+bhfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then bh_data= fltarr(5,lines)

openr, 1, bhfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, bh_data
close, 1


time= bh_data[0,*]
bhm= bh_data[1,*]
pm= bh_data[2,*]
mdot= bh_data[3,*]
mdotedd= bh_data[4,*]

mergertime= -1
idx= where(bhm le 0.0)
if idx(0) ne -1 then mergertime= time(idx(0))

end



