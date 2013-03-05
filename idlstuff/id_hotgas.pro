pro do_the_following, junk

   ;
   ; generate the ID list, all
   ; the details need to be 
   ; set manually
   ;
   generate_idlist, 1



end






;=========================================
;
;  Do the grunt work of figuring out
; what the id's are of certain particles
;
;=========================================

;
; Write a file that lists the
;  ids of all gas particles which
;  contribute to the X-ray
;  emission.
;-------------------------------
;pro generate_idlist, junk
pro generate_idlist, frun, snapnum

;if not keyword_set(junk) then begin
if not keyword_set(frun) then begin
   print, "  "
   ;print, "generate_idlist, junk"
   print, "generate_idlist, frun, snapnum"
   return
endif



;frun= "/raid4/tcox/vc3vc3e_2"
;snapnum= 107
;frun= "/raid4/tcox/vc3vc3e"
;snapnum= 30
;frun= "/raid4/tcox/vc3vc3e_no"
;snapnum= 107


idlistname= frun+"/hotgas_idlist.txt"


ok= fload_snapshot_bh(frun,snapnum)


gids= fload_gas_id(1)

xray= fload_gas_xray_luminosity(1,/diffuse_hotgas)
xraylum= xray/1.0d+22


; take only xray luminous gas particles
; --------------------------------------
idx= where(xraylum gt 0)
if idx(0) ne -1 then begin
	xr= xraylum(idx)
	gids= gids(idx)
endif
print, 'there are ', n_elements(idx), ' non-zero gas particles'


; now order these
; -----------------
xrindx= sort(xr)
gids= gids(xrindx)
xr= xr(xrindx)


idlist= gids

print, 'writing ', n_elements(idlist), ' ids'


;  Write id list to a file
; ----------------------------
cmt= 'xr emitting gas order by luminosity'
id_filename= idlistname
print, ' writing '+cmt+' to '+id_filename
write_id_list, idlist, id_filename, cmt=cmt



end




; ---------------------------------------------------------------------




;==================================
;  Write the ID list
;==================================
pro write_id_list, idlist, idfilename, cmt=cmt

if not keyword_set(cmt) then cmt=' '

; ----------------------------
;  Write id list to a file
;
get_lun, unit
openw, unit, idfilename

printf, unit, '#  file: '+idfilename
printf, unit, "#  list of gas ids "
printf, unit, "#  "
printf, unit, '#  '+cmt
printf, unit, "#  "
for i=0L,n_elements(idlist)-1 do printf, unit, idlist[i]

close, unit

end





;==================================
;  Load the ID list
;==================================
function load_id_list, idfilename


spawn, 'wc '+idfilename, result
lines= long(result)
idlist= lonarr(lines-5)

get_lun, unit
openr, unit, idfilename

textjunk= ''
readf, unit, textjunk
readf, unit, textjunk
readf, unit, textjunk
readf, unit, textjunk
readf, unit, textjunk
readf, unit, idlist
close, unit


return, idlist


end






; ------------------------------------------------------------------------






;==============================================
;
;  Plots of gas properties as a function of
;  x-ray emission in the merger remnant
;
;==============================================


;
; Generate plot showing the average gas
; metallicity as a function of time for all
; x-ray gas in the merger remnant.
;----------------------------------------------
pro plot_gas_z, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "plot_gas_z, junk"
   return
endif

;  1
; ---
frun= "/raid4/tcox/vc3vc3e_2"
;frun= "/raid4/tcox/vc3vc3e"
;frun= "/raid4/tcox/vc3vc3e_no"

process_metals_time, frun, time=time, z_avg=z_avg, z_err=z_err, $
			z_median= z_median, z_plus25=z_plus25, z_min25=z_min25

time_1=       time
zmets_1=      z_avg
zmetserr_1=   z_err
zmet_25_1=    z_min25
zmet_50_1=    z_median
zmet_75_1=    z_plus25


;  2
; ---
;frun= "/raid4/tcox/vc3vc3e_2"
;frun= "/raid4/tcox/vc3vc3e"
frun= "/raid4/tcox/vc3vc3e_no"

process_metals_time, frun, time=time, z_avg=z_avg, z_err=z_err, $
			z_median= z_median, z_plus25=z_plus25, z_min25=z_min25

time_2=       time
zmets_2=      z_avg
zmetserr_2=   z_err
zmet_25_2=    z_min25
zmet_50_2=    z_median
zmet_75_2=    z_plus25





;----------------
; Print it up
;----------------

;filename=frun+'/newstarage.eps'
filename='zenrich.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

xaxistitle = "Time (h!E-1!NGyr)"
yaxistitle = "Metallicity (Z!D!9n!3!N)"

xmax = 3.0
xmin = 0.0

ymax = 4.0
;ymin = 0.0001
ymin = 0.01

;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; average and 1 sigma
;thispsym= 3
;thiscolor= 00
;oplot, time, zmets, thick=4.0, psym=-thispsym, color=thiscolor
;oplot, time, zmets+zmetserr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
;oplot, time, zmets-zmetserr, thick=1.0, psym=-3, color=thiscolor, linestyle=1

; median and +/- 25 percentile
thispsym= 5
thiscolor= 150
oplot, time_1, zmet_50_1, thick=4.0, psym=-thispsym, color=thiscolor
oplot, time_1, zmet_25_1, thick=2.0, psym=-3, color=thiscolor, linestyle=1
oplot, time_1, zmet_75_1, thick=2.0, psym=-3, color=thiscolor, linestyle=1

xyouts, 0.25, 0.85, 'black hole', /normal, charthick=1, size=1.33, color=150

; median and +/- 25 percentile
thispsym= 5
thiscolor= 50
oplot, time_2, zmet_50_2, thick=4.0, psym=-thispsym, color=thiscolor
oplot, time_2, zmet_25_2, thick=2.0, psym=-3, color=thiscolor, linestyle=1
oplot, time_2, zmet_75_2, thick=2.0, psym=-3, color=thiscolor, linestyle=1

xyouts, 0.25, 0.80, 'no black hole', /normal, charthick=1, size=1.33, color=50


device, /close


end






;
; do the work
;---------------------------------------
pro process_metals_time, frun, $
			time= time, $
			z_avg= z_avg, $
			z_err= z_err, $
			z_median= z_median, $
			z_plus25= z_plus25, $
			z_min25= z_min25


idlistname= frun+"/hotgas_idlist.txt"

; need to load id_accounting
;   (and order by id number)
idlist_fmfile= load_id_list(idlistname)
print, "id's in file= ",n_elements(idlist_fmfile)
sidx=sort(idlist_fmfile)
idlist_fmfile= idlist_fmfile(sidx)


; this assumes they are ordered
;  0 through  x


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])



ttime= fltarr(nsnaps)

zmets=      fltarr(nsnaps)
zmetserr=   fltarr(nsnaps)

zmet_25=    fltarr(nsnaps)
zmet_50=    fltarr(nsnaps)
zmet_75=    fltarr(nsnaps)



; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"

	ok= fload_snapshot_bh(frun,i)

	        ; what time is it?
        ttime[i]= fload_time(1)


	; get current snap info
	; -----------------------
	gasz= fload_gas_metallicity(1)
	gasz= gasz/0.02
	gids= fload_gas_id(1)


	; select desired id's
	; ---------------------
	idx= intarr(n_elements(gids))
	for ii=0,n_elements(idlist_fmfile)-1 do begin
		;if idlist_fmfile[i] eq lstid then duplicates= duplicates+1
		inlst= where(gids eq idlist_fmfile[ii])
		if inlst(0) ne -1 then begin
		   idx(inlst)= 1
		endif else begin
		   print, "couldn't find id=",idlist_fmfile[ii]
		endelse
		;lstid= idlist_fmfile[i]
	endfor
	midx= where(idx eq 1)
	print, "matching gid's= ",n_elements(midx)


	idlistzs= gasz(midx)


	; average and 1 sigma
	; --------------------
	z_moment= moment(idlistzs)
	zmets[i]= z_moment[0]
	zmetserr[i]= sqrt(z_moment[1])


	; median and +/- 25%
	; -------------------
	sortidx= sort(idlistzs)
	zs_ordered= idlistzs(sortidx)
	halfidx= long(n_elements(sortidx)/2.0)
	two5= long(n_elements(sortidx)/4.0)
	seven5= halfidx+two5

	print, "median= ", zs_ordered[halfidx],"   +/- 25%= ", zs_ordered[two5], zs_ordered[seven5]
	zmet_50[i]= zs_ordered[halfidx]
	zmet_25[i]= zs_ordered[two5]
	zmet_75[i]= zs_ordered[seven5]


	print, "Z= ", zmets[i], "   +/- ",zmetserr[i]

endfor

time= ttime
z_avg= zmets
z_err= zmetserr
z_median= zmet_50
z_plus25= zmet_25
z_min25= zmet_75


end






;--------------------------------------------------------------------------------



;
; Generate plot showing the average gas
; SFR as a function of time for all
; x-ray gas in the merger remnant.
;----------------------------------------------
pro plot_gas_sfr, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "plot_gas_sfr, junk"
   return
endif

;  1
; ---
frun= "/raid4/tcox/vc3vc3e_2"
;frun= "/raid4/tcox/vc3vc3e"
;frun= "/raid4/tcox/vc3vc3e_no"

process_sfr_time, frun, time=time, sfrtot=sfrtot

time_1=       time
sfr_tot_1=    sfrtot


;  2
; ---
;frun= "/raid4/tcox/vc3vc3e_2"
;frun= "/raid4/tcox/vc3vc3e"
frun= "/raid4/tcox/vc3vc3e_no"

process_sfr_time, frun, time=time, sfrtot=sfrtot

time_2=       time
sfr_tot_2=    sfrtot





;----------------
; Print it up
;----------------

;filename=frun+'/newstarage.eps'
filename='sfrenrich.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

xaxistitle = "Time (h!E-1!NGyr)"
yaxistitle = "SFR (M!D!9n!3!N yr!E-1!N)"

xmax = 3.0
xmin = 0.0

ymax = 100.0
;ymin = 0.0001
ymin = 0.001

;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; median and +/- 25 percentile
thispsym= 5
thiscolor= 150
oplot, time_1, sfr_tot_1, thick=4.0, psym=-thispsym, color=thiscolor
xyouts, 0.25, 0.85, 'black hole', /normal, charthick=1, size=1.33, color=150

; median and +/- 25 percentile
thispsym= 6
thiscolor= 50
oplot, time_2, sfr_tot_2, thick=4.0, psym=-thispsym, color=thiscolor
xyouts, 0.25, 0.80, 'no black hole', /normal, charthick=1, size=1.33, color=50


device, /close


end






;
; do the work
;---------------------------------------
pro process_sfr_time, frun, time= time, sfrtot=sfrtot


idlistname= frun+"/hotgas_idlist.txt"

; need to load id_accounting
;   (and order by id number)
idlist_fmfile= load_id_list(idlistname)
print, "id's in file= ",n_elements(idlist_fmfile)
sidx=sort(idlist_fmfile)
idlist_fmfile= idlist_fmfile(sidx)


; this assumes they are ordered
;  0 through  x


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])



ttime= fltarr(nsnaps)

sfr_tot=      fltarr(nsnaps)



; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"

	ok= fload_snapshot_bh(frun,i)

	        ; what time is it?
        ttime[i]= fload_time(1)


	; get current snap info
	; -----------------------
	sfr= fload_gas_sfr(1)

	gids= fload_gas_id(1)


	; select desired id's
	; ---------------------
	idx= intarr(n_elements(gids))
	for ii=0,n_elements(idlist_fmfile)-1 do begin
		;if idlist_fmfile[i] eq lstid then duplicates= duplicates+1
		inlst= where(gids eq idlist_fmfile[ii])
		if inlst(0) ne -1 then begin
		   idx(inlst)= 1
		endif else begin
		   print, "couldn't find id=",idlist_fmfile[ii]
		endelse
		;lstid= idlist_fmfile[i]
	endfor
	midx= where(idx eq 1)
	print, "matching gid's= ",n_elements(midx)


	idlistsfr= sfr(midx)

	sfr_tot[i]= total(idlistsfr)

	print, "total sfr= ", total(sfr)
	print, "   hg sfr= ", sfr_tot[i]

endfor

time= ttime
sfrtot= sfr_tot

end






;--------------------------------------------------------------------------------




;
; Generate plot showing the average radius
; as a function of time for all
; x-ray emitting gas in the merger remnant.
;----------------------------------------------
pro plot_gas_r, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "plot_gas_r, junk"
   return
endif

;frun= "/raid4/tcox/vc3vc3e_2"
;frun= "/raid4/tcox/vc3vc3e"
frun= "/raid4/tcox/vc3vc3e_no"

process_radial_info, frun, $
		r_avg= r_avg, $
		r_err= r_err, $
		r_median= r_median, $
		r_plus25= r_plus25, $
		r_min25= r_min25


rs= r_avg
rerr= r_err
r_50= r_median
r_75= r_plus25
r_25= r_min25





;----------------
; Print it up
;----------------

filename='renrich.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

xaxistitle = "Time (h!E-1!NGyr)"
yaxistitle = "Radius (h!E-1!Nkpc)"

xmax = 3.0
xmin = 0.0

ymax = 100.0
;ymin = 0.0001
ymin = 0.0

;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; average and 1 sigma
thispsym= 3
thiscolor= 00
oplot, time, rs, thick=4.0, psym=-thispsym, color=thiscolor
oplot, time, rs+rerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
oplot, time, rs-rerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1

; median and +/- 25 percentile
thispsym= 3
thiscolor= 150
oplot, time, r_50, thick=3.0, psym=-thispsym, color=thiscolor
oplot, time, r_25, thick=2.0, psym=-3, color=thiscolor, linestyle=1
oplot, time, r_75, thick=2.0, psym=-3, color=thiscolor, linestyle=1


xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0

device, /close


end





;
;
;-----------------------------------------
pro process_radial_info, frun, $
		r_avg= rs, $
		r_err= rerr, $
		r_median= r_50, $
		r_plus25= r_75, $
		r_min25= r_25

idlistname= frun+"/hotgas_idlist.txt"

; need to load id_accounting
;   (and order by id number)
idlist_fmfile= load_id_list(idlistname)
print, "id's in file= ",n_elements(idlist_fmfile)
sidx=sort(idlist_fmfile)
idlist_fmfile= idlist_fmfile(sidx)


; this assumes they are ordered
;  0 through  x


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])



time= fltarr(nsnaps)

rs=      fltarr(nsnaps)
rerr=   fltarr(nsnaps)

r_25=    fltarr(nsnaps)
r_50=    fltarr(nsnaps)
r_75=    fltarr(nsnaps)



; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"

	ok= fload_snapshot_bh(frun,i)

	        ; what time is it?
        time[i]= fload_time(1)


	; get current snap info
	; -----------------------
	radius= fload_gas_xyz('r')
	gids= fload_gas_id(1)


	; select desired id's
	; ---------------------
	idx= intarr(n_elements(gids))
	for ii=0,n_elements(idlist_fmfile)-1 do begin
		;if idlist_fmfile[i] eq lstid then duplicates= duplicates+1
		inlst= where(gids eq idlist_fmfile[ii])
		if inlst(0) ne -1 then begin
		   idx(inlst)= 1
		endif else begin
		   print, "couldn't find id=",idlist_fmfile[ii]
		endelse
		;lstid= idlist_fmfile[i]
	endfor
	midx= where(idx eq 1)
	print, "matching gid's= ",n_elements(midx)


	; now we can select the radii
	; we want - and only those
	; ------------------------------
	idlistrs= radius(midx)


	; average and 1 sigma
	; --------------------
	r_moment= moment(idlistrs)
	rs[i]= r_moment[0]
	rerr[i]= sqrt(r_moment[1])


	; median and +/- 25%
	; -------------------
	sortidx= sort(idlistrs)
	rs_ordered= idlistrs(sortidx)
	halfidx= long(n_elements(sortidx)/2.0)
	two5= long(n_elements(sortidx)/4.0)
	seven5= halfidx+two5

	print, "median= ", rs_ordered[halfidx],"   +/- 25%= ", rs_ordered[two5], rs_ordered[seven5]
	r_50[i]= rs_ordered[halfidx]
	r_25[i]= rs_ordered[two5]
	r_75[i]= rs_ordered[seven5]


	print, "R= ", rs[i], "   +/- ",rerr[i]

endfor

r_avg= rs
r_err= rerr
r_median= r_50
r_plus25= r_75
r_min25= r_25


end






;--------------------------------------------------------------------------------






;
; Generate plot showing the average gas
; metallicity as a function of time for all
; x-ray gas in the merger remnant.
;----------------------------------------------
pro plot_gas_zfrac, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "plot_gas_zfrac, junk"
   return
endif

;frun= "/raid4/tcox/vc3vc3e_2"
;frun= "/raid4/tcox/vc3vc3e"
frun= "/raid4/tcox/vc3vc3e_no"


; -------------------------------------------------
; need to load id_accounting
;   (and order by id number)
idlistname= frun+"/hotgas_idlist.txt"
idlist_fmfile= load_id_list(idlistname)
print, "id's in file= ",n_elements(idlist_fmfile)
sidx=sort(idlist_fmfile)
idlist_fmfile= idlist_fmfile(sidx)
; -------------------------------------------------
; need to load id_accounting
;   (and order by id number)
idlistname= frun+"/wind_idlist.txt"
idlist_fmfile_2= load_id_list(idlistname)
print, "id's in file= ",n_elements(idlist_fmfile_2)
sidx=sort(idlist_fmfile_2)
idlist_fmfile_2= idlist_fmfile_2(sidx)
; -------------------------------------------------


; this assumes they are ordered
;  0 through  x


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])



time= fltarr(nsnaps)

metals_gas=  fltarr(nsnaps)
metals_hgas= fltarr(nsnaps)
metals_wgas= fltarr(nsnaps)

metals_ns=   fltarr(nsnaps)



; ----------------------------------------
; This part loops through the snapshots
; and compiles hot gas information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"

	ok= fload_snapshot_bh(frun,i)

	        ; what time is it?
        time[i]= fload_time(1)


	; get current snap info
	; -----------------------
	gasz= fload_gas_metallicity(1)
	gasm= fload_gas_mass(1)

	; total metals in gas
	metals_gas[i]= total(gasz*gasm)

	nsz= fload_newstars_z(1)
	nsm= fload_newstars_mass(1)

	; total metals in stars
	metals_ns[i]= total(nsz*nsm)

	gids= fload_gas_id(1)



	; Hot X-ray Gas
	; ---------------------
	; 
	; select desired id's 
	idx= intarr(n_elements(gids))
	for ii=0,n_elements(idlist_fmfile)-1 do begin
		;if idlist_fmfile[i] eq lstid then duplicates= duplicates+1
		inlst= where(gids eq idlist_fmfile[ii])
		if inlst(0) ne -1 then begin
		   idx(inlst)= 1
		endif else begin
		   print, "couldn't find id=",idlist_fmfile[ii]
		endelse
		;lstid= idlist_fmfile[i]
	endfor
	midx= where(idx eq 1)
	print, "matching gid's= ",n_elements(midx)

	metals_hgas[i]= total(gasz(midx)*gasm(midx))


        ; Wind (E>=0) Gas
        ; ---------------------
        ; 
        ; select desired id's 
        idx= intarr(n_elements(gids))
        for ii=0,n_elements(idlist_fmfile_2)-1 do begin
                inlst= where(gids eq idlist_fmfile_2[ii])
                if inlst(0) ne -1 then begin
                   idx(inlst)= 1
                endif else begin
                   print, "couldn't find id=",idlist_fmfile_2[ii]
                endelse
        endfor
        midx= where(idx eq 1)
        print, "matching gid's= ",n_elements(midx)

        metals_wgas[i]= total(gasz(midx)*gasm(midx))


endfor






;----------------
; Print it up
;----------------

filename='zfenrich.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

xaxistitle = "Time (h!E-1!NGyr)"
yaxistitle = "Metal Mass (10!E10!N M!D!9n!3!N)"

xmax = 3.0
xmin = 0.0

ymax = 0.1
ymin = 0.0001

; make it metal fraction
;metalf= 1
metalf= 0
if metalf eq 1 then begin
     yaxistitle = "Metal Fraction (%)"
     ymax= 1.0
     ymin= 0.0
     ; ylog needs manual setting
     ; yaxis exp needs manual setting
endif

;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


if metalf eq 1 then begin
	metals_total= metals_gas+metals_ns
	idx= where(metals_total le 0.0)
	if idx(0) eq -1 then metals_total(idx)= 1.0

	; metals in the gas
	thispsym= 3
	thiscolor= 00
	oplot, time, metals_gas/metals_total, thick=4.0, psym=-thispsym, color=thiscolor

	; metals in the gas (hot x-ray gas)
	thispsym= 4
	thiscolor= 150
	oplot, time, metals_hgas/metals_total, thick=4.0, psym=-thispsym, color=thiscolor

	; metals in the new stars
	thispsym= 3
	thiscolor= 50
	oplot, time, metals_ns/metals_total, thick=3.0, psym=-thispsym, color=thiscolor
endif else begin
	; metals in the gas
	thispsym= 3
	thiscolor= 00
	oplot, time, metals_gas, thick=4.0, psym=-thispsym, color=thiscolor
	xyouts, 0.65, 0.53, 'gas - total', /normal, charthick= 2.0, size=1.2, color=thiscolor

	; metals in the gas (hot x-ray gas)
	thispsym= 4
	thiscolor= 150
	oplot, time, metals_hgas, thick=4.0, psym=-thispsym, color=thiscolor
	xyouts, 0.65, 0.32, 'gas - X-ray', /normal, charthick= 2.0, size=1.2, color=thiscolor

	; metals in the gas (wind e>=0 gas)
	thispsym= 6
	thiscolor= 100
	oplot, time, metals_wgas, thick=4.0, psym=-thispsym, color=thiscolor
	xyouts, 0.54, 0.20, '(none)', /normal, charthick= 2.0, size=1.2, color=thiscolor
	xyouts, 0.65, 0.20, 'gas - unbound', /normal, charthick= 2.0, size=1.2, color=thiscolor

	; metals in the new stars
	thispsym= 3
	thiscolor= 50
	oplot, time, metals_ns, thick=3.0, psym=-thispsym, color=thiscolor
	xyouts, 0.65, 0.83, 'new stars', /normal, charthick= 2.0, size=1.2, color=thiscolor
endelse

xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0

device, /close


end








;--------------------------------------------------------------------------------




;
; Generate plot showing the age of each gas
; particles, as measured by the time when it
; was last in a star forming region.
; (still only for x-ray gas in the merger remnant)
;----------------------------------------------
pro plot_gas_zages, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "plot_gas_zages, junk"
   return
endif

;frun= "/raid4/tcox/vc3vc3e_2"
frun= "/raid4/tcox/vc3vc3e"
;frun= "/raid4/tcox/vc3vc3e_no"
process_halfz_lastsf, frun, halfz=halfz, lastsf=lastsf
time_halfZ1= halfz
time_lastSF1= lastsf

;frun= "/raid4/tcox/vc3vc3e_2"
;frun= "/raid4/tcox/vc3vc3e"
frun= "/raid4/tcox/vc3vc3e_no"
process_halfz_lastsf, frun, halfz=halfz, lastsf=lastsf
time_halfZ2= halfz
time_lastSF2= lastsf




;-------------------
; Print it up - Z
;-------------------

filename='halfz_hist.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

xaxistitle = "Half Z Time (h!E-1!NGyr)"
yaxistitle = "Number"

xmax = 3.0
xmin = 0.0

;ymax = 2000
ymax = 1.2
ymin = 0

;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	yticks= 1, $
        ;ytickformat='exp_label', $
        ytickformat='(a1)', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


levels= 20.0
step= (xmax-xmin)/levels
bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
bins=[bins,xmax]

;  1
; ---
hist_time_Z= histogram(time_halfZ1, binsize=step, max=xmax, min=xmin)
print, "max hist_time_Z #1=", max(hist_time_Z)
normalization= max(hist_time_Z)
hist_time_Z= 1.0*hist_time_Z/normalization
oplot, bins, hist_time_Z, psym=10, color=150, thick= 4.0

xyouts, 0.65, 0.85, 'black hole', /normal, charthick=1, size=1.33, color=150

;  2
; ---
hist_time_Z= histogram(time_halfZ2, binsize=step, max=xmax, min=xmin)
print, "max hist_time_Z #2=", max(hist_time_Z)
hist_time_Z= 1.0*hist_time_Z/normalization
oplot, bins, hist_time_Z, psym=10, color=50, thick= 4.0, linestyle= 1

xyouts, 0.65, 0.80, 'no black hole', /normal, charthick=1, size=1.33, color=50

device, /close




;-------------------
; Print it up - SF
;-------------------

filename='lastSF_hist.eps'
        
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

xaxistitle = "Last SF Time (h!E-1!NGyr)"
yaxistitle = "Number"
        
xmax = 3.0
xmin = 0.0
        
ymax = 1.2
;ymax = 2000
ymin = 0
        
;---------------------------
                   
!p.position= [0.18, 0.15, 0.95, 0.95]
                   
plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;/ylog, $  
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        ytickformat='(a1)', $
	yticks= 1, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata 
        

levels= 20.0
step= (xmax-xmin)/levels
bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
bins=[bins,xmax]

;  1
; ---
hist_time_SF= histogram(time_lastSF1, binsize=step, max=xmax, min=xmin)
print, "Max hist_time_SF=", max(hist_time_SF)
normalization= max(hist_time_SF)
hist_time_SF= 1.0*hist_time_SF/normalization
oplot, bins, hist_time_SF, psym=10, color=150, thick= 4.0

xyouts, 0.65, 0.85, 'black hole', /normal, charthick=1, size=1.33, color=150

;  2
; ---
hist_time_SF= histogram(time_lastSF2, binsize=step, max=xmax, min=xmin)
print, "Max hist_time_SF=", max(hist_time_SF)
hist_time_SF= 1.0*hist_time_SF/normalization
oplot, bins, hist_time_SF, psym=10, color=50, thick= 4.0, linestyle= 1

xyouts, 0.65, 0.80, 'no black hole', /normal, charthick=1, size=1.33, color=50

device, /close




end








;
; actually do the work
;-------------------------
pro process_halfz_lastsf, frun, halfz=halfz, lastsf=lastsf


idlistname= frun+"/hotgas_idlist.txt"

; need to load id_accounting
;   (and order by id number)
idlist_fmfile= load_id_list(idlistname)
print, "id's in file= ",n_elements(idlist_fmfile)
sidx=sort(idlist_fmfile)
idlist_fmfile= idlist_fmfile(sidx)


; this assumes they are ordered
;  0 through  x


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])



time= fltarr(nsnaps)

;metals_gas=  fltarr(nsnaps)
;metals_hgas= fltarr(nsnaps)
;metals_ns=   fltarr(nsnaps)


; arrays to store last time when
; gas particle was SFing and/or
; reached 1/2 its metallicity
Ngas= long(n_elements(idlist_fmfile))
time_lastSF= fltarr(Ngas)
time_halfZ= fltarr(Ngas)

store_finalZ= fltarr(Ngas)


; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"

	snapnum= nsnaps-1-i
	ok= fload_snapshot_bh(frun,snapnum)

        ; what time is it?
	ttime= fload_time(1)


	; get current snap info
	; -----------------------
	gasz= fload_gas_metallicity(1)
	gasm= fload_gas_mass(1)

	rho= fload_gas_rho(1)

	gids= fload_gas_id(1)


	; select desired id's
	; ---------------------
	idx= intarr(n_elements(gids))
	for ii=0,n_elements(idlist_fmfile)-1 do begin

		inlst= where(gids eq idlist_fmfile[ii])
		if inlst(0) ne -1 then begin
		   idx(inlst)= 1

		   ; OK, we have a matching element
		   ; -------------------------------

		   ; on first pass, set final Z
		   if i eq 0 then store_finalZ[ii]= gasz(inlst)

		   ; is the metallicity below half it's final value?
		   if (gasz(inlst) lt (0.5*store_finalZ[ii]) and time_halfZ[ii] le 0.001) then time_halfZ[ii]= ttime

		   ; is the density below SF densities?
		   if (rho(inlst) gt (0.000854924) and time_lastSF[ii] le 0.001) then time_lastSF[ii]= ttime

		endif else begin
		   print, "couldn't find id=",idlist_fmfile[ii]
		endelse
	endfor

	midx= where(idx eq 1)
	print, "matching gid's= ",n_elements(midx)


endfor



halfz= time_halfZ
lastsf= time_lastSF


end






;--------------------------------------------------------------------------------



