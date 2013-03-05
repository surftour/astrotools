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
;  become unbound from the system
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
;frun= "/raid4/tcox/vc3vc3e_no"
;snapnum= 107
;frun= "/raid4/tcox/vc3vc3e"
;snapnum= 30


idlistname= frun+"/wind_idlist.txt"
;idlistname= frun+"/wind_idlist_atfirstpassage.txt"


ok= fload_snapshot_bh(frun,snapnum)


gids= fload_gas_id(1)

ke= fload_gas_energy(1,/kinetic)
pe= fload_gas_energy(1,/potential)
the= fload_gas_energy(1,/thermal)
energy= ke + pe + the



; take only xray luminous gas particles
; --------------------------------------
idx= where(energy ge 0)
if idx(0) ne -1 then begin
	en= energy(idx)
	gids= gids(idx)
endif
print, 'there are ', n_elements(idx), ' non-zero gas particles'


; now order these
; -----------------
eindx= sort(en)
gids= gids(eindx)
en= en(eindx)


idlist= gids

print, 'writing ', n_elements(idlist), ' ids'


;  Write id list to a file
; ----------------------------
cmt= 'unbound gas ordered by energy'
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
for i=0,n_elements(idlist)-1 do printf, unit, idlist[i]

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

;frun= "/raid4/tcox/vc3vc3e_2"
frun= "/raid4/tcox/vc3vc3e"

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
        time[i]= fload_time(1)


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
thispsym= 3
thiscolor= 00
oplot, time, zmets, thick=4.0, psym=-thispsym, color=thiscolor
oplot, time, zmets+zmetserr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
oplot, time, zmets-zmetserr, thick=1.0, psym=-3, color=thiscolor, linestyle=1

; median and +/- 25 percentile
thispsym= 3
thiscolor= 150
oplot, time, zmet_50, thick=3.0, psym=-thispsym, color=thiscolor
oplot, time, zmet_25, thick=2.0, psym=-3, color=thiscolor, linestyle=1
oplot, time, zmet_75, thick=2.0, psym=-3, color=thiscolor, linestyle=1


xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0

device, /close


end








;--------------------------------------------------------------------------------






