
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Cumulative Mass Profile
;   ----------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro cum_mass, frun, smoothlen=smoothlen, $
			filename=filename, $
			snapnum=snapnum, $
			center=center
        

if not keyword_set(frun) then begin
   print, "  "
   print, " cum_mass, frun, smoothlen=smoothlen, "
   print, "             filename=filename, "
   print, "             snapnum=snapnum"
   print, "  "
   return
endif   
        
if not keyword_set(filename) then filename="cummass.eps"
if not keyword_set(smoothlen) then smoothlen=0.1


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


; ----------
;  snapshot
; ----------
if not keyword_set(snapnum) then snapnum= 25


; ----------------------------- 
; set up constants      
; -----------------------------
                        
bins = 50
                
;xmax = 300
xmax = 10.0
;xmin = 0.1
xmin= 0.0

ymax = 5e+12
ymin = 1e+8

xaxistitle= '!6R (kpc)'
;yaxistitle= 'M!Dgas!N (M!D!9n!3!N)'
yaxistitle= '!6Mass (M!D!9n!6!N)'


;---------------------------
;  Print it
;---------------------------
   
!p.position= [0.18, 0.15, 0.98, 0.98]
   
plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
        ;/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata

; ------------------------
;  determine gas profile
; ------------------------

processandplot_onemassprofile, frun, snapnum, xmin, xmax, bins, center=center



; -----------------
;  Plot Extras
; -----------------

xyouts, 0.65, 0.18, 'yuexing, grid '+strcompress(string(snapnum),/remove_all), /normal, color= 0, charthick=2.0, size=1.2


; smoothing length
;x=[smoothlen,smoothlen]
;y=[ymin,ymax]
;oplot, x, y, linestyle=1, color= 0


;--------------------------------------
;--------------------------------------

device, /close



end





;==================================================================================






; actually do the work
; --------------------------
pro processandplot_onemassprofile, frun, snapnum, xmin, xmax, bins, $
					center=center, xlog=xlog


ok=fload_snapshot_bh(frun,snapnum)
if not keyword_set(center) then center=fload_center_alreadycomp(1)


print, "using center= ", center

radius=fload_allstars_xyz('r',center=center) / 0.7
mass=fload_allstars_mass(1)
print, "total stellar mass= ", total(mass)
mass= 1.0d+10 * mass / 0.7
process_cummass_profile, radius, mass, xmin, xmax, bins, xs, ys
print, "stars max= ", ys[bins-1]
oplot, xs, ys, psym=-3, linestyle=0, color= 150, thick=4.0

xyouts, 0.65, 0.88, 'stars', /normal, color= 150, charthick=2.0, size=1.2


radius=fload_gas_xyz('r',center=center) / 0.7
mass=fload_gas_mass(1)
print, "total gas mass= ", total(mass)
mass= 1.0d+10 * mass / 0.7
process_cummass_profile, radius, mass, xmin, xmax, bins, xs, ys
print, "gas max= ", ys[bins-1]
oplot, xs, ys, psym=-3, linestyle=0, color= 50, thick=4.0

xyouts, 0.65, 0.41, 'gas', /normal, color= 50, charthick=2.0, size=1.2



radius=fload_halo_xyz('r',center=center) / 0.7
mass=fload_halo_mass(1)
print, "total halo mass= ", total(mass)
mass= 1.0d+10 * mass / 0.7
process_cummass_profile, radius, mass, xmin, xmax, bins, xs, ys
print, "halo max= ", ys[bins-1]
oplot, xs, ys, psym=-3, linestyle=2, color= 100, thick=4.0

xyouts, 0.65, 0.65, 'dark matter', /normal, color= 100, charthick=2.0, size=1.2




end





;==================================================================================




; --------------------------
; process the cumulative profile

pro process_cummass_profile, x, y, xmin, xmax, bins, $
				xs, ys, $
				xlog=xlog

xs = fltarr(bins)
ys= fltarr(bins)

tempxmax= xmax
tempxmin= xmin
if keyword_set(xlog) then begin
	;tempxmin= alog10(xmin)
	;tempxmax= alog10(xmax)
	x= alog10(x)
endif

binsize = float((tempxmax-tempxmin))/bins

for i=0, bins-1 do begin
    currentr= tempxmin + (i+1)*binsize
    idx= where(x le currentr)
    xs[i]= currentr
    if idx(0) ne -1 then ys[i]= total(y(idx)) else ys[i]= 1.0d-5
endfor



end




;==================================================================================



pro do_desika, frun

bins = 50
xmax = 10.0
xmin= 0.0

  if not keyword_set(frun) then begin
        print, "  "
        print, " do_desika, frun"
        print, "  "
        return
   endif

   bins = 50
   xmax = 10.0
   xmin= 0.0

   snapinterval= 1
   ;snapinterval= 10

   starti= 0L

   spawn, "/bin/ls "+frun+"/snap* | wc ",result                                  
   endi=long(result[0])-1


   print, "frun= ", frun
   print, "starti= ", starti
   print, "endi= ", endi
   print, "xmax= ", xmax
   print, "xmin= ", xmin
   print, "bins= ", bins

   for i=starti,endi do begin

     thisi= i
     if ((thisi mod snapinterval) eq 0) then begin

        ok=fload_snapshot_bh(frun, thisi, /nopot_in_snap)

        bhid= fload_blackhole_id(1)
        ;bhid= bhid[0]
        ;bhid= bhid[1]
        ;bhid= 200001L
        ;bhid= 280002L   ; used for z3/b4e
        ;bhid= 400002L   ; used for ds/vc3vc3e_2
        bhid= long(max(bhid))   ; uses minimum bhid - then you don't need to know the BH ID
        center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)

        print, "Blackhole ID: ", bhid
        print, "Blackhole center: ", center_bh

        ;print, "make_cummass_file, "+frun+", "+strcompress(string(thisi),/remove_all)
	make_cummass_file, frun, thisi, xmin, xmax, bins, center=center_bh, /loadedsnap
     endif

   endfor

end




;==================================================================================




pro do_desika_1, junk

bins = 50

xmax = 10.0
xmin= 0.0

;frun= "/raid4/tcox/bs/b3e_2"  &  starti= 71  &  endi= 99
;frun= "/raid4/tcox/bs/b3e"  &  starti= 11  &  endi= 12
;frun= "/raid4/tcox/ds/d5e2_2"  &  starti= 79  &  endi= 129

;frun= "/raid4/tcox/vc3vc3e_2"  &  starti= 38  &  endi= 72
;frun= "/raid4/tcox/vc3vc3e_no"  &  starti= 38  &  endi= 72
;frun= "/raid4/tcox/vc3vc3h_2"  &  starti= 38  &  endi= 72
;frun= "/raid4/tcox/sbw/sb10"  &  starti= 38  &  endi= 72
frun= "/raid4/tcox/sbw/sb10BH"  &  starti= 38  &  endi= 72

for i=starti, endi do begin
	make_cummass_file, frun, i, xmin, xmax, bins
endfor

make_cummass_file, frun, 82, xmin, xmax, bins
make_cummass_file, frun, 87, xmin, xmax, bins
make_cummass_file, frun, 97, xmin, xmax, bins

end





;==================================================================================




pro do_desika_2, junk

bins = 50

xmax = 10.0
xmin= 0.0


dirlist=["/raid4/tcox/brantsruns/full_model/z6/v4g1q0z6_v4g1q0z6_p", $
         "/raid4/tcox/brantsruns/full_model/z6/v4g1q1z6_v4g1q1z6_p", $
         "/raid4/tcox/brantsruns/full_model/z6/v4g2q0z6_v4g2q0z6_p", $
         "/raid4/tcox/brantsruns/full_model/z6/v4g2q1z6_v4g2q1z6_p", $
         "/raid4/tcox/brantsruns/full_model/z6/v5g1q0z6_v5g1q0z6_p", $
         "/raid4/tcox/brantsruns/full_model/z6/v5g1q1z6_v5g1q1z6_p", $
         "/raid4/tcox/brantsruns/full_model/z6/v5g2q0z6_v5g2q0z6_p", $
         "/raid4/tcox/brantsruns/full_model/z6/v5g2q1z6_v5g2q1z6_p", $
         "/raid4/tcox/brantsruns/full_model/z6/v6g1q0z6_v6g1q0z6_p", $
         "/raid4/tcox/brantsruns/full_model/z6/v6g1q1z6_v6g1q1z6_p", $
         "/raid4/tcox/brantsruns/full_model/z6/v6g2q0z6_v6g2q0z6_p", $
         "/raid4/tcox/brantsruns/full_model/z6/v6g2q1z6_v6g2q1z6_p"]

nd= n_elements(dirlist)-1


for i=0, nd do begin

        print, "==============================================================="
        print, "==============================================================="
        print, dirlist[i]

	make_cummass_file, dirlist[i], 20, xmin, xmax, bins

endfor


end






;==================================================================================





; produce a cumulative mass file
; -------------------------------
pro make_cummass_file, frun, snapnum, xmin, xmax, bins, $
					center=center, xlog=xlog, $
					filename=filename, loadedsnap=loadedsnap


if not keyword_set(loadedsnap) then begin
	ok=fload_snapshot_bh(frun,snapnum)
endif

exts='0000'
exts=exts+strcompress(string(snapnum),/remove_all)
elen= strlen(exts)
exts= strmid(exts,elen-3,3)

if not keyword_set(center) then center=fload_center_alreadycomp(1)
if not keyword_set(filename) then filename= frun+"/cummass/cummass_"+exts+".txt"


print, "using center= ", center
radius=fload_allstars_xyz('r',center=center) / 0.7
mass=fload_allstars_mass(1)
print, "total stellar mass= ", total(mass)
mass= 1.0d+10 * mass / 0.7
process_cummass_profile, radius, mass, xmin, xmax, bins, xs, ys
print, "stars max= ", ys[bins-1]
xs_stars= xs
ys_stars= ys



print, "using center= ", center
radius=fload_gas_xyz('r',center=center) / 0.7
mass=fload_gas_mass(1)
print, "total gas mass= ", total(mass)
mass= 1.0d+10 * mass / 0.7
process_cummass_profile, radius, mass, xmin, xmax, bins, xs, ys
print, "gas max= ", ys[bins-1]
xs_gas= xs
ys_gas= ys



print, "using center= ", center
radius=fload_halo_xyz('r',center=center) / 0.7
mass=fload_halo_mass(1)
print, "total halo mass= ", total(mass)
mass= 1.0d+10 * mass / 0.7
process_cummass_profile, radius, mass, xmin, xmax, bins, xs, ys
print, "halo max= ", ys[bins-1]
xs_dm= xs
ys_dm= ys


ys_total= ys_dm
for i=0, bins-1 do begin
   ys_total[i]= ys_dm[i] + ys_gas[i] + ys_stars[i]
endfor

ys_total= alog10(ys_total)
ys_dm= alog10(ys_dm)
ys_stars= alog10(ys_stars)
ys_gas= alog10(ys_gas)



centerlbl= string(center[0])+','+string(center[1])+','+string(center[2])



; Xrays FILE
; -----------
print, "writing: ", filename
openw, 1, filename, ERROR=err

printf, 1, "# "+filename+' frun='+frun+'  snapnum='+strcompress(string(snapnum),/remove_all)+'  center='+centerlbl
printf, 1, "#  "
printf, 1, "#  "
printf, 1, "# radius      <----- log(Mass) in M_solar -----> "
printf, 1, "#  (kpc)      Total       DM     stars       gas"
for i=0,bins-1 do begin
        printf, 1, FORMAT= '(F8.3,"    ", 4(F8.5,"  "))', $
                xs_dm[i], ys_total[i], ys_dm[i], ys_stars[i], ys_gas[i]
endfor
close, 1



end


