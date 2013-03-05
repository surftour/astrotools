;-------------------------------------------------------------
;
;
;-------------------------------------------------------------

;
; this routine measures the age histogram of stars in the tidal
; dwarfs
;
;

pro tdwarfs_nsage_hist, junk, filename=filename


if not keyword_set(junk) then begin
   print, "  "
   print, "tdwarfs_nsage_hist, junk, filename=filename"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='td_agehist.eps'





;
;  histogram parameters
;
; --------------------------------------
age_min = 0.0
age_max = 4.25

nbins= 50






; ----------------------------
;
; plot the age histogram
;
; ----------------------------



initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


xmax= age_max
xmin= age_min

;ymax= 190
ymax= 90
ymin= 0

yaxistitle="!6Number"

;xaxistitle = "!6New Star Age (!8h!6!E-1!N Gyr)"
;xaxistitle = "!6New Star Age (Gyr)"
xaxistitle = "!6New Star Formation Time (Gyr)"


!p.position= [0.16, 0.14, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; need snapshot open
;ok= fload_snapshot("/home/tcox/Sbc10x", 60)

;process_one_group, "/home/tcox/Sbc10x/idlists", 3, age_max, age_min, nbins, linecolor= 150
process_one_group, "/home/tcox/Sbc10x/idlists", 5, age_max, age_min, nbins, linecolor= 50
process_one_group, "/home/tcox/Sbc10x/idlists", 13, age_max, age_min, nbins, linecolor= 150
;process_one_group, "/home/tcox/Sbc10x/idlists", 5, age_max, age_min, nbins, linecolor= 200
;process_one_group, "/home/tcox/Sbc10x/idlists", 6, age_max, age_min, nbins, linecolor= 170
process_one_group, "/home/tcox/Sbc10x/idlists", 7, age_max, age_min, nbins, linecolor= 0
;process_one_group, "/home/tcox/Sbc10x/idlists", 8, age_max, age_min, nbins, linecolor= 110
;process_one_group, "/home/tcox/Sbc10x/idlists", 9, age_max, age_min, nbins, linecolor= 50
process_one_group, "/home/tcox/Sbc10x/idlists", 10, age_max, age_min, nbins, linecolor= 200
process_one_group, "/home/tcox/Sbc10x/idlists", 15, age_max, age_min, nbins, linecolor= 100






; we're done
; ------------
device, /close


end










; -------------------------------------------------------------------------------


pro process_one_group, groupiddir, groupnum, age_max, age_min, nbins, $
				linecolor=linecolor, $
				normalized=normalized, $
				fit_powerlaw=fit_powerlaw




;
;  Now, read the group information, and
;   compute the histogram
; --------------------------------------

grp_center= fltarr(3)

idlistfilename= groupiddir+'/group_'+strcompress(string(groupnum), /remove_all)+'_idlist.txt'
;idlist_fmfile= load_id_list(idlistname)
; --------------------------------------
;spawn, 'wc '+idlistfilename, result
;lines= long(result)
;idlist= lonarr(lines-5)

openr, 1, idlistfilename

textjunk= ''
readf, 1, textjunk
readf, 1, textjunk
readf, 1, textjunk
readf, 1, textjunk
readf, 1, textjunk
        tempjunk= strsplit(textjunk,/extract,count=count)
        grp_center(0)= float(tempjunk(3))
        grp_center(1)= float(tempjunk(4))
        grp_center(2)= float(tempjunk(5))
close, 1


print, "group center= ", grp_center(0), grp_center(1), grp_center(2)

newstarage= fload_newstars_age(1)/0.7
ns_x= fload_newstars_xyz('x', center=grp_center)
ns_y= fload_newstars_xyz('y', center=grp_center)
ns_z= fload_newstars_xyz('z', center=grp_center)

r= sqrt(ns_x*ns_x + ns_y*ns_y + ns_z*ns_z)

idx= where(r lt 6.0)
if idx(0) eq -1 then begin
	print, " "
	print, " WARNING: couldn't find any new stars near this group center"
	print, " "
	return
endif

print, "number of new stars near group= ", n_elements(idx)


; make the age histogram
;--------------------------

; compute histogram
m_bin= (age_max - age_min)/nbins
age_hist = histogram(newstarage(idx),binsize=m_bin,min=age_min,max=age_max)

; make sure it extends to 0
age_hist= [age_hist(0), age_hist, age_hist(nbins-1)]

print, "largest group= ", max(age_hist)

; define x-axis 
binsx= (age_max-age_min)*findgen(n_elements(age_hist))/float(n_elements(age_hist)-1) + age_min
binsx= [age_min, binsx, age_max]


printthis= age_hist
if keyword_set(normalized) then printthis= age_hist/total(age_hist)



oplot, binsx, printthis, psym=10, thick=5.0, color=linecolor


; fill in histogram
; ------------------
xbins= binsx+(m_bin*0.5)          ; make x coord
xbins[0]= age_min
xbins=[xbins,xbins]
xbins= xbins(sort(xbins))

ntts= fltarr(2.*nbins+ 2)
for i=1,nbins do begin 
   ntts[2.*i-1]= printthis[i]
   ntts[2.*i]= printthis[i]
endfor

if linecolor eq 150 then begin
	polyfill, xbins, ntts, /data, color= 150, /line_fill, linestyle=0, $
                                thick=3.0, orientation=45.0
endif

if linecolor eq 50 then begin
	polyfill, xbins, ntts, /data, color= 50, /fill, thick=3.0
endif





end



; -------------------------------------------------------------------------------


