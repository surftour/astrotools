;-------------------------------------------------------------
;
; Adapted from Brant's routines
;
;-------------------------------------------------------------

;
;this routine measures a mass function from a group finder output
;
;

pro tdwarfs_massf, junk, filename=filename


if not keyword_set(junk) then begin
   print, "  "
   print, "tdwarfs_massf, junk, filename=filename"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='mass_function.eps'





;
;  Mass function histogram parameters
;
; --------------------------------------
;mass_min = 4.5		;log h^-1 Msun
mass_min = 5.5		;log h^-1 Msun
;mass_max = 7.5		;log h^-1 Msun
mass_max = 10.5		;log h^-1 Msun

nbins= 20






; ----------------------------
;
; plot the mass function
;
; ----------------------------



initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


xmax= mass_max
xmin= mass_min

ymax= 20.5
ymin= 0

yaxistitle="!6Number"

;xaxistitle = "!6Log Mass (!8h!6!E-1!N M!D!9n!6!N)"
xaxistitle = "!6Log Mass (M!D!9n!6!N)"


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




;process_one_group, "/home/tcox/Sbc", 60, mass_max, mass_min, nbins, linecolor= 50
;xyouts, 0.5, 0.8, '!6Low Res, N= 28', /normal, color= 50, size=1.5, charthick=3.0

process_one_group, "/home/tcox/Sbc10x", 60, mass_max, mass_min, nbins, linecolor= 150
xyouts, 0.5, 0.7, '!6High Res, N= 53', /normal, color= 150, size=1.5, charthick=3.0






; we're done
; ------------
device, /close


end










; -------------------------------------------------------------------------------


pro process_one_group, groupdir, snapnum, mass_max, mass_min, nbins, $
				linecolor=linecolor, $
				normalized=normalized, $
				fit_powerlaw=fit_powerlaw




;
;  Now, read the group information, and
;   compute the histogram
; --------------------------------------


num = snapnum

if num gt 999 then begin
	exts='0000'
	exts=exts+strcompress(string(num),/remove_all) ; sets exts to snap number
	exts=strmid(exts,strlen(exts)-4,4)
endif else begin
	exts='000'
	exts=exts+strcompress(string(num),/remove_all) ; sets exts to snap number
	exts=strmid(exts,strlen(exts)-3,3)
endelse  


;define some file locations
; -------------------------
fpdir      = groupdir				;directory to store plots in
fdata      = fpdir+'/groups/groups_'+exts	;base name of group catalogue
fsnap      = fpdir+'/snapshot_'+exts		;base name of snapshots
run_label  = 'Sbc'				;a name for the run
snap_label = 'snapshot_'+exts			;a name for a given snapshot



;Read in the group data
;test =read_groups(fdata,ngroups,nsubgroups,npart,fofcat,fofprop,subcat,subprop,ids,mass,pos,vel,types,rho,u,nel,nHI,hsml,sfr,age,z)
test= tdwarfs_read_fofgroup(fdata,ngroups,nsubgroups,npart,fofcat,fofprop,subcat,subprop,ids,mass,pos,vel,types,rho,u,nel,nHI,hsml,sfr,age,z)


; now parse the group information
; ---------------------------------
m_group       = fltarr(ngroups)
m_group_bound = fltarr(nsubgroups)

pos_group       = fltarr(ngroups,3)

for i=0,ngroups-1 do begin
	offset = fofcat.offset(i)
	length = fofcat.length(i)
	totmass= total(mass(offset:offset+length-1))
	totmass= totmass / 0.7
	m_group(i)= totmass
	;pos_group(i,0)= mean(pos(0,offset:offset+length-1))
	;pos_group(i,1)= mean(pos(1,offset:offset+length-1))
	;pos_group(i,2)= mean(pos(2,offset:offset+length-1))
	pos_group(i,0)= median(pos(0,offset:offset+length-1))
	pos_group(i,1)= median(pos(1,offset:offset+length-1))
	pos_group(i,2)= median(pos(2,offset:offset+length-1))
	print, "group i= ", i, "  len= ", length, " ", totmass, pos_group(i,0), pos_group(i,1), pos_group(i,2)

	gids= ids(offset:offset+length-1)

	;save the group ID's (if needed/wanted)
	idfile= '/home/tcox/Sbc10x/group_'+strcompress(string(i), /remove_all)+'_idlist.txt'
	openw, 1, idfile, ERROR=err
	printf, 1, "#  file: "+idfile
	printf, 1, "#  i: "+strcompress(string(i), /remove_all)
	printf, 1, "#  N: "+strcompress(string(length), /remove_all)
	printf, 1, "#  Mass: "+strcompress(string(totmass), /remove_all)
	printf, 1, "#  x|y|z = ", pos_group(i,0), pos_group(i,1), pos_group(i,2)
	for id=0L, length-1 do printf, 1, gids[id]
	close, 1
endfor
for i=0,nsubgroups-1 do begin
	offsetb= subcat.offset(i)
	lengthb= subcat.length(i)
	totmass = totmass
	m_group_bound(i) = totmass
endfor



; ditch the primary galaxy
;--------------------------
;print, "WARNING: ditching masses above 10^9"
;idx= where(m_group lt 1.0e+10)
;if idx(0) ne -1 then m_group= m_group(idx)



;make the mass function
;------------------------

; compute histogram
m_bin= (mass_max - mass_min)/nbins
; all groups
;mass_function = histogram(alog10(m_group(0:ngroups-1)),binsize=m_bin,min=mass_min,max=mass_max)
; ditch the biggest group
mass_function = histogram(alog10(m_group(1:ngroups-1)),binsize=m_bin,min=mass_min,max=mass_max)

; make sure it extends to 0
mass_function= [mass_function(0), mass_function, mass_function(nbins-1)]

print, "largest group= ", max(mass_function)

; define x-axis 
mass_function_x= (mass_max-mass_min)*findgen(n_elements(mass_function))/float(n_elements(mass_function)-1) + mass_min
mass_function_x= [mass_min, mass_function_x, mass_max]


printthis= mass_function
if keyword_set(normalized) then printthis= mass_function/total(mass_function)



oplot, mass_function_x, printthis, psym=10, thick=5.0, color=linecolor


; fill in histogram
; ------------------
xbins= mass_function_x+(m_bin*0.5)          ; make x coord
xbins[0]= mass_min
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




; fit a power law to the mass function
; ------------------------------------
if keyword_set(fit_powerlaw) then begin
	si = where(mass_function_x gt 5.25 and mass_function gt 0,ns)
	if ns gt 0 then begin
	        s = linfit(mass_function_x(si),alog10(mass_function(si)))
	        print,s

	        ;write the best fit slope to the file
	        slabel=strcompress(string(s(1)),/remove_all) ; sets exts to snap number
	        slabel=strmid(slabel,0,5)
	        ;tlabel=tlabel+", !Mb = "+slabel
	        ;xyouts,xxo,0.45,tlabel,/data,font=1,charsize=1.3,color=c0
	        xyouts,0.20,0.85,slabel,/normal,font=1,charsize=1.3,color=c0

	endif
endif


end



; -------------------------------------------------------------------------------


