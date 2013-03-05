;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Histogram of "Earth's"  position during merger with M31
;   -----------------------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro ehist, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "ehist, junk"
   print, "  "
   return
endif

;filename='lgehist.eps'
filename='lgehist2.eps'
snapnum=350
;snapnum=150


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; radius
xaxistitle= "!6Radius (kpc)"
xmax = 200.0
xmin = 0.0

; number (histogram)
yaxistitle= ' '
ymax = 1.2
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.05
y0= 0.15
x1= 0.98
y1= 0.95

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
	color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


; -----------------------------------------------


; ------------------
; compute histogram
; ------------------



frun= "/raid4/tcox/localgroup/bhires"
mw_bhid= 817955
m31_bhid= 2145299

ok= fload_snapshot_bh(frun,snapnum)

mw_bh_xyz= fload_blackhole_xyz('xyz',idtofollow=mw_bhid,center=[0,0,0])
print, "MW blackhole ID= ", mw_bhid
print, " position = ", mw_bh_xyz


orig_center= mw_bh_xyz



;
; load id's from list
; -----------------------
;idlist_fmfile= fload_id_list('/raid4/tcox/localgroup/v4/earth_idlist.txt')
idlist_fmfile= fload_id_list('/raid4/tcox/localgroup/bhires/earth_idlist.txt')
;frun= fload_frun(1)
;idlist_fmfile= fload_id_list(frun+'/id_list_test.txt')
print, "id's in file= ",n_elements(idlist_fmfile)
sidx=sort(idlist_fmfile)
idlist_fmfile= idlist_fmfile(sidx)



; allstars
; ----------
;x= fload_allstars_xyz('x',center=[0,0,0])
;y= fload_allstars_xyz('y',center=[0,0,0])
;z= fload_allstars_xyz('z',center=[0,0,0])
x= fload_allstars_xyz('x',center=orig_center)
y= fload_allstars_xyz('y',center=orig_center)
z= fload_allstars_xyz('z',center=orig_center)
gid= fload_allstars_ids(1)



; grab appropriate particles
; -----------------------------
idx= intarr(n_elements(gid))
lstid= -1
duplicates= 0
for i=0,n_elements(idlist_fmfile)-1 do begin
    if idlist_fmfile[i] eq lstid then begin
	duplicates= duplicates+1
	print, lstid, idlist_fmfile[i]
    endif
    inlst= where(gid eq idlist_fmfile[i])
    if inlst(0) ne -1 then begin 
	idx(inlst)= 1 
    endif else begin
	;print, "couldn't find id=",idlist_fmfile[i]
    endelse
    lstid= idlist_fmfile[i]
endfor
midx= where(idx eq 1)
print, "matching id's= ",n_elements(midx)
print, "duplicate id's= ",duplicates




if midx(0) ne -1 then begin
	ex= x(midx)
	ey= y(midx)
	ez= z(midx)
endif else begin
	print, " "
	print, " "
	print, " PROBLEM: no matching ID particles"
	print, " "
	print, " "
	return
endelse



earthr= sqrt(ex*ex + ey*ey + ez*ez)

print, "earth radius max/min= ", max(earthr), min(earthr)


temp= process_histogram(earthr, xmax=xmax, xmin=xmin, levels=100, oplotit=50)
print, min(temp), max(temp)

; -----------------------------------------------

;avg_metallicity= mean(newstar_metals)
;print, "avg= ", avg_metallicity
;avglbl= strcompress(string(avg_metallicity),/remove_all)
;xyouts, 0.10, 0.75, "<Z(Z!D!9n!6!N)> = "+avglbl, size=1.2, color=0, /normal, charthick=3.0



;xyouts, 0.75, 0.90, fload_timelbl(1,2), size=1.2, color=0, /normal, charthick=4.0
xyouts, 0.75, 0.90, fload_timelbl(0.7,2), size=1.2, color=0, /normal, charthick=4.0
;xyouts, 0.55, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0
;xyouts, 0.75, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0

;xyouts, 0.10, 0.85, "!6V!D200!N = 50 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0
;xyouts, 0.10, 0.85, "!6V!D200!N = 500 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0
;xyouts, 0.10, 0.85, "!6V!D200!N = 160 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0



; print extras
; -------------


; done
; -----
device, /close


end






;====================================================================================






pro ehist_comp, msg, snap, filename=filename

if not keyword_set(msg) then begin
   print, "  "
   print, "ehist_comp, msg, snap(array), filename=filename"
   print, "  "
   return
endif

;filename='lgehist.eps'
if not keyword_set(filename) then filename='lgehistall.eps'
if not keyword_set(msg) then msg= 'post-merger'
if not keyword_set(snap) then snap= [317, 100, 120, 120, 120]


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; radius
xaxistitle= "!6Radius (kpc)"
xmax = 100.0
xmin = 0.0

; number (histogram)
yaxistitle= ' '
ymax = 1.2
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.05
y0= 0.15
x1= 0.95
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
	color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


; -----------------------------------------------


; ------------------
; compute histogram
; ------------------


;do_earth_hist_general, "/raid4/tcox/localgroup/bhires", 317, 817955L, 2145299L, $
do_earth_hist_general, "/raid4/tcox/localgroup/bhires", snap[0], 817955L, 2145299L, $
			'/raid4/tcox/localgroup/bhires/earth_idlist.txt.dr=0.1', plotcolor=150, $
                        hist_xmax=xmax, hist_xmin=xmin, $
                        orig_center=orig_center


;do_earth_hist_general, "/raid4/tcox/localgroup/v4", 100, 500001L, 1300002L, $
do_earth_hist_general, "/raid4/tcox/localgroup/v4", snap[1], 500001L, 1300002L, $
			'/raid4/tcox/localgroup/v4/earth_idlist.txt', plotcolor=50, $
                        hist_xmax=xmax, hist_xmin=xmin, $
                        orig_center=orig_center


;do_earth_hist_general, "/raid4/tcox/localgroup/v7", 120, 500001L, 1300002L, $
do_earth_hist_general, "/raid4/tcox/localgroup/v7", snap[2], 500001L, 1300002L, $
			'/raid4/tcox/localgroup/v7/earth_idlist.txt', plotcolor=0, $
                        hist_xmax=xmax, hist_xmin=xmin, $
                        orig_center=orig_center


;do_earth_hist_general, "/raid4/tcox/localgroup/v17", 120, 500001L, 1300002L, $
do_earth_hist_general, "/raid4/tcox/localgroup/v17", snap[3], 500001L, 1300002L, $
			'/raid4/tcox/localgroup/v17/earth_idlist.txt', plotcolor=100, $
                        hist_xmax=xmax, hist_xmin=xmin, $
                        orig_center=orig_center


;do_earth_hist_general, "/raid4/tcox/localgroup/v16", 120, 500001L, 1300002L, $
do_earth_hist_general, "/raid4/tcox/localgroup/v16", snap[4], 500001L, 1300002L, $
			'/raid4/tcox/localgroup/v16/earth_idlist.txt', plotcolor=200, $
                        hist_xmax=xmax, hist_xmin=xmin, $
                        orig_center=orig_center






; -------------
; print extras
; -------------

;xyouts, 0.55, 0.90, 'post-merger', size=1.8, color=0, /normal, charthick=4.0
xyouts, 0.48, 0.90, msg, size=1.8, color=0, /normal, charthick=4.0




; -----
; done
; -----

device, /close


end


















;========================================================================================
;
;
;   This produces the earth radius histograms for various interesting times.
;
;     event              time in v5   time in bhires
;                           (Gyr)     (Gyr), snapnum
;     -----------------   ------       ------------
;     today                  5.0         6.48, 150  
;     first passage          6.8         8.49, 197
;     post-first             7.2         9.0,  209
;     second passage         8.5         12.0, 279
;     post-second passage    9.0         12.5, 291
;     merger                 9.5         13.6, 317
;
;



pro do_6, junk

; careful
frun= "/raid4/tcox/localgroup/bhires"


ehist2, ' ', snapnum=150, filename=frun+'/earthist_150.eps', msg='Today'
ehist2, ' ', snapnum=154, filename=frun+'/earthist_154.eps', msg='+0.15 Gyr'
ehist2, ' ', snapnum=197, filename=frun+'/earthist_197.eps', msg='+1.8 Gyr, first passage'
ehist2, ' ', snapnum=209, filename=frun+'/earthist_209.eps', msg='+2.2 Gyr, post-first passage'
ehist2, ' ', snapnum=279, filename=frun+'/earthist_279.eps', msg='+3.5 Gyr, second passage'
ehist2, ' ', snapnum=291, filename=frun+'/earthist_291.eps', msg='+4.0 Gyr, post-second passage'
ehist2, ' ', snapnum=317, filename=frun+'/earthist_317.eps', msg='+4.5 Gyr, post-merger'

end





;========================================================================================
;
;
;   ----------------------
;   |                    |
;   |                    |
;   |   stellar image    |
;   |                    |
;   |                    |
;   |    + earth pts.    |
;   |                    |
;   ----------------------
;   |                    |
;   |                    |
;   |    radial hist     |
;   |                    |
;   |                    |
;   |                    |
;   |                    |
;   ----------------------
;
;
;





pro ehist2, junk, snapnum=snapnum, filename=filename, msg=msg

if not keyword_set(junk) then begin
   print, "  "
   print, "ehist2, junk"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='lgehist.eps'



; ------------------
; Plot this up
; ------------------
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=0, newxsize= 13.0, newysize= 23.0



x0= 0.02
x1= 0.98
xs= (x1-x0)

y0= 0.10
y1= 0.45
y2= 0.99



; -----------------------------
; set up constants
; -----------------------------


; radius
xaxistitle= "!6Radius (kpc)"
;hist_xmax = 200.0
hist_xmax = 100.0
;hist_xmax = 50.0
hist_xmin = 0.0

; number (histogram)
yaxistitle= ' '
hist_ymax = 1.2
hist_ymin = 0.0


frun= "/raid4/tcox/localgroup/bhires"
mw_bhid= 817955
m31_bhid= 2145299


earthidlist= '/raid4/tcox/localgroup/bhires/earth_idlist.txt'
; also set manually in lg_contour_add_idlist


;snapnum=350
;snapnum=209
if not keyword_set(snapnum) then snapnum=209


xlen= 100.0


; ----------------------------------------
;   1111:   Stellar Image + Earth points
; ----------------------------------------


ok= fload_snapshot_bh(frun,snapnum)

;do_allstars= 0
do_allstars= 1

if do_allstars eq 1 then begin
    ndisk= fload_npart(2)                               ; disk particles
    nbulge= fload_npart(3)                              ; bulge particles
    nstars= fload_npart(4)
    npart= long(ndisk) + long(nbulge) + long(nstars)
    print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars

    x= fload_allstars_xyz('x',center=[0,0,0])
    y= fload_allstars_xyz('y',center=[0,0,0])
    z= fload_allstars_xyz('z',center=[0,0,0])
    m= fload_allstars_mass(1)
    ;hsml= fload_allstars_hsml(1)
    hsml= 0.2 + m*0.0


    ;x= x - center[0]
    ;y= y - center[1]
    ;z= z - center[2]

endif

	;center= [0.0, 0.0, 0.0]
	center= fload_blackhole_xyz('xyz',idtofollow=mw_bhid,center=[0,0,0])
	orig_center= center
	print, "MW blackhole ID= ", mw_bhid
	print, " using center = ", center

	

	;set_maxden= 10.0
	set_maxden= 1.0    ; 10^4 msolar/pc2
	;set_dynrng= 1.0e+5
	;set_dynrng= 1.0e+6
	;set_maxden= 0.5e-1
	set_dynrng= 1.0e+4

	;pixels= 256L


	contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
	        hsml=hsml, $
	        filename=filename, fitstoo=fitstoo, $
	        xthickness=xthickness, ythickness=ythickness, $
	        pixels=pixels, zthickness=zthickness, $
	        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	        crude=crude, center=center, $
	        set_maxden= set_maxden, $
	        set_dynrng= set_dynrng, $
	        NxNImage=NxNImage


	tv, NxNImage, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal


        !p.position=[x0,y1,x1,y2]
        !p.ticklen=0.03

        xmin= -xlen
        xmax=  xlen
        ymin= -xlen
        ymax=  xlen

        plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
                      xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'



                ;drawcenters= 1
                if keyword_set(drawcenters) then begin

                        ; warning, need to load process_directory for this
                        read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=fload_frun(1)+"/centers.txt"
                        cp_cen2_x= cp_cen2[0,*]
                        cp_cen2_y= cp_cen2[1,*]
                        cp_cen2_z= cp_cen2[2,*]
                        idx= where(cp_time le fload_time(1))
                        cp_cen2_x= cp_cen2_x(idx)
                        cp_cen2_y= cp_cen2_y(idx)

                        oplot, cp_cen2_x, cp_cen2_y, psym=-3, linestyle=0, color= 150, thick=4.0


                        ; double up baby
                        ;read_centerpositions, cp_time, cp_cen1, cp_cen2, filename="Data/nodrag/centers.txt"
                        ;cp_cen2_x= cp_cen2[0,*]
                        ;cp_cen2_y= cp_cen2[1,*]
                        ;idx= where(cp_time le fload_time(1))
                        ;cp_cen2_x= cp_cen2_x(idx)
                        ;cp_cen2_y= cp_cen2_y(idx)
                        ;oplot, cp_cen2_x, cp_cen2_y, psym=-3, linestyle=0, color= 100, thick=4.0


                endif




                ; -----------------------
                ;
                ;   Overplot ID's
                ;
                ; -----------------------
                ;overplot_ids= 0
                overplot_ids= 1
                if overplot_ids eq 1 then begin
                        ;fancy_center= 1
                        ;if fancy_center eq 1 then center= fload_center_alreadycomp(1)
			center= orig_center
                        lg_contour_add_idlist, 'junk', center=center
                        ;fload_newcolortable, 1
                        fload_newcolortable, 4
                endif







; ----------------------------------------
;   2222:   Histogram of Earth Radii
; ----------------------------------------
!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[hist_xmin,hist_xmax], yrange=[hist_ymin,hist_ymax], xstyle=1, ystyle=1, $
	color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


; -----------------------------------------------


; ------------------
; compute histogram
; ------------------



mw_bh_xyz= fload_blackhole_xyz('xyz',idtofollow=mw_bhid,center=[0,0,0])
print, "MW blackhole ID= ", mw_bhid
print, " position = ", mw_bh_xyz


orig_center= mw_bh_xyz



;do_earth_hist, earthidlist, plotcolor=50, $
;			hist_xmax=hist_xmax, hist_xmin=hist_xmin, $
;			orig_center=orig_center


;do_earth_hist, '/raid4/tcox/localgroup/bhires/earth_idlist.txt.dr=0.01', plotcolor=50, $
;                       hist_xmax=hist_xmax, hist_xmin=hist_xmin, $
;                       orig_center=orig_center

;do_earth_hist, '/raid4/tcox/localgroup/bhires/earth_idlist.txt.dr=0.05', plotcolor=150, $
;                       hist_xmax=hist_xmax, hist_xmin=hist_xmin, $
;                       orig_center=orig_center

do_earth_hist, '/raid4/tcox/localgroup/bhires/earth_idlist.txt.dr=0.1', plotcolor=150, $
                       hist_xmax=hist_xmax, hist_xmin=hist_xmin, $
                       orig_center=orig_center

;do_earth_hist, '/raid4/tcox/localgroup/bhires/earth_idlist.txt.dr=0.2', plotcolor=100, $
;                       hist_xmax=hist_xmax, hist_xmin=hist_xmin, $
;                       orig_center=orig_center

;do_earth_hist, '/raid4/tcox/localgroup/bhires/earth_idlist.txt.dr=0.5', plotcolor=200, $
;                       hist_xmax=hist_xmax, hist_xmin=hist_xmin, $
;                       orig_center=orig_center



; -----------------------------------------------
; print extras
;

;avg_metallicity= mean(newstar_metals)
;print, "avg= ", avg_metallicity
;avglbl= strcompress(string(avg_metallicity),/remove_all)
;xyouts, 0.10, 0.75, "<Z(Z!D!9n!6!N)> = "+avglbl, size=1.2, color=0, /normal, charthick=3.0



;xyouts, 0.75, 0.90, fload_timelbl(1,2), size=1.2, color=0, /normal, charthick=4.0
;xyouts, 0.75, 0.90, fload_timelbl(0.7,2), size=1.2, color=0, /normal, charthick=4.0
;xyouts, 0.08, 0.94, fload_timelbl(0.7,2,/noteq), size=1.5, color=0, /normal, charthick=4.0
;xyouts, 0.55, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0
;xyouts, 0.75, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0

;xyouts, 0.10, 0.85, "!6V!D200!N = 50 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0
;xyouts, 0.10, 0.85, "!6V!D200!N = 500 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0
;xyouts, 0.10, 0.85, "!6V!D200!N = 160 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0


xyouts, 0.08, 0.94, msg, size=2.0, color=0, /normal, charthick=2.0

; ----------------------
; done
;
device, /close


end






;====================================================================================





pro do_earth_hist, earthidlist, plotcolor=plotcolor, $
			hist_xmax=hist_xmax, hist_xmin=hist_xmin, $
			orig_center=orig_center

;earthidlist= '/raid4/tcox/localgroup/bhires/earth_idlist.txt'
; also set manually in lg_contour_add_idlist

;
; load id's from list
; -----------------------
idlist_fmfile= fload_id_list(earthidlist)
print, "id's in file= ",n_elements(idlist_fmfile)
sidx=sort(idlist_fmfile)
idlist_fmfile= idlist_fmfile(sidx)



; allstars
; ----------
;x= fload_allstars_xyz('x',center=[0,0,0])
;y= fload_allstars_xyz('y',center=[0,0,0])
;z= fload_allstars_xyz('z',center=[0,0,0])
x= fload_allstars_xyz('x',center=orig_center)
y= fload_allstars_xyz('y',center=orig_center)
z= fload_allstars_xyz('z',center=orig_center)
gid= fload_allstars_ids(1)



; grab appropriate particles
; -----------------------------
idx= intarr(n_elements(gid))
lstid= -1
duplicates= 0
for i=0,n_elements(idlist_fmfile)-1 do begin
	; -----
	if (i mod 100) eq 0 then print, i
	; -----
    if idlist_fmfile[i] eq lstid then begin
	duplicates= duplicates+1
	print, lstid, idlist_fmfile[i]
    endif
    inlst= where(gid eq idlist_fmfile[i])
    if inlst(0) ne -1 then begin 
	idx(inlst)= 1 
    endif else begin
	;print, "couldn't find id=",idlist_fmfile[i]
    endelse
    lstid= idlist_fmfile[i]
endfor
midx= where(idx eq 1)
print, "matching id's= ",n_elements(midx)
print, "duplicate id's= ",duplicates




if midx(0) ne -1 then begin
	ex= x(midx)
	ey= y(midx)
	ez= z(midx)
endif else begin
	print, " "
	print, " "
	print, " PROBLEM: no matching ID particles"
	print, " "
	print, " "
	return
endelse



earthr= sqrt(ex*ex + ey*ey + ez*ez)

print, "xxxxxxxxxxxxxxxx"
print, "earth radius max/min= ", max(earthr), min(earthr)
print, "n= ", n_elements(earthr)
print, "frac > 10 kpc= ", 1.0*n_elements(where(earthr gt 10.0))/n_elements(earthr)
print, "frac > 20 kpc= ", 1.0*n_elements(where(earthr gt 20.0))/n_elements(earthr)
print, "frac > 30 kpc= ", 1.0*n_elements(where(earthr gt 30.0))/n_elements(earthr)


m31_bhid= 2145299
m31_bh_xyz= fload_blackhole_xyz('xyz',idtofollow=m31_bhid,center=orig_center)
print, "M31 blackhole ID= ", m31_bhid
print, " position = ", m31_bh_xyz

mx= ex - m31_bh_xyz[0]
my= ey - m31_bh_xyz[1]
mz= ez - m31_bh_xyz[2]
m31_r= sqrt(mx*mx + my*my + mz*mz)

print, "n_m31_r= ", n_elements(m31_r)
print, "frac < 10 of m31= ", 1.0*n_elements(where(m31_r lt 10.0))/n_elements(m31_r)


temp= process_histogram(earthr, xmax=hist_xmax, xmin=hist_xmin, levels=100, oplotit=plotcolor)
print, min(temp), max(temp)

; -----------------------------------------------

end






;====================================================================================





pro do_earth_hist_general, frun, snapnum, mw_bhid, m31_bhid, $
			earthidlist, plotcolor=plotcolor, $
			hist_xmax=hist_xmax, hist_xmin=hist_xmin, $
			orig_center=orig_center


ok= fload_snapshot_bh(frun,snapnum)

bhid= fload_blackhole_id(1)
if n_elements(bhid) eq 1 then mw_bhid= bhid

mw_bh_xyz= fload_blackhole_xyz('xyz',idtofollow=mw_bhid,center=[0,0,0])
print, "MW blackhole ID= ", mw_bhid
print, " position = ", mw_bh_xyz


orig_center= mw_bh_xyz


;
; load id's from list
; -----------------------
idlist_fmfile= fload_id_list(earthidlist)
print, "id's in file= ",n_elements(idlist_fmfile)
sidx=sort(idlist_fmfile)
idlist_fmfile= idlist_fmfile(sidx)



; allstars
; ----------
;x= fload_allstars_xyz('x',center=[0,0,0])
;y= fload_allstars_xyz('y',center=[0,0,0])
;z= fload_allstars_xyz('z',center=[0,0,0])
x= fload_allstars_xyz('x',center=orig_center)
y= fload_allstars_xyz('y',center=orig_center)
z= fload_allstars_xyz('z',center=orig_center)
gid= fload_allstars_ids(1)



; grab appropriate particles
; -----------------------------
idx= intarr(n_elements(gid))
lstid= -1
duplicates= 0
for i=0,n_elements(idlist_fmfile)-1 do begin
	; -----
	if (i mod 100) eq 0 then print, i
	; -----
    if idlist_fmfile[i] eq lstid then begin
	duplicates= duplicates+1
	print, lstid, idlist_fmfile[i]
    endif
    inlst= where(gid eq idlist_fmfile[i])
    if inlst(0) ne -1 then begin 
	idx(inlst)= 1 
    endif else begin
	;print, "couldn't find id=",idlist_fmfile[i]
    endelse
    lstid= idlist_fmfile[i]
endfor
midx= where(idx eq 1)
print, "matching id's= ",n_elements(midx)
print, "duplicate id's= ",duplicates




if midx(0) ne -1 then begin
	ex= x(midx)
	ey= y(midx)
	ez= z(midx)
endif else begin
	print, " "
	print, " "
	print, " PROBLEM: no matching ID particles"
	print, " "
	print, " "
	return
endelse



earthr= sqrt(ex*ex + ey*ey + ez*ez)

print, "xxxxxxxxxxxxxxxx"
print, "earth radius max/min= ", max(earthr), min(earthr)
print, "n= ", n_elements(earthr)
print, "frac > 10 kpc= ", 1.0*n_elements(where(earthr gt 10.0))/n_elements(earthr)
print, "frac > 20 kpc= ", 1.0*n_elements(where(earthr gt 20.0))/n_elements(earthr)
print, "frac > 30 kpc= ", 1.0*n_elements(where(earthr gt 30.0))/n_elements(earthr)


m31_bhid= 2145299
m31_bh_xyz= fload_blackhole_xyz('xyz',idtofollow=m31_bhid,center=orig_center)
print, "M31 blackhole ID= ", m31_bhid
print, " position = ", m31_bh_xyz

mx= ex - m31_bh_xyz[0]
my= ey - m31_bh_xyz[1]
mz= ez - m31_bh_xyz[2]
m31_r= sqrt(mx*mx + my*my + mz*mz)

print, "n_m31_r= ", n_elements(m31_r)
print, "frac < 10 of m31= ", 1.0*n_elements(where(m31_r lt 10.0))/n_elements(m31_r)


temp= process_histogram(earthr, xmax=hist_xmax, xmin=hist_xmin, levels=100, oplotit=plotcolor)
print, min(temp), max(temp)


end






;====================================================================================





