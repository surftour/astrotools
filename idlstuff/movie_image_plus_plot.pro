;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;  
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro do_seq, frun

   starti= 0
   ;starti= 1
   ;starti= 35

   ;endi= 11

   spawn, "/bin/ls "+frun+"/snap* | wc ",result
   endi=long(result[0])-1

   print, "frun= ", frun
   print, "starti= ", starti
   print, "endi= ", endi

   spawn, "mkdir "+frun+"/movie_plot"

   for i=starti,endi do begin

	thisi= i

	exts='0000'+strcompress(string(thisi),/remove_all)

	;thisfile= frun+'/movie_plot/img_sfr_'+strmid(exts,strlen(exts)-3,3)+'.eps'
	;thisfile= frun+'/movie_plot/img_sfr_'+strmid(exts,strlen(exts)-3,3)+'.jpg'
	;img2_sfr, frun, thisi, filename=thisfile
	;img2_sfr, frun, thisi, filename='temp.eps'
	;img4_sfr, frun, thisi, filename='temp.eps'

	;thisfile= frun+'/movie_plot/img_bhacc_'+strmid(exts,strlen(exts)-3,3)+'.jpg'
	;img2_bhacc, frun, thisi, filename='temp2.eps'

	thisfile= frun+'/movie_plot/img_both_'+strmid(exts,strlen(exts)-3,3)+'.jpg'
	img4_both, frun, thisi, filename='temp.eps'

	cmd= "convert -quality 95 temp.eps "+thisfile
	spawn, cmd

   endfor

end












;==================================================================================
;
;
;   ----------------
;   |              |
;   |              |
;   |              |    ---------------------
;   |              |    |                   |
;   |              |    |                   |
;   |              |    |                   |
;   ----------------    |      SFR          |
;   |              |    |                   |
;   |              |    |                   |
;   |              |    |                   |
;   |              |    ---------------------
;   |              |
;   |              |
;   ----------------
; 
;
;
;
;
;

pro img2_sfr, frun, snapnum, filename=filename
        

if not keyword_set(frun) then begin
   print, "  "
   print, " img2_sfr, frun, snapnum, filename=filename, "
   print, "  "
   print, "  "
   print, "  "
   return
endif   
        
if not keyword_set(snapnum) then snapnum=12
if not keyword_set(filename) then filename="img_sfr.eps"




;--------------------------------------------

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=30, newysize=22

;--------------------------------------------



; now images
;-----------------

xlen= 50.0

x0= 0.03
x1= 0.52

y0= 0.03
y1= 0.50
y2= 0.97



ok=fload_snapshot_bh(frun,snapnum)

; option 1 - use computed center
;center=fload_center_alreadycomp(1)
;orig_center= center
;print, "alreadycomp center= ", center

; option 2 - use BH center
;bhid= fload_blackhole_id(1)
;bhid= bhid[0]
;bhid= bhid[1]
;bhid= 200001L    ; used for ds/vc3vc3h_2 
;bhid= 280002L   ; used for z3/b4e
;bhid= 400002L   ; used for ds/vc3vc3e_2
;center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)

;print, "Blackhole ID: ", bhid
;print, "Blackhole center: ", center_bh
;orig_center= center_bh


orig_center= [0,0,0]


;do_gas= 0
do_gas= 1

if do_gas eq 1 then begin
    npart= fload_npart(0)
    N= long(npart)

    center= orig_center
    x= fload_gas_xyz('x',center=center) / 0.7
    y= fload_gas_xyz('y',center=center) / 0.7
    z= fload_gas_xyz('z',center=center) / 0.7
    m= fload_gas_mass(1) / 0.7
    hsml= fload_gas_hsml(1) / 0.7

endif



; ----------------------------------------------
;  Call our generic contour plotting procedure
; ----------------------------------------------

;set_maxden= 10.0  
set_maxden= 1.0    ; 10^4 msolar/pc2
;set_dynrng= 1.0e+3
;set_dynrng= 1.0e+4
set_dynrng= 1.0e+5
;set_dynrng= 1.0e+6
;set_maxden= 0.5e-1
;set_dynrng= 1.0e+6


center= [0.0, 0.0, 0.0]

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
        hsml=hsml, $
        xthickness=xthickness, ythickness=ythickness, $
        pixels=pixels, zthickness=zthickness, $
        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        crude=crude, center=center, $
        set_maxden= set_maxden, $
        set_dynrng= set_dynrng, $
        NxNImage=NxNImage

        XYImage= NxNImage


xz= 1
center= [0.0, 0.0, 0.0]

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
        hsml=hsml, $
        xthickness=xthickness, ythickness=ythickness, $
        pixels=pixels, zthickness=zthickness, $
        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        crude=crude, center=center, $
        set_maxden= set_maxden, $
        set_dynrng= set_dynrng, $
        NxNImage=NxNImage

        XZImage= NxNImage



tv, XYImage, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal
tv, XZImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


xmin= -xlen
xmax=  xlen
ymin= -xlen
ymax=  xlen

; 
; xy
;
;------
!p.position=[x0,y1,x1,y2]
!p.ticklen=0.03

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'

xyouts, x0+0.02, y2-0.04, '!6xy', size=1.2, /normal, color= 0, charthick=1.5

n_bh= fload_npart(5)
bhid= fload_blackhole_id(1)

for i=0, n_bh-1 do begin
	;usersym,2.0*cos(findgen(49)/49*2*!pi),2.0*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[1]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[1]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor


; 
; xy
;
;------
!p.position=[x0,y0,x1,y1]

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'


xyouts, x0+0.02, y1-0.04, '!6xz', size=1.2, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[2]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[2]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor





;==================================================================================






; and SFR
;-----------------

; needs .run sfr_multi

x0= 0.65
x1= 0.95
y0= 0.20
y1= 0.80

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

xmin= 0.0
xmax= 3.0
ymax= 300.0
ymin= 0.05


!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, /ylog, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        ytickformat='exp_label', xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase, /normal



open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac

; physical units
sfrtime = sfrtime / 0.7

oplot, sfrtime, sfrsfr, psym=-3, color= 0, linestyle= 0, thick= 4.0

idx=where(sfrtime ge (fload_time(1)/0.7))
if idx(0) eq -1 then idx[0]= n_elements(sfrtime)-1
oplot, [sfrtime[idx[0]]], [sfrsfr[idx[0]]], psym= 4, color= 150, thick=10.0, symsize=1.6


sfrtot= total(fload_gas_sfr(1))
lblz = strcompress(string(sfrtot),/remove_all)
digz= 5
if sfrtot lt 100.0 then digz= 4
if sfrtot lt 10.0 then digz= 3
lblz = strmid(lblz,0,digz)        ; 0.xxx
xyouts, x0+0.04, y0+0.08, lblz+" M!D!9n!6!N Yr!E-1!N", /normal, color= 0, charthick=1.5, size=1.5




;------------------
device, /close

end








;==================================================================================
;
;
;   -------------------------------
;   |              |              |
;   |              |              |
;   |              |              |    ---------------------
;   |              |              |    |                   |
;   |              |              |    |                   |
;   |              |              |    |                   |
;   -------------------------------    |      SFR          |
;   |              |              |    |                   |
;   |              |              |    |                   |
;   |              |              |    |                   |
;   |              |              |    ---------------------
;   |              |              |
;   |              |              |
;   -------------------------------
; 
;
;
;
;
;

pro img4_sfr, frun, snapnum, filename=filename
        

if not keyword_set(frun) then begin
   print, "  "
   print, " img4_sfr, frun, snapnum, filename=filename, "
   print, "  "
   print, "  "
   print, "  "
   return
endif   
        
if not keyword_set(snapnum) then snapnum=12
if not keyword_set(filename) then filename="img_sfr.eps"




;--------------------------------------------

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=36, newysize=24

;--------------------------------------------



; now images
;-----------------

xlen= 50.0

x0= 0.01
x1= 0.33
x2= 0.65

y0= 0.015
y1= 0.50
y2= 0.985



ok=fload_snapshot_bh(frun,snapnum)

n_bh= fload_npart(5)
bhid= fload_blackhole_id(1)

orig_center= [0,0,0]



;==================================================================================


fload_newcolortable, 1


;
;   Gas
;
; ----------------------------------------------

npart= fload_npart(0)
N= long(npart)

center= orig_center
x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
m= fload_gas_mass(1) / 0.7
hsml= fload_gas_hsml(1) / 0.7



;set_maxden= 10.0  
set_maxden= 1.0    ; 10^4 msolar/pc2
;set_maxden= 0.5e-1

;set_dynrng= 1.0e+3
;set_dynrng= 1.0e+4
set_dynrng= 1.0e+5
;set_dynrng= 1.0e+6
;set_dynrng= 1.0e+6


center= orig_center

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, hsml=hsml, center=center, $
        pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        set_maxden= set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

        XYImage= NxNImage


xz= 1
center= orig_center

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, hsml=hsml, center=center, $
        pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        set_maxden= set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

        XZImage= NxNImage



tv, XYImage, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal
tv, XZImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


xmin= -xlen
xmax=  xlen
ymin= -xlen
ymax=  xlen


; 
; xy
;
;---------------------
!p.position=[x0,y1,x1,y2]
!p.ticklen=0.03

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'

xyouts, x0+0.02, y2-0.04, '!6gas - xy', size=1.5, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	;usersym,2.0*cos(findgen(49)/49*2*!pi),2.0*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[1]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[1]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor


; 
; xz
;
;---------------------
!p.position=[x0,y0,x1,y1]

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'


xyouts, x0+0.02, y1-0.04, '!6gas - xz', size=1.5, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[2]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[2]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor





;==================================================================================


x0= x1 & x1= x2

fload_newcolortable, 4

;
;   Stars
;
; ----------------------------------------------

npart= fload_npart(2)+fload_npart(3)+fload_npart(4)
N= long(npart)

center= orig_center
x= fload_allstars_xyz('x',center=center) / 0.7
y= fload_allstars_xyz('y',center=center) / 0.7
z= fload_allstars_xyz('z',center=center) / 0.7
m= fload_allstars_mass(1) / 0.7
;hsml= fload_allstars_hsml(1) / 0.7
hsml= m*0.0 + 0.2



center= orig_center

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, hsml=hsml, center=center, $
        pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        set_maxden= set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

        XYImage= NxNImage


xz= 1
center= orig_center

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, hsml=hsml, center=center, $
        pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        set_maxden= set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

        XZImage= NxNImage



tv, XYImage, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal
tv, XZImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


xmin= -xlen
xmax=  xlen
ymin= -xlen
ymax=  xlen


; 
; xy
;
;---------------------
!p.position=[x0,y1,x1,y2]
!p.ticklen=0.03

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'

xyouts, x0+0.02, y2-0.04, '!6stars - xy', size=1.5, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	;usersym,2.0*cos(findgen(49)/49*2*!pi),2.0*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[1]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[1]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor


; 
; xz
;
;---------------------
!p.position=[x0,y0,x1,y1]

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'


xyouts, x0+0.02, y1-0.04, '!6stars - xz', size=1.5, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[2]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[2]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor







;==================================================================================






; and SFR
;-----------------

; needs .run sfr_multi

x0= 0.75
x1= 0.98
y0= 0.30
y1= 0.70

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

xmin= 0.0
xmax= 3.0
ymax= 300.0
ymin= 0.05


!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, /ylog, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        ytickformat='exp_label', xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase, /normal



open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac

; physical units
sfrtime = sfrtime / 0.7

oplot, sfrtime, sfrsfr, psym=-3, color= 0, linestyle= 0, thick= 4.0

idx=where(sfrtime ge (fload_time(1)/0.7))
if idx(0) eq -1 then idx[0]= n_elements(sfrtime)-1
oplot, [sfrtime[idx[0]]], [sfrsfr[idx[0]]], psym= 4, color= 150, thick=5.0, symsize=1.6


sfrtot= total(fload_gas_sfr(1))
lblz = strcompress(string(sfrtot),/remove_all)
digz= 5
if sfrtot lt 100.0 then digz= 4
if sfrtot lt 10.0 then digz= 3
lblz = strmid(lblz,0,digz)        ; 0.xxx
xyouts, x0+0.04, y0+0.08, lblz+" M!D!9n!6!N Yr!E-1!N", /normal, color= 0, charthick=1.5, size=1.5




;------------------
device, /close

end








;==================================================================================
;
;
;   ----------------
;   |              |
;   |              |
;   |              |    ---------------------
;   |              |    |                   |
;   |              |    |                   |
;   |              |    |                   |
;   ----------------    |      BH Acc       |
;   |              |    |                   |
;   |              |    |                   |
;   |              |    |                   |
;   |              |    ---------------------
;   |              |
;   |              |
;   ----------------
; 
;
;
;
;
;

pro img2_bhacc, frun, snapnum, filename=filename
        

if not keyword_set(frun) then begin
   print, "  "
   print, " img2_bhacc, frun, snapnum, filename=filename, "
   print, "  "
   print, "  "
   print, "  "
   return
endif   
        
if not keyword_set(snapnum) then snapnum=12
if not keyword_set(filename) then filename="img_bhacc.eps"




;--------------------------------------------

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=30, newysize=22

;--------------------------------------------



; now images
;-----------------

xlen= 50.0

x0= 0.03
x1= 0.52

y0= 0.03
y1= 0.50
y2= 0.97



ok=fload_snapshot_bh(frun,snapnum)

; option 1 - use computed center
;center=fload_center_alreadycomp(1)
;orig_center= center
;print, "alreadycomp center= ", center

; option 2 - use BH center
;bhid= fload_blackhole_id(1)
;bhid= bhid[0]
;bhid= bhid[1]
;bhid= 200001L    ; used for ds/vc3vc3h_2 
;bhid= 280002L   ; used for z3/b4e
;bhid= 400002L   ; used for ds/vc3vc3e_2
;center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)

;print, "Blackhole ID: ", bhid
;print, "Blackhole center: ", center_bh
;orig_center= center_bh


orig_center= [0,0,0]


;do_gas= 0
do_gas= 1

if do_gas eq 1 then begin
    npart= fload_npart(0)
    N= long(npart)

    center= orig_center
    x= fload_gas_xyz('x',center=center) / 0.7
    y= fload_gas_xyz('y',center=center) / 0.7
    z= fload_gas_xyz('z',center=center) / 0.7
    m= fload_gas_mass(1) / 0.7
    hsml= fload_gas_hsml(1) / 0.7

endif



; ----------------------------------------------
;  Call our generic contour plotting procedure
; ----------------------------------------------

;set_maxden= 10.0  
set_maxden= 1.0    ; 10^4 msolar/pc2
;set_dynrng= 1.0e+3
;set_dynrng= 1.0e+4
set_dynrng= 1.0e+5
;set_dynrng= 1.0e+6
;set_maxden= 0.5e-1
;set_dynrng= 1.0e+6


center= [0.0, 0.0, 0.0]

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
        hsml=hsml, $
        xthickness=xthickness, ythickness=ythickness, $
        pixels=pixels, zthickness=zthickness, $
        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        crude=crude, center=center, $
        set_maxden= set_maxden, $
        set_dynrng= set_dynrng, $
        NxNImage=NxNImage

        XYImage= NxNImage


xz= 1
center= [0.0, 0.0, 0.0]

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
        hsml=hsml, $
        xthickness=xthickness, ythickness=ythickness, $
        pixels=pixels, zthickness=zthickness, $
        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        crude=crude, center=center, $
        set_maxden= set_maxden, $
        set_dynrng= set_dynrng, $
        NxNImage=NxNImage

        XZImage= NxNImage



tv, XYImage, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal
tv, XZImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


xmin= -xlen
xmax=  xlen
ymin= -xlen
ymax=  xlen

; 
; xy
;
;------
!p.position=[x0,y1,x1,y2]
!p.ticklen=0.03

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'

xyouts, x0+0.02, y2-0.04, '!6xy', size=1.2, /normal, color= 0, charthick=1.5

n_bh= fload_npart(5)
bhid= fload_blackhole_id(1)

for i=0, n_bh-1 do begin
	;usersym,2.0*cos(findgen(49)/49*2*!pi),2.0*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[1]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[1]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor


; 
; xy
;
;------
!p.position=[x0,y0,x1,y1]

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'


xyouts, x0+0.02, y1-0.04, '!6xz', size=1.2, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[2]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[2]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor





;==================================================================================




; and BH properties
;--------------------

; needs .run bh_multi

x0= 0.65
x1= 0.95
y0= 0.20
y1= 0.80

;yaxistitle="!6Black Hole Mass (!8h!6!E-1!N M!D!9n!6!N)"
yaxistitle="!6Accretion Rate (M!D!9n!6!N Yr!E-1!N)"
;yaxistitle="!6Accretion Rate (Eddington)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

xmin= 0.0
xmax= 3.0

ymax= 5.0e+2
ymin= 1.0e-5


!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, /ylog, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        ytickformat='exp_label', xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase, /normal


;-------------------------------------
;   Get BH Info from fid.bh file
;-------------------------------------
open_blackholes_file, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass


; physical units
bhtime = bhtime / 0.7
bh_mass= 1.0e+10*bh_mass
bh_totalmass= 1.0e+10*bh_totalmass


;if lineary then bh_mass= bh_mass/1.0e7
;oplot, bhtime, bh_mass, psym=-3, color= zcolor, linestyle= lstyle, thick= 4.0

oplot, bhtime, bh_mdot_sunyr, thick=4.0, psym=-3, color= 0, linestyle= 0

idx=where(bhtime ge (fload_time(1)/0.7))
if idx(0) eq -1 then idx[0]= n_elements(sfrtime)-1
oplot, [bhtime[idx[0]]], [bh_mdot_sunyr[idx[0]]], psym= 4, color= 150, thick=5.0, symsize=1.6


;oplot, bhtime, bh_mdot_edd, thick=4.0, psym=-3, color= zcolor, linestyle= lstyle


; print bh_mass
;
idx=where(bhtime ge fload_time(1))
bhmass= alog10(bh_mass(idx(0)))
lblz = strcompress(string(bhmass),/remove_all)
digz= 4
lblz = strmid(lblz,0,digz)        ; 0.xxx
xyouts, x0+0.04, y1-0.12, "Log M!DBH!N= "+lblz+" M!D!9n!6!N", /normal, color= 0, charthick=1.5, size=1.5




;------------------
device, /close

end




;==================================================================================




;==================================================================================
;
;
;   -------------------------------
;   |              |              |
;   |              |              |    ---------------------
;   |              |              |    |                   |
;   |              |              |    |                   |
;   |              |              |    |     SFR           |
;   |              |              |    |                   |
;   -------------------------------    |-------------------|
;   |              |              |    |                   |
;   |              |              |    |     BH Acc        |
;   |              |              |    |                   |
;   |              |              |    |                   |
;   |              |              |    ---------------------
;   |              |              |
;   -------------------------------
; 
;
;
;
;
;

pro img4_both, frun, snapnum, filename=filename
        

if not keyword_set(frun) then begin
   print, "  "
   print, " img4_both, frun, snapnum, filename=filename, "
   print, "  "
   print, "  "
   print, "  "
   return
endif   
        
if not keyword_set(snapnum) then snapnum=12
if not keyword_set(filename) then filename="img_both.eps"




;--------------------------------------------

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=36, newysize=24

;--------------------------------------------



; now images
;-----------------

xlen= 50.0

x0= 0.01
x1= 0.33
x2= 0.65

y0= 0.015
y1= 0.50
y2= 0.985



ok=fload_snapshot_bh(frun,snapnum)

n_bh= fload_npart(5)
bhid= fload_blackhole_id(1)

orig_center= [0,0,0]



;==================================================================================


fload_newcolortable, 1


;
;   Gas
;
; ----------------------------------------------

npart= fload_npart(0)
N= long(npart)

center= orig_center
x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
m= fload_gas_mass(1) / 0.7
hsml= fload_gas_hsml(1) / 0.7



;set_maxden= 10.0  
set_maxden= 1.0    ; 10^4 msolar/pc2
;set_maxden= 0.5e-1

;set_dynrng= 1.0e+3
;set_dynrng= 1.0e+4
set_dynrng= 1.0e+5
;set_dynrng= 1.0e+6
;set_dynrng= 1.0e+6


xz= 0
center= orig_center

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, hsml=hsml, center=center, $
        pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        set_maxden= set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

        XYImage= NxNImage


xz= 1
center= orig_center

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, hsml=hsml, center=center, $
        pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        set_maxden= set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

        XZImage= NxNImage



tv, XYImage, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal
tv, XZImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


xmin= -xlen
xmax=  xlen
ymin= -xlen
ymax=  xlen


; 
; xy
;
;---------------------
fload_newcolortable, 4
!p.position=[x0,y1,x1,y2]
!p.ticklen=0.03

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'

xyouts, x0+0.02, y2-0.04, '!6gas - xy', size=1.5, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	;usersym,2.0*cos(findgen(49)/49*2*!pi),2.0*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[1]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[1]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor


; 
; xz
;
;---------------------
!p.position=[x0,y0,x1,y1]

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'


xyouts, x0+0.02, y1-0.04, '!6gas - xz', size=1.5, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[2]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[2]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor





;==================================================================================


x0= x1 & x1= x2

fload_newcolortable, 4

;
;   Stars
;
; ----------------------------------------------

npart= fload_npart(2)+fload_npart(3)+fload_npart(4)
N= long(npart)

center= orig_center
x= fload_allstars_xyz('x',center=center) / 0.7
y= fload_allstars_xyz('y',center=center) / 0.7
z= fload_allstars_xyz('z',center=center) / 0.7
m= fload_allstars_mass(1) / 0.7
;hsml= fload_allstars_hsml(1) / 0.7
hsml= m*0.0 + 0.2



xz= 0
center= orig_center

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, hsml=hsml, center=center, $
        pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        set_maxden= set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

        XYImage= NxNImage


xz= 1
center= orig_center

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, hsml=hsml, center=center, $
        pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        set_maxden= set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

        XZImage= NxNImage



tv, XYImage, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal
tv, XZImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


xmin= -xlen
xmax=  xlen
ymin= -xlen
ymax=  xlen


; 
; xy
;
;---------------------
!p.position=[x0,y1,x1,y2]
!p.ticklen=0.03

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'

xyouts, x0+0.02, y2-0.04, '!6stars - xy', size=1.5, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	;usersym,2.0*cos(findgen(49)/49*2*!pi),2.0*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[1]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[1]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor


; 
; xz
;
;---------------------
!p.position=[x0,y0,x1,y1]

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'


xyouts, x0+0.02, y1-0.04, '!6stars - xz', size=1.5, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[2]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[2]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor







;==================================================================================






; and SFR
;-----------------

; needs .run sfr_multi

x0= 0.72
x1= 0.98
y0= 0.50
y1= 0.90

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

xmin= 0.0
xmax= 3.0
ymax= 300.0
ymin= 0.05


!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, /ylog, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        ytickformat='exp_label', xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase, /normal
        ;ytickformat='exp_label', xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase, /normal



open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac

; physical units
sfrtime = sfrtime / 0.7

oplot, sfrtime, sfrsfr, psym=-3, color= 0, linestyle= 0, thick= 4.0

idx=where(sfrtime ge (fload_time(1)/0.7))
if idx(0) eq -1 then idx[0]= n_elements(sfrtime)-1
oplot, [sfrtime[idx[0]]], [sfrsfr[idx[0]]], psym= 4, color= 150, thick=5.0, symsize=1.6


sfrtot= total(fload_gas_sfr(1))
lblz = strcompress(string(sfrtot),/remove_all)
digz= 5
if sfrtot lt 100.0 then digz= 4
if sfrtot lt 10.0 then digz= 3
lblz = strmid(lblz,0,digz)        ; 0.xxx
xyouts, x0+0.04, y0+0.04, lblz+" M!D!9n!6!N Yr!E-1!N", /normal, color= 0, charthick=1.5, size=1.5





; and BH properties
;--------------------

; needs .run bh_multi

x0= 0.72
x1= 0.98
y0= 0.10
y1= 0.50

;yaxistitle="!6Black Hole Mass (!8h!6!E-1!N M!D!9n!6!N)"
;yaxistitle="!6Accretion Rate (M!D!9n!6!N Yr!E-1!N)"
;yaxistitle="!6Accretion Rate (Eddington)"
yaxistitle="!6Accretion Rate"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

xmin= 0.0
xmax= 3.0

ymax= 5.0e+0
ymin= 1.0e-5


!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, /ylog, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        ytickformat='exp_label', xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase, /normal


;-------------------------------------
;   Get BH Info from fid.bh file
;-------------------------------------
open_blackholes_file, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass, bh_mdot_edd


; physical units
bhtime = bhtime / 0.7
bh_mass= 1.0e+10*bh_mass
bh_totalmass= 1.0e+10*bh_totalmass

n=n_elements(bhtime)

bhtime= (transpose(bhtime))[2:n-1]
bh_mdot_sunyr= (transpose(bh_mdot_sunyr))[2:n-1]
bh_mdot_edd= (transpose(bh_mdot_edd))[2:n-1]

oplot, bhtime, bh_mdot_sunyr, thick=4.0, psym=-3, color= 100, linestyle= 0
xyouts, x0+0.04, y1-0.22, "M!D!9n!6!N Yr!E-1!N", /normal, color= 0, charthick=1.5, size=1.5

idx=where(bhtime ge (fload_time(1)/0.7))
if idx(0) eq -1 then idx[0]= n_elements(sfrtime)-1
oplot, [bhtime[idx[0]]], [bh_mdot_sunyr[idx[0]]], psym= 4, color= 0, thick=5.0, symsize=1.6



oplot, bhtime, bh_mdot_edd, thick=4.0, psym=-3, color= 200, linestyle= 0
xyouts, x0+0.04, y1-0.06, "Eddington", /normal, color= 0, charthick=1.5, size=1.5

idx=where(bhtime ge (fload_time(1)/0.7))
if idx(0) eq -1 then idx[0]= n_elements(sfrtime)-1
oplot, [bhtime[idx[0]]], [bh_mdot_edd[idx[0]]], psym= 4, color= 0, thick=5.0, symsize=1.6



; print bh_mass
;
idx=where(bhtime ge fload_time(1))
bhmass= alog10(bh_mass(idx(0)))
lblz = strcompress(string(bhmass),/remove_all)
digz= 4
lblz = strmid(lblz,0,digz)        ; 0.xxx
xyouts, x0+0.04, y0+0.04, "Log M!DBH!N= "+lblz+" M!D!9n!6!N", /normal, color= 0, charthick=1.5, size=1.5

;------------------
device, /close

end








;==================================================================================





;==================================================================================
;
;
;   -------------------------------
;   |              |              |    ---------------------
;   |              |              |    |     SFR           |
;   |              |              |    |                   |
;   |              |              |    |                   |
;   |              |              |    |-------------------|
;   |              |              |    |                   |
;   -------------------------------    |    BH Mass        |
;   |              |              |    |                   |
;   |              |              |    |-------------------|
;   |              |              |    |                   |
;   |              |              |    |     BH Acc        |
;   |              |              |    |                   |
;   |              |              |    ---------------------
;   -------------------------------
; 
;
;
;
;
;

pro img4_3, frun, snapnum, filename=filename
        

if not keyword_set(frun) then begin
   print, "  "
   print, " img4_3, frun, snapnum, filename=filename, "
   print, "  "
   print, "  "
   print, "  "
   return
endif   
        
if not keyword_set(snapnum) then snapnum=12
if not keyword_set(filename) then filename="img_both.eps"




;--------------------------------------------

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=36, newysize=24

;--------------------------------------------



; now images
;-----------------

xlen= 70.0

x0= 0.01
x1= 0.33
x2= 0.65

y0= 0.015
y1= 0.50
y2= 0.985



ok=fload_snapshot_bh(frun,snapnum,/nopot_in_snap)

n_bh= fload_npart(5)
bhid= fload_blackhole_id(1)

orig_center= [0,0,0]


fancy_center1= 1
;fancy_center1= 0
if fancy_center1 eq 1 then begin

        center_bh= [0,0,0]
        ;center_bh= fload_center_alreadycomp(1)

        ; two black holes
        if fload_npart(5) gt 1 then begin
                bhid= fload_blackhole_id(1)
                bhid1= bhid[0]
                bhid2= bhid[1]
                center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
                print, "Blackhole : ", bhid1, center1
                center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
                print, "Blackhole : ", bhid2, center2
                idx= where(bhid eq min(bhid))
                if idx(0) eq 0 then center_bh= center1 else center_bh= center2
        endif

        ; one black hole
        if fload_npart(5) eq 1 then begin
                bhid= fload_blackhole_id(1)
                center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
                print, "Blackhole : ", bhid, center_bh
        endif

        orig_center= center_bh
        print, "using this center= ", orig_center

endif





;==================================================================================


fload_newcolortable, 1


;
;   Gas
;
; ----------------------------------------------

npart= fload_npart(0)
N= long(npart)

center= orig_center
x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
;x= fload_gas_xyz('x',center=[0,0,0]) / 0.7
;y= fload_gas_xyz('y',center=[0,0,0]) / 0.7
;z= fload_gas_xyz('z',center=[0,0,0]) / 0.7
m= fload_gas_mass(1) / 0.7
hsml= fload_gas_hsml(1) / 0.7



;set_maxden= 10.0  
set_maxden= 1.0    ; 10^4 msolar/pc2
;set_maxden= 0.5e-1

;set_dynrng= 1.0e+3
;set_dynrng= 1.0e+4
set_dynrng= 1.0e+5
;set_dynrng= 1.0e+6
;set_dynrng= 1.0e+6


xz= 0
center= [0,0,0]     ; if account for center when loading
;center= orig_center
contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, hsml=hsml, center=center, $
        pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        set_maxden= set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

        XYImage= NxNImage


xz= 1
center= [0,0,0]     ; if account for center when loading
;center= orig_center
contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, hsml=hsml, center=center, $
        pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        set_maxden= set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

        XZImage= NxNImage



tv, XYImage, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal
tv, XZImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


xmin= -xlen
xmax=  xlen
ymin= -xlen
ymax=  xlen


; 
; xy
;
;---------------------
fload_newcolortable, 4
!p.position=[x0,y1,x1,y2]
!p.ticklen=0.03

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'

xyouts, x0+0.02, y2-0.04, '!6gas - xy', size=1.5, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	;usersym,2.0*cos(findgen(49)/49*2*!pi),2.0*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[1]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[1]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor

;sidelen= 2.0 * xlen
sidelen= 2.0 * xlen / 0.7
xlenlbl= strcompress(string(sidelen),/remove_all)
if (sidelen) ge 1.0 then digs= 1
if (sidelen) ge 10.0 then digs= 2
if (sidelen) ge 100.0 then digs= 3
;xlenlbl = '!94!6'+strmid(xlenlbl,0,digs)+' kpc/h !96!6'        ; T=0.x (4+digits after decimal)
xlenlbl = '!94!6'+strmid(xlenlbl,0,digs)+' kpc !96!6'        ; T=0.x (4+digits after decimal)
xyouts, x1-0.10, y1+0.03, xlenlbl, /normal, size= 1.4, charthick=3.0, color= 0

oplot_bh_separation, x0+0.02, y1+0.055


; 
; xz
;
;---------------------
!p.position=[x0,y0,x1,y1]

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'


xyouts, x0+0.02, y1-0.04, '!6gas - xz', size=1.5, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[2]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[2]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor

oplot_bh_separation, x0+0.02, y0+0.055, /xz, /projonly




;==================================================================================


x0= x1 & x1= x2

fload_newcolortable, 4

;
;   Stars
;
; ----------------------------------------------

npart= fload_npart(2)+fload_npart(3)+fload_npart(4)
N= long(npart)

center= orig_center
x= fload_allstars_xyz('x',center=center) / 0.7
y= fload_allstars_xyz('y',center=center) / 0.7
z= fload_allstars_xyz('z',center=center) / 0.7
;x= fload_allstars_xyz('x',center=[0,0,0]) / 0.7
;y= fload_allstars_xyz('y',center=[0,0,0]) / 0.7
;z= fload_allstars_xyz('z',center=[0,0,0]) / 0.7
m= fload_allstars_mass(1) / 0.7
;hsml= fload_allstars_hsml(1) / 0.7
hsml= m*0.0 + 0.2



xz= 0
center= [0,0,0]     ; if account for center when loading
;center= orig_center
contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, hsml=hsml, center=center, $
        pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        set_maxden= set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

        XYImage= NxNImage


xz= 1
center= [0,0,0]     ; if account for center when loading
;center= orig_center
contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, hsml=hsml, center=center, $
        pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        set_maxden= set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

        XZImage= NxNImage



tv, XYImage, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal
tv, XZImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


xmin= -xlen
xmax=  xlen
ymin= -xlen
ymax=  xlen


; 
; xy
;
;---------------------
!p.position=[x0,y1,x1,y2]
!p.ticklen=0.03

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'

xyouts, x0+0.02, y2-0.04, '!6stars - xy', size=1.5, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	;usersym,2.0*cos(findgen(49)/49*2*!pi),2.0*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[1]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[1]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor


; 
; xz
;
;---------------------
!p.position=[x0,y0,x1,y1]

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'


xyouts, x0+0.02, y1-0.04, '!6stars - xz', size=1.5, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[2]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[2]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor






;==================================================================================
;
;
;    OK, now we do the data panels on the right
;
;
;
;==================================================================================


snaplbl=strcompress(string(snapnum),/remove_all)
xyouts, 0.72, 0.95, fload_fid(1)+' ('+snaplbl+')', /normal, color= 0, charthick=2.0, size=1.5




; and SFR
;-----------------

; needs .run sfr_multi

x0= 0.72
x1= 0.98
y0= 0.65
y1= 0.93

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

xmin= 0.0
xmax= 3.0
;ymax= 900.0
ymax= 600.0
;ymax= 300.0
ymin= 0.05

open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac

; physical units
sfrtime = sfrtime / 0.7


;xmax= float(long(max(sfrtime)))    ; adjust for simulation
;xmax= 0.5 + float(long(max(sfrtime)))    ; adjust for simulation
xmax= 1.0 + float(long(max(sfrtime)))    ; adjust for simulation

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, /ylog, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        ytickformat='exp_label', xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase, /normal
        ;ytickformat='exp_label', xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase, /normal



oplot, sfrtime, sfrsfr, psym=-3, color= 0, linestyle= 0, thick= 4.0

idx=where(sfrtime ge (fload_time(1)/0.7-1.0e-6))
if idx(0) eq -1 then idx[0]= n_elements(sfrtime)-1
oplot, [sfrtime[idx[0]]], [sfrsfr[idx[0]]], psym= 4, color= 150, thick=10.0, symsize=1.6


sfrtot= total(fload_gas_sfr(1))
lblz = strcompress(string(sfrtot),/remove_all)
digz= 5
if sfrtot lt 100.0 then digz= 4
if sfrtot lt 10.0 then digz= 3
lblz = strmid(lblz,0,digz)        ; 0.xxx
xyouts, x0+0.04, y0+0.04, lblz+" M!D!9n!6!N Yr!E-1!N", /normal, color= 0, charthick=2.0, size=1.5





;--------------------
;  BH Masses
;

; needs .run bh_details

x0= 0.72
x1= 0.98
y0= 0.37
y1= 0.65

;yaxistitle="!6Black Hole Mass (!8h!6!E-1!N M!D!9n!6!N)"
;yaxistitle="!6Black Hole Mass (M!D!9n!6!N)"
yaxistitle="!6Log M!DBH!N (M!D!9n!6!N)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

xmin= 0.0
;xmax= 3.0   ; set above

;ymax= 8.0e+8
;ymin= 6.0e+5
;ymax= 9.7
ymax= 9.2
;ymax= 8.7
ymin= 6.9
;ymin= 6.3


!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $ ;/ylog, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase, /normal
        ;ytickformat='exp_label', xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase, /normal


;-------------------------------------
;   Get BH Info from bh_#.txt files
;-------------------------------------

n_bh= fload_npart(5)
bhid= fload_blackhole_id(1)

;
; bh #1
;
bhidlbl1= strcompress(string(bhid[0]),/remove_all)
read_bh_file, frun, bhid[0], time1, bhm1, pm1, mdot1, mdotedd1, bhfile=frun+"/bh_"+bhidlbl1+".txt", mergertime=mergertime
time1= time1 / 0.7
bhm1= alog10(bhm1 * 1.0e+10 / 0.7)
oplot, time1, bhm1, thick=4.0, psym=-3, color= 0, linestyle= 1
idx=where(time1 ge (fload_time(1)/0.7-0.01))
if idx(0) eq -1 then idx[0]= n_elements(time1)-1
oplot, [time1[idx[0]]], [bhm1[idx[0]]], psym= 4, color= 150, thick=10.0, symsize=1.6

;
; bh #2
;
if n_bh gt 1 then begin
	bhidlbl2= strcompress(string(bhid[1]),/remove_all)
	read_bh_file, frun, bhid[1], time2, bhm2, pm2, mdot2, mdotedd2, bhfile=frun+"/bh_"+bhidlbl2+".txt", mergertime=mergertime
	time2= time2 / 0.7
	bhm2= alog10(bhm2 * 1.0e+10 / 0.7)
	oplot, time2, bhm2, thick=4.0, psym=-3, color= 0, linestyle= 2
	idx=where(time2 ge (fload_time(1)/0.7-0.01))
	if idx(0) eq -1 then idx[0]= n_elements(time2)-1
	oplot, [time2[idx[0]]], [bhm2[idx[0]]], psym= 4, color= 150, thick=10.0, symsize=1.6
endif



;---------------
;  *.bh file
;---------------
open_blackholes_file, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass, bh_mdot_edd
; physical units
bhtime = bhtime / 0.7
bh_mass= alog10(1.0e+10*bh_mass/0.7)
bh_totalmass= 1.0e+10*bh_totalmass

n=n_elements(bhtime)
bhtime= (transpose(bhtime))[2:n-1]
bh_mdot_sunyr= (transpose(bh_mdot_sunyr))[2:n-1]
bh_mdot_edd= (transpose(bh_mdot_edd))[2:n-1]

;oplot, bhtime, bh_mass, thick=4.0, psym=-3, color= 100, linestyle= 0

idx=where(bhtime ge (fload_time(1)/0.7-1.0e-6))
if idx(0) eq -1 then idx[0]= n_elements(sfrtime)-1
;oplot, [bhtime[idx[0]]], [bh_mass[idx[0]]], psym= 4, color= 50, thick=5.0, symsize=1.6






;----------------------
;   BH Accretion Rate
;

; needs .run bh_details

x0= 0.72
x1= 0.98
y0= 0.09
y1= 0.37

;yaxistitle="!6Accretion Rate (M!D!9n!6!N Yr!E-1!N)"
;yaxistitle="!6Accretion Rate (Eddington)"
yaxistitle="!6Accretion Rate"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

xmin= 0.0
; xmax= 3.0  ; set above

ymax= 5.0e+0
ymin= 1.0e-5


!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, /ylog, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        ytickformat='exp_label', xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase, /normal



n_bh= fload_npart(5)
bhid= fload_blackhole_id(1)

;
; bh #1
;
bhidlbl1= strcompress(string(bhid[0]),/remove_all)
read_bh_file, frun, bhid[0], time1, bhm1, pm1, mdot1, mdotedd1, bhfile=frun+"/bh_"+bhidlbl1+".txt", mergertime=mergertime
time1= time1 / 0.7
oplot, time1, mdotedd1, thick=4.0, psym=-3, color= 0, linestyle= 1
idx=where(time1 ge (fload_time(1)/0.7-1.0e-6))
if idx(0) eq -1 then idx[0]= n_elements(time1)-1
oplot, [time1[idx[0]]], [mdotedd1[idx[0]]], psym= 4, color= 150, thick=10.0, symsize=1.6

;
; bh #2
;
if n_bh gt 1 then begin
	bhidlbl2= strcompress(string(bhid[1]),/remove_all)
	read_bh_file, frun, bhid[1], time2, bhm2, pm2, mdot2, mdotedd2, bhfile=frun+"/bh_"+bhidlbl2+".txt", mergertime=mergertime
	time2= time2 / 0.7
	oplot, time2, mdotedd2, thick=4.0, psym=-3, color= 0, linestyle= 2
	idx=where(time2 ge (fload_time(1)/0.7-1.0e-6))
	if idx(0) eq -1 then idx[0]= n_elements(time2)-1
	oplot, [time2[idx[0]]], [mdotedd2[idx[0]]], psym= 4, color= 150, thick=10.0, symsize=1.6
endif





;-------------------------------------
;   Get BH Info from fid.bh file
;-------------------------------------
open_blackholes_file, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass, bh_mdot_edd


; physical units
bhtime = bhtime / 0.7
bh_mass= 1.0e+10*bh_mass
bh_totalmass= 1.0e+10*bh_totalmass

n=n_elements(bhtime)

bhtime= (transpose(bhtime))[2:n-1]
bh_mdot_sunyr= (transpose(bh_mdot_sunyr))[2:n-1]
bh_mdot_edd= (transpose(bh_mdot_edd))[2:n-1]

;oplot, bhtime, bh_mdot_sunyr, thick=4.0, psym=-3, color= 100, linestyle= 0
;xyouts, x0+0.04, y1-0.22, "M!D!9n!6!N Yr!E-1!N", /normal, color= 0, charthick=1.5, size=1.5

idx=where(bhtime ge (fload_time(1)/0.7))
if idx(0) eq -1 then idx[0]= n_elements(sfrtime)-1
;oplot, [bhtime[idx[0]]], [bh_mdot_sunyr[idx[0]]], psym= 4, color= 0, thick=5.0, symsize=1.6



;oplot, bhtime, bh_mdot_edd, thick=4.0, psym=-3, color= 200, linestyle= 0
;xyouts, x0+0.04, y1-0.06, "Eddington", /normal, color= 0, charthick=1.5, size=1.5

idx=where(bhtime ge (fload_time(1)/0.7))
if idx(0) eq -1 then idx[0]= n_elements(sfrtime)-1
;oplot, [bhtime[idx[0]]], [bh_mdot_edd[idx[0]]], psym= 4, color= 0, thick=5.0, symsize=1.6




;------------------
device, /close

end








;==================================================================================
