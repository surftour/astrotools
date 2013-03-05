pro bh_multi, junk, filename=filename, suppressfid=suppressfid, $
		addthumbnails=addthumbnails, fillinregion=fillinregion, $
		h=h, bhmsg=bhmsg

if not keyword_set(junk) then begin
   print, "  "
   print, "bh, junk, filename=filename, /suppressfid,"
   print, "             /addthumbnails (you must set from within), "
   print, "             /fillinregion, /h, bhmsg=bhmsg"
   print, "  "
   print, "  * automatically goes to bhmass.eps and bhar.eps"
   return
endif



if not keyword_set(bhmsg) then bhmsg= ' '
if not keyword_set(filename) then filename='bhmass.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;-------------------------------------
;  Which runs do we do up?
;-------------------------------------

;fruns= ["Sbc201a_wBH", "Sbc201b_wBH"]
;       lbls= ['q=0.05, w/BH','q=0.25, w/BH']
;       msg= ' '


; -------------------------------------
;  Volker's Runs

;fruns= ["aA1w","bA1w"] & lbls=['coll_1/vc1/bh/','coll_2/vc1/bh/']
;fruns= ["aA2w","bA2w"] & lbls=['coll_1/vc2/bh/','coll_2/vc2/bh/']
;fruns= ["aA3w","bA3w"] & lbls=['coll_1/vc3/bh/','coll_2/vc3/bh/']
;fruns= ["aA4w","bA4w"] & lbls=['coll_1/vc4/bh/','coll_2/vc4/bh/']

;fruns= ["vc3bvc3b_wBH","vc3bvc3_wBH","vc3bvc2_wBH","vc3bvc1_wBH","vc3b_wBH"]
;fruns= ["vc3vc3_wBH","vc3vc2_wBH","vc3vc1_wBH","vc3_wBH"]
;fruns= ["vc2vc2_wBH","vc2vc1_wBH","vc2_wBH"]
;fruns= ["vc1vc1_wBH","vc1_wBH"]
;msg=''

;fruns= ["As/A3", "As/A3_10x"]
;fruns= ["As/A3", "As/A3_1x", "As/A3_10x"]
;fruns= ["As/A3_1x", "As/A3_10x","As/A3_10xa"]
msg=' '



; ---------------------------------------
; bh mass runs

;fruns= ["vc3avc3a_1", "vc3avc3a_2", "vc3avc3a_3"]
;fruns= ["vc3avc3a_3","vc3avc3a_2","vc3avc3a_1","vc3avc3a_4", "vc3avc3a_6"]
;fruns= ["vc3ba_2_1", "vc3ba_2_2", "vc3ba_2_3"]
;fruns= ["vc3ba_2_3","vc3ba_2_2","vc3ba_2_1","vc3ba_2_4","vc3ba_2_5"]
        ;tmerg=[1.1,1.1,1.1,1.1,1.1]
        ;lbls=['1e3','1e4','1e5','1e6','1e7']
        ;lbls=['1e3 (89%)','1e4 (77%)','1e5 (65%)','1e6 (59%)','1e7 (-)']    ; vc3avc3a_* with gas consumption
        ;lbls=['1e3 (94%)','1e4 (93%)','1e5 (78%)','1e6 (69%)','1e7 (42%)']    ; vc3ba_2_* with gas consumption
        ;msg=''


; remnant bh's

;fruns= ["vc3vc3_wBH", "vc3rem_vc3", "vc3rem_vc3rem"]
;fruns= ["vc3vc3", "vc3rem_vc3","vc3rem_vc2"]
;fruns= ["vc3vc3", "vc3rem_vc3","vc3rem_vc2","vc3rem_vc1","vc3rem_vc1a","vc3rem_vc1b","vc3rem_vc1c"]
;        tmerg=[1.1,1.1,1.1,1.1,1.1,1.1,1.1]
;        msg=''

; vc3 - orbits
; --------------
;fruns= ["vc3vc3b", "vc3vc3c","vc3vc3d","vc3vc3e","vc3vc3f","vc3vc3g","vc3vc3h"]
;fruns= ["vc3vc3i", "vc3vc3j","vc3vc3k","vc3vc3l","vc3vc3m","vc3vc3n","vc3vc3o","vc3vc3p"]
;        tmreg=[1.1,1.1,1.1,1.1,1.1,1.1,1.1]
;        msg=''

;fruns= ["vc3vc3e","vc3vc3h"]
;lbls= [" "]
;msg=' '


; larger vcs
;fruns= ["vc4avc4a", "vc5avc5a", "vc6avc6a"]
;        tmerg=[1.1,1.1,1.1]
;        msg=''


; reposition on pot min - does it matter
;fruns= ["vc4vc4a","vc4vc4b"]
;msg= ' '



; with halogas
;fruns= ["v3c","v3h","v3hh"]
;fruns= ["vc3vc3e","vc3HGe"] & lbls= ["with halo gas", "std"]
;fruns=["vc3HGh","vc3vc3h"] & lbls= ["with halo gas", "std"]



; varying gas fractions (with N excursion)
;fruns= ["vc3vc3e","vc3vc3e_gf6","vc3vc3e_gf8","vc3vc3e_gf8a"]


; -------------------------------------------

;fruns= ["vc3vc3e_2", "vc3vc3e_sb8BH", "vc3vc3e_sb10BH", "vc3vc3e_sb13BH"]
;lbls= ["std", "sb8", "sb10", "sb13"]

;fruns= ["vc3vc3e_2", "sbw/sb10BH", "sbw/sb10BHnofb"]
;lbls= ["std", "sb10", "sb10, no fb"]

;fruns= ["vc3vc3e_2", "sbw/sb22BH", "sbw/sb22BHtr1", "sbw/sb22BHtr2","sbw/sb22BHtr3"]
;lbls= ['std', '20.0', '2.0', '0.2', '0.02']

;fruns= ["vc3vc3e_2", "sbw/sb13BH", "sbw/sb13BHtr1", "sbw/sb13BHtr2","sbw/sb13BHtr3"]
;lbls= ['std', '20.0', '2.0', '0.2', '0.02']

;fruns= ["vc3vc3e_2", "sbw/sb8BH", "sbw/sb8BHtr1", "sbw/sb8BHtr2","sbw/sb8BHtr3"]
;lbls= ['std', '20.0', '2.0', '0.2', '0.02']

;fruns= ["vc3vc3e_2", "sbw/sb10BH", "sbw/sb10BHtr1", "sbw/sb10BHtr2","sbw/sb10BHtr3"]
;lbls= ['std', '20.0', '2.0', '0.2', '0.02']

; -------------------------------------------

;fruns= ["bs/b5e", "bs/b5e_igm1"]
;lbls= ['std', 'with IGM']

; -------------------------------------------

;fruns= ["bs/b5e", "sb10_mass/b5e", "sb8_mass/b5e", "sb13_mass/b5e"]
;lbls= ['no wind', 'sb10', 'sb8', 'sb13']

; -------------------------------------------

;fruns= ["bs/b2e_igm3", "minor/min_30"]
;lbls= ['major merger (1:1)', 'minor merger (1:8)']

; -------------------------------------------


; group A
;fruns= ["grpA_1","grpA_2","grpA_3","grpA_4","grpA_5", "grpA_no","grpA_no_1"]
;lbls= ["1e5","1e2","1e1","1e9","1e8", "grpA_no","grpA_no_1"]
;fruns= ["grpA_3","grpA_2","grpA_1","grpA_5","grpA_4"]
;lbls= ['1e1','1e2','1e5','1e8','1e9']
;lbls= [' ',' ',' ',' ',' ']
;msg= ' '

; -------------------------------------


;process_one_sfr, "priya/v3b1v3b1_e", lcolor=0, lthick= 2.0, msg='v3b1v3b1_e', x0= 0.70, y0= 0.94
;process_one_sfr, "priya/v3b1v3b1_e1", lcolor=0, lthick= 2.0, msg='v3b1v3b1_e1', x0= 0.70, y0= 0.91
;process_one_sfr, "tst/tst1", lcolor=100, lthick= 2.0, msg='tst1', x0= 0.70, y0= 0.85
;process_one_sfr, "tst/tst2", lcolor=50,  lthick= 2.0, msg='tst2', x0= 0.70, y0= 0.82
;process_one_sfr, "tst/tst3", lcolor=200, lthick= 2.0, msg='tst3', x0= 0.70, y0= 0.79
;process_one_sfr, "tst/tst3_hdf", lcolor=100, lthick= 2.0, msg='tst3_hdf', x0= 0.70, y0= 0.76
;process_one_sfr, "tst/tst4", lcolor=150, lthick= 2.0, msg='tst4', x0= 0.70, y0= 0.73
;process_one_sfr, "tst/tst5", lcolor=50, lthick= 2.0, msg='tst5', x0= 0.70, y0= 0.70

;process_one_sfr, "priya/v3b1v3b1_e1", lcolor=0, lthick= 2.0, msg='Gadget', x0= 0.30, y0= 0.40
;process_one_sfr, "tst/tst3", lcolor=150, lthick= 2.0, msg='Arepo', x0= 0.30, y0= 0.35

;fruns= ["priya/v3b1v3b1_e1", "tst/tst1", "tst/tst2", "tst/tst3", "tst/tst3_hdf", "tst/tst4", "tst/tst5"]
;lbls= ["v3b1v3b1_e1", "tst1", "tst2", "tst3", "tst3_hdf", "tst4", "tst5"]
;fruns= ["priya/v3b1v3b1_e1", "tst/tst3", "tst/tst3_hdf", "tst/tst4"]
;lbls= ["v3b1v3b1_e1", "tst3", "tst3_hdf", "tst4"]
;fruns= ["priya/v3b1v3b1_e1", "tst/tst3", "tst/tst3_2", "tst/tst3_3"]
;lbls= ["v3b1v3b1_e1", "tst3", "tst3_2", "tst3_3"]

;fruns= ["priya/v3b1v3b1_e1", "tst/tst3"]
;lbls= ["Gadget", "Arepo"]

; -------------------------------------
; std
;fruns= ["priya/v3b1v3b1_e", $
;	"priya/v3b1v3b1_f", $
;	"priya/v3b1v3b1_k"]
;---------------------------
; runs with correct BH
;fruns= ["priya/v3b1v3b1_e", $
;	"priya/v3b1v3b1_f", $
;	"priya/v3b1v3b1_k", $
;	"priya/v3b1v3b4_e", $
;	"priya/v3b1v3b4_f", $
;	"priya/v3b1v3b4_k", $
;	"priya/v3b2v3b1_e", $
;	"priya/v3b2v3b1_f", $
;	"priya/v3b2v3b1_k"]
;---------------------------
; runs with undermassive BH
;fruns= ["priya/v3b2v3b2_e", $
;	"priya/v3b2v3b2_f", $
;	"priya/v3b2v3b2_k", $
;	"priya/v3b2v3b4_e", $
;	"priya/v3b2v3b4_f", $
;	"priya/v3b2v3b4_k"]
;---------------------------
; runs with overmassive BH
;fruns= ["priya/v3b3v3b1_e", $
;	"priya/v3b3v3b1_f", $
;	"priya/v3b3v3b1_k", $
;	"priya/v3b3v3b2_e", $
;	"priya/v3b3v3b2_f", $
;	"priya/v3b3v3b2_k", $
;	"priya/v3b3v3b3_e", $
;	"priya/v3b3v3b3_f", $
;	"priya/v3b3v3b3_k", $
;	"priya/v3b3v3b4_e", $
;	"priya/v3b3v3b4_f", $
;	"priya/v3b3v3b4_k"]


;fruns= ["priya/v3b1v3b1_e", $
;	"priya/v3b2v3b2_e", $
;	"priya/v3b3v3b3_e"]


;fruns= ["priya/iso_v3b1", $
;	"priya/iso_v3b1_2", $
;	"priya/iso_v3b2", $
;	"priya/iso_v3b3", $
;	"priya/iso_v3b4", $
;	"tst/iso_1", $
;	"tst/iso_2"]

fruns= ["priya/iso_v3b1", $
	"priya/iso_v3b1_2", $
	"tst/iso_1", $
	"tst/iso_2"]


; --------------------------------------------

; local group

;fruns= ["localgroup/v3","localgroup/v3_noigm"]
;lbls= ["!6The Local Group", "Isolated MW + M31"]
;msg= ' '




;-------------------------------------
;  Setup various plot info
;-------------------------------------



; physical units
xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
if keyword_set(h) then begin
	xaxistitle = "!6Time (Gyr)"
endif


plottype= 1
;plottype= 2
;plottype= 3

case plottype of
	1: yaxistitle="!6Black Hole Mass (!8h!6!E-1!N M!D!9n!6!N)"
	2: yaxistitle="!6Accretion Rate (M!D!9n!6!N Yr!E-1!N)"
	3: yaxistitle="!6Accretion Rate (Eddington)"
endcase


;xmax= 12.0
;xmax= 10.0
;xmax= 4.0
;xmax= 3.0
xmax= 2.5
;xmax= 2.0
;xmax= 1.5
;xmax= 0.5
xmin= 0.0

; good for mass - exp_label OK for these
if plottype eq 1 then begin
	;ymax= 7.0e+10
	;ymax= 7.0e+9
	;ymax= 3.0e+9
	;ymax= 3.0e+8
	ymax= 1.2e+8

	;ymin= 2.0e+7
	ymin= 7.0e+6
	;ymin= 6.0e+5
	;ymin= 2.0e+5
	;ymin= 7.0e+4
	;ymin= 1.0e+3
	;ymin= 1.0e+0

	lineary= 0     ; need to flip off ylog, exp_label by hand
	;lineary= 1
	if lineary then begin
	   ymax= 1.1
	   ymin= 0.0
	   yaxistitle="!6Black Hole Mass (10!E7!N h!E-1!N M!D!9n!6!N)"
	endif
endif


; acc. rate (msolar/yr)
if plottype eq 2 then begin
	;ymax= 5.0e+2
	;ymax= 1.0e+1
	ymax= 2.0e+0
	;ymax= 1.0e+0
	ymin= 1.0e-4
	;ymin= 1.0e-5
	;ymin= 1.0e-8
	;ymin= 1.0e-9
endif


; acc. rate (edd)
if plottype eq 3 then begin
	;ymax= 10
	;ymax= 7
	ymax= 4
	;ymin= 1.0e-3
	ymin= 1.0e-4
	;ymin= 1.0e-5
endif



zcolor= 50
deltacolor= (250-zcolor)/n_elements(fruns)
;deltacolor= 50
lstyle= 0

;
;zcolor= 10
;deltacolor= 5
;zcolor= 80
;deltacolor= 10
;zcolor= 150
;deltacolor= 10

zcolor_orig= zcolor

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $ 
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata         



; loop through the files
; -------------------------

ifinal= n_elements(fruns)-1
for i=0,ifinal do begin

	;-------------------------------------
	;   Get BH Info from fid.bh file
	;-------------------------------------
	open_blackholes_file, fruns[i], bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass, bh_mdot_edd

;stop


	; physical units
	if keyword_set(h) then begin
		xaxistitle = "Time (Gyr)"
		h = fload_cosmology('h')
		bhtime = bhtime / h
		bh_mass = bh_mass / h
		bh_totalmass = bh_totalmass / h
	endif

	bh_mass= 1.0e+10*bh_mass
	bh_totalmass= 1.0e+10*bh_totalmass

	;---------------------------
	; determine time of black hole merger
	idx= where(bh_num eq 1)
	if idx(0) ne -1 then timemerge= bhtime(idx(0))

	if plottype eq 1 then begin
		;ytoparrow= 3.0e+8
		;ybottomarrow= 8.0e+8
		ytoparrow= 8.0e+5
		ybottomarrow= 4.0e+5
	endif

	if plottype eq 2 then begin
		ytoparrow= 1.0
		ybottomarrow= 3.0
	endif

	if plottype eq 3 then begin
		ytoparrow= 2.0
		ybottomarrow= 6.0
	endif

	if timemerge gt 0 then begin
            arrow, timemerge, ybottomarrow, timemerge, ytoparrow, $
			COLOR=zcolor, THICK=3.0, hthick=3.0, /data
	endif
	timemerge= 0

	;---------------------------
	; plot the various things
	if plottype eq 1 then begin
	   if lineary then bh_mass= bh_mass/1.0e7
	   oplot, bhtime, bh_mass, psym=-3, color= zcolor, linestyle= lstyle, thick= 4.0
	endif

	if plottype eq 2 then begin
	   oplot, bhtime, bh_mdot_sunyr, thick=4.0, psym=-3, color= zcolor, linestyle= lstyle
	endif 

	if plottype eq 3 then begin
	   oplot, bhtime, bh_mdot_edd, thick=4.0, psym=-3, color= zcolor, linestyle= lstyle
	endif

	lstyle= lstyle+1
	zcolor= zcolor+deltacolor
endfor

; done with loop
; ---------------

xyouts, 0.75, 0.90, bhmsg, /normal, charthick=3.0, size=1.5, color= 0


; std
;std= 0
std= 1

if plottype eq 1 then begin
    x0= 0.23
    y0= 0.94

    if lineary eq 1 then begin
	;x0= 0.65
	;y0= 0.52
	x0= 0.23
	y0= 0.94
    endif
endif

if plottype eq 2 then begin
    ;x0= 0.70
    x0= 0.25
    y0= 0.92
endif

if plottype eq 3 then begin
    x0= 0.66
    y0= 0.92
endif

if std eq 1 then begin
    zcolor= zcolor_orig
    for i=0, n_elements(fruns)-1 do begin
        ;y0= y0-0.04
        y0= y0-0.03
        if n_elements(lbls) gt 0 then begin
                xyouts, x0, y0, lbls[i], /normal, charthick=3, size=1.33, color=zcolor
                    zcolor= zcolor+deltacolor
        endif else begin
                ;xyouts, x0, y0, fruns[i], /normal, charthick=3, size=1.33, color=zcolor
                ;    zcolor= zcolor+deltacolor
                xyouts, x0, y0, fruns[i], /normal, charthick=3, size=1.0, color=zcolor
                    zcolor= zcolor+deltacolor
        endelse
    endfor
endif



;--------------------------------------
;--------------------------------------

device, /close




end










; ================================================================================



pro open_blackholes_file, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass, bh_mdot_edd


	spawn, 'echo $HOME', result
	homedir= strcompress(result,/remove_all)

	;-------------------------------------
	;   Get BH Info from fid.bh file
	;-------------------------------------
	; get bh data
	    loopnum= 0
            repeat begin

		if strmid(frun,1,4) eq 'raid' then begin
			spawn, '/bin/ls '+frun+'/*.bh ',result
		endif else begin
                	case loopnum of
                	  0: datadir='/n/circelfs/hernquist_lab/tcox/'
                	  1: datadir=homedir+'/'
                	  2: datadir='/n/scratch/hernquist_lab/tcox/'
                	  3: datadir='data/'
                	  4: datadir='data1/'
                	  else: break
                	endcase

			nfrun= fload_getid(frun)
			;spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/*.bh ',result
			spawn, '/bin/ls '+datadir+nfrun+'/*.bh ',result
		endelse
    
                bhfile=strcompress(result[0],/remove_all)
    

                get_lun, unit
                openr, unit, bhfile, ERROR=err
                close, unit
                free_lun, unit

                if (err NE 0) then begin
                   ERR= 0
                   foundit= 0
                endif else begin
                   foundit= 1
                endelse

                loopnum= loopnum+1

            endrep until (foundit eq 1)


            if foundit eq 0 then begin
                msg= 'Problem: file not found  file='+bhfile
                print, "  "
                print, msg
                print, "  "
                return
            endif

		; open it!
		print, "opening: ",bhfile
		bhdata= read_ascii(bhfile)
		bhtime= bhdata.field1[0,*]
		bh_num= bhdata.field1[1,*]
		bh_mass= bhdata.field1[2,*]
		bh_mdot_gu= bhdata.field1[3,*]
		bh_mdot_sunyr= bhdata.field1[4,*]
		bh_totalmass= bhdata.field1[5,*]
		bh_mdot_edd= bhdata.field1[6,*]

end




; ================================================================================








pro bh_at_snaptimes, frun, h=h


if not keyword_set(frun) then begin
   print, "  "
   print, "bh_at_snaptimes, frun, /h"
   print, "  "
   print, "  "
   return
endif



;-------------------------------------------
;-------------------------------------------


; determine the number of
; snapshots in frun directory
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snapshot* | wc ",result
nsnaps=long(result[0])




; ---------------------------------------------------


time= fltarr(nsnaps)

bhmass_inst= fltarr(nsnaps)
bhmass_avg= fltarr(nsnaps)
bhar_inst= fltarr(nsnaps)
bhar_avg= fltarr(nsnaps)
bhared_inst= fltarr(nsnaps)
bhared_avg= fltarr(nsnaps)


ymin = 0 


; physical units
if keyword_set(h) then begin
	h = fload_cosmology('h')
endif


; ---------------------------------------------------

	    spawn, 'echo $HOME', result
	    homedir= strcompress(result,/remove_all)


	;-------------------------------------------
	;   Get BH information from txt
	;-------------------------------------------
	    ; get bh data (from file)
	    ; -------------------------
	    loopnum= 0
	    repeat begin
	    ;bhfile= '/data/tcox/bh/'+fload_getid(fruns[i])+'.bh'        ; twopiter
	    ;bhfile= '/home/tcox/data/bh/'+fload_getid(fruns[i])+'.bh'    ; harvard
	    ;bhfile= '/raid4/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.bh'    ; harvard
	    ;bhfile= '/raid2/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.bh'    ; harvard
	    ;bhfile= '/data6/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.bh'    ; harvard

                ;if strmid(frun,1,4) eq 'raid' then begin
                ;print, strmid(frun,0,7)
                if strmid(frun,0,7) eq '/n/data' then begin
                        spawn, '/bin/ls '+frun+'/*.bh ',result
                endif else begin
                        case loopnum of
                          0: datadir=homedir+'/data'
                          1: datadir='/n/circelfs/hernquist_lab'
                          2: datadir=homedir
                          3: datadir='/n/home'
                          4: datadir='/data6'
                          5: datadir='/data7'
                          6: datadir=''
                          else: break
                        endcase

                        nfrun= fload_getid(frun)
                        spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/*.bh ',result
                endelse


		;bhfile= datadir+'/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.bh'

                bhfile=strcompress(result[0],/remove_all)
    

		openr, 1, bhfile, ERROR=err
		close, 1

		if (err NE 0) then begin
		   ERR= 0
		   foundit= 0
		endif else begin
		   foundit= 1
		endelse

		loopnum= loopnum+1

	    endrep until (foundit eq 1)

	    if foundit eq 0 then begin
		msg= 'Problem: bh file no found  file=[data]/tcox/'+fload_getid(frun)+'/'+fload_getid(frun)+'.bh'
	        print, "  "
	        print, msg
	        print, "  "
		return
	    endif

		; open it!
		print, "opening: ",bhfile
		bhdata= read_ascii(bhfile)
		bhtime= bhdata.field1[0,*]
		bh_num= bhdata.field1[1,*]
		bh_mass= bhdata.field1[2,*]
		bh_mdot_gu= bhdata.field1[3,*]
		bh_mdot_sunyr= bhdata.field1[4,*]
		bh_totalmass= bhdata.field1[5,*]
		bh_mdot_edd= bhdata.field1[6,*]



	bh_mass= 1.0e+10*bh_mass
	bh_totalmass= 1.0e+10*bh_totalmass

	;---------------------------


	    ; physical units
	    if keyword_set(h) then bhtime = bhtime / h


; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"


        ; open snapshot
        ;ok=fload_snapshot(frun,snapnum)
        ;ok=fload_snapshot_bh(frun,i,/nopot_in_snap)
        ok=fload_snapshot_bh(frun,i)


        ; what time is it?
        time[i]= fload_time(1)

	idx_gtr_snaptime= where(bhtime ge time[i])
	idx_curr_snaptime= idx_gtr_snaptime[0]


	; instantaneous bh properties
	bhmass_inst[i]= bh_mass(idx_curr_snaptime)
	bhar_inst[i]= bh_mdot_sunyr(idx_curr_snaptime)
	bhared_inst[i]= bh_mdot_edd(idx_curr_snaptime)


	; set time window
	dt= 0.01
	idx_window= where((bhtime ge (time[i]-dt)) and (bhtime le (time[i]+dt)))
	bh_dt= bh_mass(idx_window)
	bhmass_avg[i]= mean(bh_dt)
	bh_dt= bh_mdot_sunyr(idx_window)
	bhar_avg[i]= mean(bh_dt)
	bh_dt= bh_mdot_edd(idx_window)
	bhared_avg[i]= mean(bh_dt)

	print, "T= ", time[i], bhmass_inst[i], bhmass_avg[i]
	print, "          ", bhar_inst[i], bhar_avg[i]
	print, "          ", bhared_inst[i], bhared_avg[i]

endfor


bhmass_avg= alog10(bhmass_avg)
bhar_avg= alog10(bhar_avg)
bhared_avg= alog10(bhared_avg)



;--------------------------------------
;--------------------------------------


; SFR(snaptimes) FILE
; --------------------
openw, 1, frun+'/bh_snaptimes.txt', ERROR=err

printf, 1, "#   bh_snaptimes.txt"
printf, 1, "#   frun: "+frun
printf, 1, "#           inst.    accretion   accretion "
printf, 1, "#   time       mass       rate        rate "
printf, 1, "# (Gyr/h)  Log(M_sun)  Log(Mo /Yr)   Log(Edd) "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F8.3,"     ",3(F8.4,"    "))', $
                time[i], bhmass_avg[i], bhar_avg[i], bhared_avg[i]
endfor
close, 1




; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            "
print, "---------------------------------"
print, "---------------------------------"





end






;=================================================================







pro eachbh_at_snaptimes, frun, h=h


if not keyword_set(frun) then begin
   print, "  "
   print, "eachbh_at_snaptimes, frun, /h"
   print, "  "
   print, "  "
   return
endif



;-------------------------------------------
;-------------------------------------------


; determine the number of
; snapshots in frun directory
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snapshot* | wc ",result
nsnaps=long(result[0])




; ---------------------------------------------------


time= fltarr(nsnaps)

bhmass_inst= fltarr(nsnaps)
bhmass_avg= fltarr(nsnaps)
bhar_inst= fltarr(nsnaps)
bhar_avg= fltarr(nsnaps)
bhared_inst= fltarr(nsnaps)
bhared_avg= fltarr(nsnaps)


ymin = 0 


; physical units
if keyword_set(h) then begin
	h = fload_cosmology('h')
endif


; ---------------------------------------------------


	;-------------------------------------------
	;   Get BH information from txt
	;-------------------------------------------
	    ; get bh data (from file)
	    ; -------------------------
	    loopnum= 0
	    repeat begin
	    ;bhfile= '/data/tcox/bh/'+fload_getid(fruns[i])+'.bh'        ; twopiter
	    ;bhfile= '/home/tcox/data/bh/'+fload_getid(fruns[i])+'.bh'    ; harvard
	    ;bhfile= '/raid4/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.bh'    ; harvard
	    ;bhfile= '/raid2/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.bh'    ; harvard
	    ;bhfile= '/data6/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.bh'    ; harvard

                ;if strmid(frun,1,4) eq 'raid' then begin
                ;print, strmid(frun,0,7)
                if strmid(frun,0,7) eq '/n/data' then begin
                        spawn, '/bin/ls '+frun+'/*.bh ',result
                endif else begin
                        case loopnum of
                          0: datadir='/n/home/tcox/data'
                          1: datadir='/n/circelfs/hernquist_lab'
                          2: datadir='/n/home/tcox'
                          3: datadir='/n/home'
                          4: datadir='/data6'
                          5: datadir='/data7'
                          6: datadir=''
                          else: break
                        endcase

                        nfrun= fload_getid(frun)
                        spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/*.bh ',result
                endelse


		;bhfile= datadir+'/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.bh'

                bhfile=strcompress(result[0],/remove_all)
    

		openr, 1, bhfile, ERROR=err
		close, 1

		if (err NE 0) then begin
		   ERR= 0
		   foundit= 0
		endif else begin
		   foundit= 1
		endelse

		loopnum= loopnum+1

	    endrep until (foundit eq 1)

	    if foundit eq 0 then begin
		msg= 'Problem: bh file no found  file=[data]/tcox/'+fload_getid(frun)+'/'+fload_getid(frun)+'.bh'
	        print, "  "
	        print, msg
	        print, "  "
		return
	    endif

		; open it!
		print, "opening: ",bhfile
		bhdata= read_ascii(bhfile)
		bhtime= bhdata.field1[0,*]
		bh_num= bhdata.field1[1,*]
		bh_mass= bhdata.field1[2,*]
		bh_mdot_gu= bhdata.field1[3,*]
		bh_mdot_sunyr= bhdata.field1[4,*]
		bh_totalmass= bhdata.field1[5,*]
		bh_mdot_edd= bhdata.field1[6,*]



	bh_mass= 1.0e+10*bh_mass
	bh_totalmass= 1.0e+10*bh_totalmass

	;---------------------------


	    ; physical units
	    if keyword_set(h) then bhtime = bhtime / h


; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"


        ; open snapshot
        ;ok=fload_snapshot(frun,snapnum)
        ;ok=fload_snapshot_bh(frun,i,/nopot_in_snap)
        ok=fload_snapshot_bh(frun,i)


        ; what time is it?
        time[i]= fload_time(1)

	idx_gtr_snaptime= where(bhtime ge time[i])
	idx_curr_snaptime= idx_gtr_snaptime[0]


	; instantaneous bh properties
	bhmass_inst[i]= bh_mass(idx_curr_snaptime)
	bhar_inst[i]= bh_mdot_sunyr(idx_curr_snaptime)
	bhared_inst[i]= bh_mdot_edd(idx_curr_snaptime)


	; set time window
	dt= 0.01
	idx_window= where((bhtime ge (time[i]-dt)) and (bhtime le (time[i]+dt)))
	bh_dt= bh_mass(idx_window)
	bhmass_avg[i]= mean(bh_dt)
	bh_dt= bh_mdot_sunyr(idx_window)
	bhar_avg[i]= mean(bh_dt)
	bh_dt= bh_mdot_edd(idx_window)
	bhared_avg[i]= mean(bh_dt)

	print, "T= ", time[i], bhmass_inst[i], bhmass_avg[i]
	print, "          ", bhar_inst[i], bhar_avg[i]
	print, "          ", bhared_inst[i], bhared_avg[i]

endfor


bhmass_avg= alog10(bhmass_avg)
bhar_avg= alog10(bhar_avg)
bhared_avg= alog10(bhared_avg)



;--------------------------------------
;--------------------------------------


; SFR(snaptimes) FILE
; --------------------
openw, 1, frun+'/bh_snaptimes.txt', ERROR=err

printf, 1, "#   bh_snaptimes.txt"
printf, 1, "#   frun: "+frun
printf, 1, "#           inst.    accretion   accretion "
printf, 1, "#   time       mass       rate        rate "
printf, 1, "# (Gyr/h)  Log(M_sun)  Log(Mo /Yr)   Log(Edd) "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F8.3,"     ",3(F8.4,"    "))', $
                time[i], bhmass_avg[i], bhar_avg[i], bhared_avg[i]
endfor
close, 1




; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            "
print, "---------------------------------"
print, "---------------------------------"





end








;======================================================================




pro open_blackhole_txt, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, $
				bh_totalmass, bh_mdot_edd


	    spawn, 'echo $HOME', result
	    homedir= strcompress(result,/remove_all)

	    ; get sfr data (from file)
	    ; -------------------------
	    loopnum= 0
	    repeat begin
	    ;bhfile= '/data/tcox/bh/'+fload_getid(fruns[i])+'.bh'        ; twopiter
	    ;bhfile= '/home/tcox/data/bh/'+fload_getid(fruns[i])+'.bh'    ; harvard
	    ;bhfile= '/raid4/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.bh'    ; harvard
	    ;bhfile= '/raid2/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.bh'    ; harvard
	    ;bhfile= '/data6/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.bh'    ; harvard
  
                if strmid(frun,1,4) eq 'raid' then begin
                        spawn, '/bin/ls '+frun+'/*.bh ',result
                endif else begin
                        case loopnum of
                          0: datadir='/raid4'
                          1: datadir='/home'
                          2: datadir='/raid2'
                          3: datadir='/data'
                          4: datadir='/data6'
                          5: datadir='/data7'
                          else: break
                        endcase

                        nfrun= fload_getid(frun)
                        spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/*.bh ',result
                endelse

		bhfile= strcompress(result[0],/remove_all)

		get_lun, unit
		openr, unit, bhfile, ERROR=err
		close, unit
		free_lun, unit

		if (err NE 0) then begin
		   ERR= 0
		   foundit= 0
		endif else begin
		   foundit= 1
		endelse

		loopnum= loopnum+1

	    endrep until (foundit eq 1)

	    if foundit eq 0 then begin
		msg= 'Problem: bh file no found  file=[data]/tcox/'+fload_getid(frun)+'/'+fload_getid(frun)+'.bh'
	        print, "  "
	        print, msg
	        print, "  "
		return
	    endif

		; open it!
		print, "opening: ",bhfile
		bhdata= read_ascii(bhfile)
		bhtime= bhdata.field1[0,*]
		bh_num= bhdata.field1[1,*]
		bh_mass= bhdata.field1[2,*]
		bh_mdot_gu= bhdata.field1[3,*]
		bh_mdot_sunyr= bhdata.field1[4,*]
		bh_totalmass= bhdata.field1[5,*]
		bh_mdot_edd= bhdata.field1[6,*]



	bh_mass= 1.0e+10*bh_mass
	bh_totalmass= 1.0e+10*bh_totalmass

	;---------------------------


end











