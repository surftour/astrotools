;-------------------------------------
;-------------------------------------

pro setup_plot_stuff, sendto, filename=filename, colortable=colortable, $
				newxsize=newxsize, newysize=newysize


COMMON PlotInfo


if not keyword_set(colortable) then colortable= 0

if (sendto EQ 'ps') then begin
        if not keyword_set(filename) then begin
                filename= 'fp.eps'
                ans= ''
                read, ans, PROMPT='eps filename ['+filename+']:'
                if strlen(ans) GT 0 then filename= ans
        endif

	if keyword_set(newxsize) then nxsize=newxsize else nxsize=14
	if keyword_set(newysize) then nysize=newysize else nysize=14

        set_plot, 'ps'
        device, filename= filename, /encapsulated,/color,bits_per_pixel=8
        device, SET_CHARACTER_SIZE=[200,300], xsize=nxsize, ysize=nysize
        ;!p.ticklen=0.03
        ;device,/times,/italic,font_index=20

	; this is the whole crux of doing this, only have to
	;  load the color table once and then everything uses this.
        if ct_is_loaded eq 0 then begin
		;loadct, 0
		case colortable of
		   1: loadct, 1               ; blue/white
		   2: loadct, 2               ; grn-red-blu-wht
		   3: loadct, 3               ; red temperature
		   4: loadct, 4               ; blue/green/red/yellow (std)
		   else: loadct, 0            ; b-w linear
		endcase
                tvlct,r,g,b,/get
                v1=[0,255]
                v2=[0,255]
                v3=[0,255]
                tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white
                ct_is_loaded= 1
		ct_current= colortable
        endif

endif else begin
        ;set_plot, 'x'
        ;thetitle = "Fix this in idlstuff/plot/setup_plot_stuff.pro"
        ;window,13,xsize=600,ysize=600, title=thetitle
        ;!p.background= getcolor('white')
	print, " "
	print, " sorry, it's only possible to make an eps file."
	print, " please go surfing now."
	print, " "
	print, " "
endelse



end

