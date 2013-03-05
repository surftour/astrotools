
;===========================================================================
;
; TJ adapted from code written by Patrik, Chris, and Desika
;
;==========================================================================




mcrxfile= 'mcrx_032.fits'



!P.FONT = -1
set_plot,'ps'
device,/encaps,/color,bits_per_pixel=8, $
	preview=0,pre_depth=8,pre_x=3.5,pre_y=3,/inches
device, SET_CHARACTER_SIZE=[300,350],xsize=14,ysize=12
!X.THICK=3
!Y.THICK=3
!P.THICK=2
!X.MARGIN=[8,3]

!P.MULTI=0
colr=[0,255,255,0  ,0  ,255,215,  0,155,  0,155,235]
colg=[0,255,0  ,0  ,235,0  ,215,255,155,130,  0, 95]
colb=[0,255,0  ,255,0  ,255,  0,255,155,186,215,  0]
tvlct,colr,colg,colb


fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7




;; figure out which HDU's

;debug
print,'mcrxfile='+mcrxfile

scatlam_HDU = find_HDU (mcrxfile, "scattering_lambdas")
lambda_HDU = find_HDU (mcrxfile, "lambda")
iq_HDU = find_HDU (mcrxfile, "integrated_quantities")
mcrx_HDU = find_HDU (mcrxfile, "MCRX")

;; get number of cameras
h = HEADFITS (mcrxfile, exten= mcrx_HDU)
n_camera = FXPAR (h, 'n_camera')

print,'number of cameras = '+strcompress(n_camera)

;; get wavelengths
ftab_ext, mcrxfile, ext = iq_HDU, "lambda", lambda
;; the nonscatter SED's have lambda like the images, hokey
ftab_ext, mcrxfile, ext = lambda_HDU, "lambda", lambda_ns

;ftab_ext, mcrxfile, ext = scatlam_HDU, "entry", lambda_entries
;lambda_entries = lambda_entries (1:*)
junk=where(lambda ge 5e-4)
last_lambda=junk[0]-1


; number of camers will set colors, legends, etc.
cols=indgen(n_camera)


colors=indgen(n_elements(cols))+3
colors[0:3]=[c1,c2,c3,c4]



legtxt='w/o dust'
for cc=0,n_elements(cols) - 1 do legtxt = [legtxt, 'camera ' + strcompress (string (cols[cc] ),/remove_all)]


lines=fltarr(n_camera)



;-----------------------------------------------------

    h = HEADFITS (mcrxfile, exten= iq_HDU)
    print, ' Loading ', mcrxfile

    ; load INTEGRATED_QUANTITIES
    d=mrdfits(mcrxfile,iq_HDU)
    title_base= 'SED: '+mcrxfile
                                ;for c = 0,n_camera - 1 do begin
    ; loop through cameras
    for cc = 0,n_elements(cols)-1 do begin
        c=cols[cc]
        ;; get data from fits file
        colstring = strcompress (string (c, format = '(i 2.1)' ),/remove_all)
        colname_ns ="L_LAMBDA_NONSCATTER"+colstring
        colname="L_LAMBDA_SCATTER"+colstring
	colname_ir="L_LAMBDA_IR"+colstring
                                ;ftab_ext, mcrxfile, ext = iq_HDU, colname_ns, L_lambda_ns
                                ;ftab_ext, mcrxfile, ext = iq_HDU, colname, L_lambda
        llnse=where(tag_names(d) eq colname_ns)
        lle=where(tag_names(d) eq colname)

        L_lambda_ns=d.(llnse)
        L_lambda=d.(lle)

	; IR 
	llire=where(tag_names(d) eq colname_ir)
	L_lambda_ir=d.(llire)
	L_lambda_tot=L_lambda_ir+L_lambda
	; set IR to zero
	; L_lambda_ir=0*L_lambda_ns
	; L_lambda_tot=L_lambda


        ;; plot SED
	; if not plotting all cameras on one plot or current is camera 0
        if (c eq 0) then begin
	    plot_name= mcrxfile+"SED-"+ string (c, format = '(i 2.2)' )
            device,filename=plot_name+'.eps'
            yrange = [min(L_lambda* lambda), max([L_lambda_ns*lambda,L_lambda*lambda,L_lambda_ir*lambda])]
            title = title_base +" camera "+ string (c, format = '(i 2.1)' )

	    ; plot the axes etc.
            plot,[1],[1], xrange = [1e-8,1e-3], yrange=yrange,$
              ystyle = 16 ,xtitle=textoidl('\lambda/m'), $
              ytitle=textoidl("\lambdaL_\lambda/W"), /xlog,/ylog,$
              title=title,$
              background = bgc,color=fgc,/nodata

            oplot,lambda_ns, lambda_ns*L_lambda_ns,color=c1

            legend,legtxt, colors=[c1,colors], textcolor=fgc,linestyle=0,/right,box=0
        endif


        oplot,lambda, lambda*L_lambda, color = colors[cc],lines=lines[cc]
;        oplot,lambda, lambda*L_lambda, color = colors[1],lines=lines[1]
	oplot,lambda, lambda*L_lambda_ir, color = c3, lines=lines[2]
	oplot,lambda, lambda*L_lambda_tot, color = c4, lines=lines[3]


    endfor


    print, "Closing PostScript file"
    device,/close



;======================
; Done
;======================

end

