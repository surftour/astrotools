;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Determine Disk and Bulge Components
;   -----------------------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------







pro bd, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "bd, junk"
   print, "  "
   return
endif


filename='jybd.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; mass ratio
xaxistitle= "!6M!Dpri!N / M!Dsat!N"
xmax = 0.1
xmin = 9.0

; bulge-to-disk
yaxistitle= 'B/D'
ymax = 2.0
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.15
y0= 0.15
x1= 0.98
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
	color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle


; -----------------------------------------------

massr= [95.2/95.2, 95.2/51.0, 95.2/23.25, 95.2/11.9]

fruns= ['Sbfg0.4Sbfg0.4_000', $
	'Sbfg0.4Scfg0.4_000', $
	'Sbfg0.4Sdfg0.4_000', $
	'Sbfg0.4Imfg0.4_000']

bd= grab_bd_ratio(fruns)
oplot, massr, bd, psym=-2, color= 50, thick=3.0


fruns= ['Sbfg0.4Sbfg0.4_030', $
        'Sbfg0.4Scfg0.4_030', $
        'Sbfg0.4Sdfg0.4_030', $
        'Sbfg0.4Imfg0.4_030']

bd= grab_bd_ratio(fruns)
oplot, massr, bd, psym=-7, color= 0, thick=3.0


fruns= ['Sbfg0.4Sbfg0.4_090', $
        'Sbfg0.4Scfg0.4_090', $
        'Sbfg0.4Sdfg0.4_090', $
        'Sbfg0.4Imfg0.4_090']

bd= grab_bd_ratio(fruns)
oplot, massr, bd, psym=-6, color= 200, thick=3.0


fruns= ['Sbfg0.4Sbfg0.4_150', $
        'Sbfg0.4Scfg0.4_150', $
        'Sbfg0.4Sdfg0.4_150', $
        'Sbfg0.4Imfg0.4_150']

bd= grab_bd_ratio(fruns)
oplot, massr, bd, psym=-5, color= 150, thick=3.0


fruns= ['Sbfg0.4Sbfg0.4_180', $
        'Sbfg0.4Scfg0.4_180', $
        'Sbfg0.4Sdfg0.4_180', $
        'Sbfg0.4Imfg0.4_180']

bd= grab_bd_ratio(fruns)
oplot, massr, bd, psym=-1, color= 100, thick=3.0





;------------------------------------

x=[1.0,1.0]
y=[ymin,ymax]
oplot, x, y, psym=-3, linestyle= 1, color= 0




; done
; -----
device, /close


end



;====================================================================================



function grab_bd_ratio, fruns

   numfrun= n_elements(fruns)

   bdr= fltarr(numfrun)

   for i=0, numfrun-1 do begin

	filename='/raid4/tcox/minor/'+fruns[i]+'/bulge.txt'

	read_bulge_file, filename, m, m_thin, m_other, m_nonrot, $
                        g1osm, g1_os_m_thin, g1_os_m_other, g1_os_m_nonrot, $
                        g1nsm, g1_ns_m_thin, g1_ns_m_other, g1_ns_m_nonrot, $
                        g2osm, g2_os_m_thin, g2_os_m_other, g2_os_m_nonrot, $
                        g2nsm, g2_ns_m_thin, g2_ns_m_other, g2_ns_m_nonrot

	bdr[i]= m_nonrot / (m_thin + m_other)
	;bdr[i]= m_nonrot / (m_thin)
	;bdr[i]= m_nonrot / m

   endfor

   return, bdr

end





;====================================================================================



function grab_bt_ratio, fruns

   numfrun= n_elements(fruns)

   btr= fltarr(numfrun)

   for i=0, numfrun-1 do begin

        filename='/raid4/tcox/minor/'+fruns[i]+'/bulge.txt'

        read_bulge_file, filename, m, m_thin, m_other, m_nonrot, $
                        g1osm, g1_os_m_thin, g1_os_m_other, g1_os_m_nonrot, $
                        g1nsm, g1_ns_m_thin, g1_ns_m_other, g1_ns_m_nonrot, $
                        g2osm, g2_os_m_thin, g2_os_m_other, g2_os_m_nonrot, $
                        g2nsm, g2_ns_m_thin, g2_ns_m_other, g2_ns_m_nonrot

        btr[i]= m_nonrot / m

   endfor

   return, btr


end







;====================================================================================
;
;
;   Read the txt file
;
;
;
;
pro read_bulge_file, filename, m, m_thin, m_other, m_nonrot, $
			g1osm, g1_os_m_thin, g1_os_m_other, g1_os_m_nonrot, $
			g1nsm, g1_ns_m_thin, g1_ns_m_other, g1_ns_m_nonrot, $
			g2osm, g2_os_m_thin, g2_os_m_other, g2_os_m_nonrot, $
			g2nsm, g2_ns_m_thin, g2_ns_m_other, g2_ns_m_nonrot


openr, 1, filename
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & m= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & m_thin= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & m_other= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & m_nonrot= float(tempjunk(0))
readf, 1, junk
readf, 1, junk
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1osm= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1_os_m_thin= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1_os_m_other= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1_os_m_nonrot= float(tempjunk(0))
readf, 1, junk
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1nsm= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1_ns_m_thin= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1_ns_m_other= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1_ns_m_nonrot= float(tempjunk(0))
readf, 1, junk
readf, 1, junk
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2osm= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2_os_m_thin= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2_os_m_other= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2_os_m_nonrot= float(tempjunk(0))
readf, 1, junk
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2nsm= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2_ns_m_thin= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2_ns_m_other= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2_ns_m_nonrot= float(tempjunk(0))

close, 1




end










;====================================================================================






;====================================================================================



