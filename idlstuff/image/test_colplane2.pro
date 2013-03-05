

pro gen_datfile, junk


        ;base_filename='test_colplane_flip.jpg'    ; change 0, 1, 2 ->  1, 0, 2
        base_filename='colplane_flip2'    ; change 0, 1, 2 ->  1, 2, 0



	; this is volker's original color plane
	;        (Black,Pink/White/Blue)
        cols=256
        ColPlane=fltarr(cols,cols,3)
        openr,1,"/home2/tcox/Tools/idlstuff/image/colplane2.dat"
        readu,1,ColPlane
        close,1

	; flip
	New_ColPlane= 0.0 * ColPlane

	New_ColPlane(*,*,0)= ColPlane(*,*,1)
	New_ColPlane(*,*,1)= ColPlane(*,*,2)
	New_ColPlane(*,*,2)= ColPlane(*,*,0)


        openw,1,'/home2/tcox/Tools/idlstuff/image/'+base_filename+'.dat'
        writeu,1,New_ColPlane
        close,1


	write_jpg_fromdat, '/home2/tcox/Tools/idlstuff/image/'+base_filename+'.jpg', base_filename+'.dat'


end


   
;======================================================================



pro test_v, junk

   ; Test
   ; Volker's Pink/White/Blue Colors
   ; ----------------------------------
	filename='test_colplane2.jpg'
        cols=256
        ColPlane=fltarr(cols,cols,3)

        openr,1,"/home2/tcox/Tools/idlstuff/image/colplane2.dat"
        readu,1,ColPlane
        close,1


        pixels= 3*cols
        JpegPic=bytarr(pixels,pixels,3)

        for i=0, pixels-1 do begin 
              for j=0, pixels-1 do begin
                 JPegPic(i,j,0)= ColPlane(i/3,j/3,0)
                 JPegPic(i,j,1)= ColPlane(i/3,j/3,1)
                 JPegPic(i,j,2)= ColPlane(i/3,j/3,2)
              endfor
        endfor




   write_jpeg, filename , Jpegpic, true=3, quality=90


   cmd = 'mogrify -font arial -antialias -pointsize 58 -fill white -draw '
   cmd = cmd + "'"
   cmd = cmd + 'text 40,72 '
   cmd = cmd + """
   cmd = cmd + 'colplane2.dat'
   cmd = cmd + """
   cmd = cmd + "'"
   cmd = cmd + ' -quality 95 '
   cmd = cmd + filename
   print, cmd
   spawn, cmd


end




;======================================================================





pro test_vv, datfile

   ; Test
   ; Volker's Pink/White/Blue Colors
   ; ----------------------------------
	;filename='test_colplane_flip.jpg'    ; change 0, 1, 2 ->  1, 0, 2
	filename='test_colplane_flip2.jpg'    ; change 0, 1, 2 ->  1, 2, 0
        cols=256
        ColPlane=fltarr(cols,cols,3)

        openr,1,"/home2/tcox/Tools/idlstuff/image/colplane2.dat"
        readu,1,ColPlane
        close,1


        pixels= 3*cols
        JpegPic=bytarr(pixels,pixels,3)

        for i=0, pixels-1 do begin 
              for j=0, pixels-1 do begin
                 ;JPegPic(i,j,0)= ColPlane(i/3,j/3,0)
                 ;JPegPic(i,j,1)= ColPlane(i/3,j/3,1)
                 ;JPegPic(i,j,2)= ColPlane(i/3,j/3,2)
                 JPegPic(i,j,0)= ColPlane(i/3,j/3,1)
                 JPegPic(i,j,1)= ColPlane(i/3,j/3,2)
                 JPegPic(i,j,2)= ColPlane(i/3,j/3,0)
              endfor
        endfor




   write_jpeg, filename , Jpegpic, true=3, quality=90


   cmd = 'mogrify -font arial -antialias -pointsize 58 -fill white -draw '
   cmd = cmd + "'"
   cmd = cmd + 'text 40,72 '
   cmd = cmd + """
   cmd = cmd + 'colplane_flip2'
   cmd = cmd + """
   cmd = cmd + "'"
   cmd = cmd + ' -quality 95 '
   cmd = cmd + filename
   print, cmd
   spawn, cmd


end







;======================================================================





pro write_jpg_fromdat, jpgfile, datfile

        cols=256
        ColPlane=fltarr(cols,cols,3)
        openr,1,'/home2/tcox/Tools/idlstuff/image/'+datfile
        readu,1,ColPlane
        close,1


	;multfac= 1
	;multfac= 2
	multfac= 3
        pixels= multfac*cols
        JpegPic=bytarr(pixels,pixels,3)

        for i=0, pixels-1 do begin 
              for j=0, pixels-1 do begin
                 JPegPic(i,j,0)= ColPlane(i/multfac,j/multfac,0)
                 JPegPic(i,j,1)= ColPlane(i/multfac,j/multfac,1)
                 JPegPic(i,j,2)= ColPlane(i/multfac,j/multfac,2)
              endfor
        endfor




   write_jpeg, jpgfile , Jpegpic, true=3, quality=90


   cmd = 'mogrify -font arial -antialias -pointsize 58 -fill white -draw '
   cmd = cmd + "'"
   cmd = cmd + 'text 40,72 '
   cmd = cmd + """
   cmd = cmd + datfile
   cmd = cmd + """
   cmd = cmd + "'"
   cmd = cmd + ' -quality 95 '
   cmd = cmd + jpgfile
   print, cmd
   spawn, cmd


end





;======================================================================



pro test_2, junk

   ; Test
   ; Volker's Pink/White/Blue Colors
   ; ----------------------------------
	;filename='test_colplane3.jpg'  & testid= 'test3'   ; i, j, 100
	;filename='test_colplane4.jpg'  & testid= 'test4' ; i, j, 200
	;filename='test_colplane5.jpg'  & testid= 'test5'  ; i, j, i
	;filename='test_colplane6.jpg'  & testid= 'test6'  ; i, j, 255-i
	filename='test_colplane7.jpg'  & testid= 'test7'  ; complex function


        cols=256
        ColPlane=fltarr(cols,cols,3)

        for i=0, cols-1 do begin 
              for j=0, cols-1 do begin
                 ;ColPlane(i,j,0)= i
                 ;ColPlane(i,j,1)= j
                 ;ColPlane(i,j,2)= 255-i
                 ColPlane(i,j,0)= i
                 ColPlane(i,j,1)= (cols-1) * (1- cos(!PI/2 * i/cols))
                 ColPlane(i,j,2)= (j/cols) * (cols-1) * tan(!PI / 4. * i/cols)
              endfor
        endfor

	; 255, 255, 255 = white
	; 0, 0, 0 = black

	;ColPlane[0:40,0:10,*]= 255

	;multfac= 1
	;multfac= 2
	multfac= 3
        pixels= multfac*cols
        JpegPic=bytarr(pixels,pixels,3)

        for i=0, pixels-1 do begin 
              for j=0, pixels-1 do begin
                 JPegPic(i,j,0)= ColPlane(i/multfac,j/multfac,0)
                 JPegPic(i,j,1)= ColPlane(i/multfac,j/multfac,1)
                 JPegPic(i,j,2)= ColPlane(i/multfac,j/multfac,2)
                 ;JPegPic(i,j,0)= i/3
                 ;JPegPic(i,j,1)= j/3
                 ;JPegPic(i,j,2)= 200
              endfor
        endfor



   write_jpeg, filename , Jpegpic, true=3, quality=90


   cmd = 'mogrify -font arial -antialias -pointsize 58 -fill white -draw '
   cmd = cmd + "'"
   cmd = cmd + 'text 40,72 '
   cmd = cmd + """
   cmd = cmd + testid
   cmd = cmd + """
   cmd = cmd + "'"
   cmd = cmd + ' -quality 95 '
   cmd = cmd + filename
   print, cmd
   spawn, cmd


end
