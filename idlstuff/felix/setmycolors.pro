PRO setmycolors
  
   ; define my color scheme
   R      = bytarr(256)
   G      = bytarr(256)
   B      = bytarr(256)
  
   ; last color (background) 
   R(255) = 255
   G(255) = 255
   B(255) = 255
   
   ; 1 = red
   R(1)   = 255
   G(1)   = 0
   B(1)   = 0

   ; 2 = blue
   R(2)   = 0
   G(2)   = 0
   B(2)   = 255

   ; 3 = darkgreen
   R(3)   = 0
   G(3)   = 250
   B(3)   = 87

   ; 4 = very light red/orange
   R(4)   = 255
   G(4)   = 255
   B(4)   = 50

   ; 5 = very light blue
   R(5)   = 220
   G(5)   = 220
   B(5)   = 255

   ; 6 = bright blue
   R(6)   = 90
   G(6)   = 120
   B(6)   = 255

   ; 7 = light green
   R(7)   = 0
   G(7)   = 160 ; 250
   B(7)   = 0  ; 87

   ; 8 = lila
   R(8)   = 175 ;148
   G(8)   = 0 ;0
   B(8)   = 230;211

   ; 9 = turquoise
   ;R(9)   = 0
   ;G(9)   = 220
   ;B(9)   = 220

   R(9)   = 139
   G(9)   = 35
   B(9)   = 35

   ; 10 = light red

   R(10)   = 255
   G(10)   = 200; 240
   B(10)   = 200; 240

   tvlct,R,G,B   

END
