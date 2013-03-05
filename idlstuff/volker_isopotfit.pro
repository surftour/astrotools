;Object_File= '/home/volker/tmp/ComputeIsoPhotFitIDL/isophotfit.so'
Object_File= '/home/tcox/Tools/C-Routines_for_IDL/ComputeIsoPhotFitIDL/isophotfit.so'
Base=     "/home/tcox/"
Num= 6
fout = "isophot.eps"

NLev = 10


gr= 1
ex=strcompress(string(gr),/remove_all)
exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)


fname=base+"potprojdm_"+exts+"_"+ex+".dat"
openr,1,fname
Pixels=0L 
readu,1,Pixels
M200=0.0
R200=0.0
Size=0.0
readu,1,R200,M200,Size
PotMapZ=       fltarr(Pixels,Pixels)
PotMapX=       fltarr(Pixels,Pixels)
PotMapY=       fltarr(Pixels,Pixels)
readu,1,PotMapX
readu,1,PotMapY
readu,1,PotMapZ
close,1


p1= Pixels / (2*size) * (size-R200)
p2= Pixels / (2*size) * (size+r200)

potmean = (PotMapX(p1, PIxels/2) +  PotMapX(p2, PIxels/2) + $
           PotMapY(PIxels/2, p1) +  PotMapY(PIxels/2, p2) + $
           PotMapZ(p1, PIxels/2) +  PotMapZ(p2, PIxels/2)) /6

PotMapX = PotMapX - potmean
PotMapY = PotMapY - potmean
PotMapZ = PotMapZ - potmean

Pot = transpose(PotMapZ)
MinPot= min(Pot)
PotInt = -MinPot/NLev

Map= -Pot

ma= max(Map)+PotInt
mi= 0 - 2*PotInt

ind=where(Map gt ma)
if ind(0) ne -1 then begin
    Map(ind)=ma
endif
ind=where(Map lt mi)
if ind(0) ne -1 then begin
    Map(ind)=mi
endif
image= byte((Map-mi)/(ma-mi)* (!D.table_size-3))+2




set_plot,"PS"
device,filename= fout,  /encapsulated,  /color ,bits_per_pixel=8
device,xsize=16,ysize=16
!p.font=0
device,/times,/italic,font_index=20
loadct, 26
tvlct,r,g,b,/get

v1=[255,  0,  0,255,0]
v2=[  0,255,  0,000,255]
v3=[  0,  0,255,255,255]

tvlct,v1,v2,v3,1

!p.position=[0.03,0.03,0.97,0.97]
!p.ticklen=0.025


plot,[-1],[0],psym=3,xrange=[0,PIXELS-1],yrange=[0,PIXELS-1],xstyle=1+4,ystyle=1+4, $
  /noerase

a0=!x
b0=!y


plot,[0],[0],psym=3,xrange=[-size,size],yrange=[-size,size], $
  xstyle=1,ystyle=1+4, /noerase, xthick=2.5, ythick=2.5, $
  xtitle = "!20X!7  [ !20h!7!U-1!Nkpc ]!3",$
  ytitle = "!20Y!7  [ !20h!7!U-1!Nkpc ]!3"


axis,yaxis=1, xrange=[-size,size],yrange=[-size,size], $
  xstyle=1,ystyle=1, /noerase, xthick=2.5, ythick=2.5, $
  xtitle = "!20X!7  [ !20h!7!U-1!Nkpc ]!3",$
  ytitle = "!20Y!7  [ !20h!7!U-1!Nkpc ]!3"


axis,yaxis=0, xrange=[-size,size],yrange=[-size,size], $
  xstyle=1,ystyle=1, /noerase, xthick=2.5, ythick=2.5, $
  xtitle = "!20X!7  [ !20h!7!U-1!Nkpc ]!3",$
  ytitle = "!20Y!7  [ !20h!7!U-1!Nkpc ]!3", charsize=0.001

a1=!x
b1=!y

lev= (indgen(NLEV) + 1 )*PotInt + MinPot

!x=a0
!y=b0

contour,  Pot, /overplot, /follow, LEVELS=lev, thick=4.0, $
           c_annotation=['','','','','','','','','',''] ,c_charsize=0.0001


for i=0,n_elements(lev)-1 do begin


    PIXELS=long(PIXELS)
    V=-float(Pot)
    Scale=float(2*size)
    Thresh= -float(lev(i))

    Xfit=fltarr(360)
    YFit=fltarr(360)

    Xcon=fltarr(360)
    Ycon=fltarr(360)

    AFit=0.0
    BFit=0.0
    X0Fit=0.0
    Y0Fit=0.0
    PhiFit=0.0
    a4Fit=0.0


    S = CALL_EXTERNAL(object_file, $
                      'fit_isophot', $
                      PIXELS, $
                      V,$
                      Thresh,$
                      Scale, $
                      Xfit, Yfit,$
                      PhiFit, AFit, BFit, X0Fit, Y0Fit, a4Fit, $
                      Xcon, Ycon)
    !x=a1
    !y=b1

    oplot, [Xfit, xfit(0)], [Yfit, yfit(0)], color=1

    print, "a4/a=", a4Fit/AFit

endfor

device,/close
set_plot,"X"

end



