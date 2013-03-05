
print, "running: volker_showprojdenfield"
;Base = "../../"
Base = "./"

;Num = 0
Num = 200


loadct, 15
tvlct, r, g, b, /get

;openr,1,"rainbow.clt"
openr,1,"/n/home/tcox/idlstuff/rainbow.clt"
readf,1,r
readf,1,g
readf,1,b
close,1



if num ge 1000 then begin
    exts='0000'
    exts=exts+strcompress(string(Num),/remove_all)
    exts=strmid(exts,strlen(exts)-4,4)
endif else begin
    exts='000'
    exts=exts+strcompress(string(Num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)
endelse

f= Base + "/proj_density_field_"+exts 
f= Strcompress(f, /remove_all)
XPixels = 0L
YPixels = 0L

openr,1,f
readu,1,XPixels
readu,1,YPixels
dens= fltarr(YPixels, XPixels)
readu,1,dens
temp= fltarr(YPixels, XPixels)
readu,1,temp
close,1

dens= transpose(dens)
temp= transpose(temp)

mi = min(dens)
ma = max(dens)


;;; Now need to map the result onto a color table

Map= Dens     ; select mass projection


ma= max(Map)
mi= ma/1000.0


;; now do clipping in Map

ind=where(Map lt mi) 
if ind(0) ne -1 then Map(ind)=mi
ind=where(Map gt ma)
if ind(0) ne -1 then Map(ind)=ma

;; now map the log of the field linearly onto the color table

cols= 255

Map= (alog10(Map)-alog10(mi))/(alog10(ma/mi)) * (cols-2) +2


MapT= Temp

maT= max(MapT)
miT= maT/100.0




ind=where(MapT lt miT) 
if ind(0) ne -1 then MapT(ind)=miT
ind=where(MapT gt maT)
if ind(0) ne -1 then MapT(ind)=maT

;; now map the log of the field linearly onto the color table

cols= 255

MapT= (alog10(MapT)-alog10(miT))/(alog10(maT/miT)) * (cols-2) +2







;loadct, 15                       ;  load color table
;tvlct,red,green, blue, /get


cols=256
ColPlane=fltarr(cols,cols,3)

;openr,1,"colplane2.dat"
openr,1,"/n/home/tcox/idlstuff/colplane2.dat"
readu,1,ColPlane
close,1

Pic=bytarr(XPixels,YPixels,3)

for i=0, XPixels-1 do begin
   for j=0, YPixels-1 do begin
      Pic(i,j,0)= ColPlane(Map(i,j),MapT(i,j),0)
      Pic(i,j,1)= ColPlane(Map(i,j),MapT(i,j),1)
      Pic(i,j,2)= ColPlane(Map(i,j),MapT(i,j),2)
   endfor
endfor





  
exts='0000'
exts=exts+strcompress(string(Num),/remove_all)
exts=strmid(exts,strlen(exts)-4,4)

fname =Base+"/pic_"+exts+".jpg"
write_jpeg, fname, pic, true=3, quality=98

tv, pic, true=3


end





 









