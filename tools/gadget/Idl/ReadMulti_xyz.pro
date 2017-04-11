

Base= "/gpfs/mpa/vrs/zoom/test1/snapdir_016/"
SnapBase="snap"

Num = 16  ; number of snapshot files


exts='000'
exts=exts+strcompress(string(Num),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)

f= Base + "/"+ snapbase+"_"+exts +".0"
f= Strcompress(f, /remove_all)

falt = Base + "/"+ snapbase+"_"+exts 
falt= Strcompress(falt, /remove_all)


npart=lonarr(6)		
massarr=dblarr(6)
time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedback=0L
npartall=lonarr(6)	
flag_cooling= 0L
Nsubfiles = 0L
BoxSize = 0.0D

openr,1,f,/f77_unformatted, error = err
if err ne 0 then begin
  openr,1,falt,/f77_unformatted, error = err
endif

readu,1, npart, massarr, time, redshift, flag_sfr, flag_feedback, $
         npartall, flag_cooling, Nsubfiles, BoxSize 
close,1


print, npartall
print,time,redshift
print, boxsize


count1=0L
count2=0L
for subfile=0, Nsubfiles-1 do begin

    f=base + "/"+ snapbase+"_"+exts +"." + string(subfile)
    f=strcompress(f,/remove_all)
   
    if Nsubfiles le 1 then begin
       f=base + "/"+ snapbase+"_"+exts 
       f=strcompress(f,/remove_all)
    endif

    openr,1,f,/f77_unformatted
    readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall
    print,npart,massarr
    print
    N1= npart(1)
    N2= npart(2) + npart(3) + npart(4) + npart(5)
    pos1 = fltarr(3,N1)
    pos2 = fltarr(3,N2)
    readu,1,pos1,pos2
    close,1

    if subfile eq 0 then begin
       PosA= fltarr(3, npartall(1))
       PosB= fltarr(3, npartall(2) + npartall(3) + npartall(4) + npartall(5))
    endif

    PosA(*, count1:count1+N1-1)= Pos1(*,*)
    PosB(*, count2:count2+N2-1)= Pos2(*,*)
    count1=count1+N1
    count2=count2+N2
endfor

window, xsize=800, ysize=800, retain=2 
     
posA -= 250
posB -= 250




len = 3.0

ind=where(abs(posA(2,*)) lt len)
ind=where(abs(posB(2,*)) lt len)

plot, posB(0,ind), posB(1,ind), xrange=[-len,len], yrange=[-len,len],psym=3





end

















 
