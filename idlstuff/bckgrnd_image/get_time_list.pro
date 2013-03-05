pro get_time_list, path, snapbase, NumList, TiList


NumList=[0]
TiList=[0]


num=0

repeatit:

exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)
f=path + snapbase+ "_"+exts
f=strcompress(f,/remove_all)

npart=lonarr(6)	
massarr=dblarr(6)
time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedback=0L
npartall=lonarr(6)	
bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4
la=intarr(bytesleft/2)

openr,1,f,/f77_unformatted, error=err ;, /swap_endian

if err ne 0 then goto,skip

readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,la
close,1

NumList=[NumList,Num]
TiList=[TiList,time]

num=num+1
goto,repeatit


skip:

if n_elements(NumList) gt 1 then begin
    NumList=NumList(1:*)
    TiList=TiList(1:*)
endif else begin
    print,"No files found"
endelse

end
