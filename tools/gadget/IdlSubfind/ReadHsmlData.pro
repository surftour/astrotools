
Num = 5  

Base= "/ptmp/vrs/Billennium/res0/"

sfrflag = 0
bhflag = 0

;--------------------------------------------------

exts='000'
exts=exts+strcompress(string(Num),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)


skip = 0L
fnr = 0L

repeat begin

    f = Base + "/hsmldir_" + exts +"/hsml_"+exts +"."+strcompress(string(fnr),/remove_all)

    Nhsml = 0L
    Nprevious = 0L
    Ntotal = 0LL
    Ntask = 0L

    openr,1,f
    readu,1, Nhsml, Nprevious, Ntotal, NTask

        print, Nhsml, Nprevious, Ntotal

    if fnr eq 0 then begin
        Hsml = fltarr(Ntotal)
        Density = fltarr(Ntotal)
        Veldisp = fltarr(Ntotal)
    endif

    if Nhsml gt 0 then begin
        
        locHsml = fltarr(Nhsml)
        locDensity = fltarr(Nhsml)
        locVeldisp = fltarr(Nhsml)

        readu,1, locHsml
        readu,1, locDensity
        readu,1, locVelDisp

        Hsml(skip:skip+Nhsml-1) = locHsml(*)
        Density(skip:skip+Nhsml-1) = locDensity(*)
        Veldisp(skip:skip+Nhsml-1) = locVeldisp(*)
        skip+= Nhsml

    endif

    close, 1

    fnr++

endrep until fnr eq NTask

print
print, "Total particles in density lists =", NTotal
print



;plot, hsml, density, psym=3,/xlog, /ylog
;plot, hsml, veldisp, psym=3,/xlog, /ylog

end
