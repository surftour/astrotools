;------------------------------------------------------------------------
;
;    Random Procedures related to Hot Gas
;
;
;
;
;------------------------------------------------------------------------










;------------------------------------------------------------------------
;
;     Histogram of Gas Cooling Time
;     ---------------------------------------------
;
;   plots a histogram of gas cooling times, calculated from the tabulated
;   cooling curves of sutherland and dopita
;
;------------------------------------------------------------------------

pro gas_coolingtime, frun, snapnum, filename=filename


if not keyword_set(frun) then begin
   print, "  "
   print, "gas_coolingtime, frun, snapnum, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='coolt.eps'
if not keyword_set(snapnum) then snapnum=25


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename 


;--------------------------------------
;--------------------------------------

;xmax = 100
;xmin = 0
xmax= 3.0
xmin= -1

ymax = 1000
ymin = 0


levels= 100.0


;if not keyword_set(r_g) then begin
;        rho= fload_gas_rho(1)
;        u= fload_gas_u(1)
;        utp= fload_gas_utp(1)
;endif



; load snapshot
ok=fload_snapshot_bh(frun,snapnum)


; ------------------------
;  Calculate Cooling Time
; ------------------------
m_p= 1.6726d-24         ; proton mass (in cgs, i.e. g)
k_b= 1.3806e-16         ; boltzmann constant (in cgs)
t_cgs= 3.08568e+16          ; convert time to seconds (cgs units)
UnitEnergy_in_cgs=        1.989d+53
UnitMass_in_g    =        1.989d+43
UnitDensity_in_cgs =      6.7699112d-22


; load gas quantities
; --------------------
u= fload_gas_u(1)
rho= fload_gas_rho(1)


; determine cooling rate from sd files
; -------------------------------------
;ok=load_cooling_table('zero')
ok=fload_cooling_table('-10')           ; load cooling file
temp= u*fload_gas_mu(1)/0.012381322    ; temp in kelvin
coolrate= fload_lambdan_array(temp)    ; cooling rate (ergs cm3 s-1)
zeroidx= where(coolrate le 0)
if zeroidx(0) ne -1 then coolrate(zeroidx)= 1.0


; convert to physical units
; --------------------------
rho_cgs= rho*UnitDensity_in_cgs         ; g/cm3
u_cgs= u*UnitEnergy_in_cgs/UnitMass_in_g  ; in ergs/g per unit volume


cooltime= u_cgs*m_p*m_p/(rho_cgs*coolrate)     ; cooling time (s)
cooltime= cooltime/3.08568d+16         ; cooling time (Gadget Units - Gyr)
if zeroidx(0) ne -1 then cooltime(zeroidx)= 0.0



; ------------------------
; universal step size
; ------------------------
step= (xmax-xmin)/levels


; ---------------------
;  do statistics
; ---------------------
; actual values
idx= where(cooltime le 100)
if idx(0) ne -1 then cooltime_trimmed= cooltime(idx) else cooltime_trimmed= cooltime
cooltime_moment= moment(cooltime_trimmed)
sigma= sqrt(cooltime_moment[1])
ct_avg= cooltime_moment[0]



; ------------------
; compute histogram
; ------------------

xaxislog= 1
if xaxislog eq 1 then begin
        xaxistitle= 'Log Cooling Time (Gyr)'

        if zeroidx(0) ne -1 then cooltime(zeroidx)= 1.0
        cooltime_log= alog10(cooltime)
        if zeroidx(0) ne -1 then cooltime_log(zeroidx)= -5
        idx= where(cooltime_log GT xmax)
        if idx(0) NE -1 then cooltime_log(idx)= xmax
        idx= where(cooltime_log LT xmin)
        if idx(0) NE -1 then cooltime_log(idx)= xmin

        cooltimetohist= cooltime_log

        ; now take avg of log values
        ;dx= where(cooltime le 100)
        ;if idx(0) ne -1 then cooltime_trimmed= cooltime(idx) else cooltime_trimmed= cooltime
        cooltime_moment= moment(cooltime_log)
        sigma= sqrt(cooltime_moment[1])
        ct_avg= 10^(cooltime_moment[0])
endif else begin
        xaxistitle= 'Cooling Time (Gyr)'

        idx= where(cooltime GT xmax)
        if idx(0) NE -1 then cooltime(idx)= xmax
        idx= where(cooltime LT xmin)
        if idx(0) NE -1 then cooltime(idx)= xmin

        cooltimetohist= cooltime
endelse



hist_ct= histogram(cooltimetohist, binsize=step, max=xmax, min=xmin)



;  make it a log plot
; -----------------------
ytit= 'Number'
makeitlog= 0
if makeitlog EQ 1 then begin
        ymin= -1
        ymax= 3
        ytit= 'Log Number'
        idx= where(hist_ct le 0) & hist_ct(idx)= 1e6
        hist_ct= alog10(hist_ct)
        hist_ct(idx)= -1
endif


; -----------
; x axis
; -----------
bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
bins=[bins,xmax]



        !p.position= [0.20, 0.15, 0.95, 0.95]


        plot, bins, hist_ct, psym=10, $
                xstyle=1,ystyle=1, $
                xrange=[xmin,xmax],$
                yrange=[ymin,ymax],$
                color=0, $
                xcharsize=1.50, ycharsize=1.50, $
                xthick=4.0, ythick=4.0, $
                charthick=2.0, $
                xtitle=xaxistitle, $
                ytitle=ytit





        xyouts, 0.22, 0.90, /normal, fload_timelbl(1), size=1.0, color=0
        xyouts, 0.65, 0.90, /normal, fload_fid(1), size=1.33, color=0

        ;print average velocity
        ctavglbl='t!Dcool!N='+strcompress(string(ct_avg),/remove_all)
        ctavglbl=strmid(ctavglbl,0,14)+' Gyr'
        xyouts, 0.65, 0.85, ctavglbl,/normal,size= 1.33, color=0

        ;print sigma
        siglbl='!7r!3!N='+strcompress(string(sigma),/remove_all)
        siglbl=strmid(siglbl,0,14)
        xyouts, 0.65, 0.80, siglbl,/normal,size= 1.33, color=0




device, /close



end


