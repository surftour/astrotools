;-------------------------------------------------------------------



pro load_stellar_magnitudes, Nstars, Time, Mass_Stars, Age_Stars, Z_Stars, $
				BolMag, UMag, BMag, VMag, KMag

;  Use brant's code to look up Bruzual & Charlot 2003
; stellar SED's to determine luminosities 
; of new (and eventually old) stellar particles
;

; Stellar Variables are already passed


; ---------------------------------------
; Actually determine the luminosities
; ---------------------------------------
if(Nstars gt 0) then begin

        NewStarLuminosities= fltarr(14,Nstars)      ;luminosities

        S = CALL_EXTERNAL('/n/home03/tcox/Tools/C-Routines_for_IDL/ComputeColors/colors', $
                'colors', $
                Nstars, $
                Time, $
                Mass_Stars, $
                Age_Stars, $
                Z_Stars, $
                NewStarLuminosities)
endif else begin
        print,'No stars, no luminosities (or magnitudes, in this case)!'
	return
endelse

; returns magnitudes
;  1. bolometric
;  2. U
;  3. B
;  4. V
;  5. R
;  6. I
;  7. J
;  8. H
;  9. K
; 10. u
; 11. g
; 12. r
; 13. i
; 14. z


BolMag= NewStarLuminosities(0,*)
UMag= NewStarLuminosities(1,*)
BMag= NewStarLuminosities(2,*)
VMag= NewStarLuminosities(3,*)
KMag= NewStarLuminosities(9,*)



end




;-------------------------------------------------------------------



;=====================================
; luminosity version of the above
;=====================================
pro load_stellar_luminosities, Nstars, Time, Mass_Stars, Age_Stars, Z_Stars, $
                                Lum_Bol, Lum_u, Lum_b, Lum_v, Lum_k, notsolar=notsolar



load_stellar_magnitudes, Nstars, Time, Mass_Stars, Age_Stars, Z_Stars, $
			BolMag, UMag, BMag, Vmag, KMag


; -----------------------
;  Convert to Luminosity
; -----------------------
if keyword_set(notsolar) then begin
	Lum_Bol= (10^(-0.4*BolMag)) * Mass_Stars
	Lum_u= (10^(-0.4*UMag)) * Mass_Stars
	Lum_b= (10^(-0.4*BMag)) * Mass_Stars
	Lum_v= (10^(-0.4*VMag)) * Mass_Stars
	Lum_k= (10^(-0.4*Kmag)) * Mass_Stars
endif else begin
	Lum_Bol= (10^(-0.4*BolMag)) * Mass_Stars
	Lum_u= (10^(-0.4*(UMag-5.60))) * Mass_Stars
	Lum_b= (10^(-0.4*(BMag-5.51))) * Mass_Stars
	Lum_v= (10^(-0.4*(VMag-4.84))) * Mass_Stars
	Lum_k= (10^(-0.4*(KMag-3.33))) * Mass_Stars
endelse


end


;-------------------------------------------------------------------



;=======================================
;
; The above was just for 5 fields, below
; we do it for all the bands.
;
;=======================================
pro load_all_stellar_magnitudes, Nstars, Time, Mass_Stars, Age_Stars, Z_Stars, $
				BolMag, $
				UMag, $
				BMag, $
				VMag, $
				RMag, $
				IMag, $
				JMag, $
				HMag, $
				KMag, $
				usdssMag, $
				gsdssMag, $
				rsdssMag, $
				isdssMag, $
				zsdssMag

;  Use brant's code to look up Bruzual & Charlot 2003
; stellar SED's to determine luminosities 
; of new (and eventually old) stellar particles
;

; Stellar Variables are already passed


; ---------------------------------------
; Actually determine the luminosities
; ---------------------------------------
print, "N_stars= ",Nstars

if(Nstars gt 0) then begin

        NewStarLuminosities= fltarr(14,Nstars)      ;luminosities
        MassInfo= fltarr(2,Nstars)      ;luminosities

;  old version of colors code
;------------------------------
;        S = CALL_EXTERNAL('/home03/tcox/Tools/C-Routines_for_IDL/ComputeColors/colors', $
;                'colors', $
;                Nstars, $
;                Time, $
;                Mass_Stars, $
;                Age_Stars, $
;                Z_Stars, $
;                NewStarLuminosities)


;  new version of colors code (see below)
;------------------------------
;model= 0   ; salpeter
model= 1   ; chabrier
silent= 0
;silent= 1
;S = CALL_EXTERNAL('/home03/brant/code/idl/colors/colors', $
S = CALL_EXTERNAL('/n/home03/tcox/Tools/C-Routines_for_IDL/colors/colors', $
        'main', $
        long(Nstars), $
        Age_Stars, $           ;age lin gyr
        Z_Stars, $   ;metallicity in solar units
        NewStarLuminosities, $  ;luminosities in abs magnitudes
	MassInfo, $
;        long(model))
        long(model), $   ;model==0 -> Salpeter, model==1 Chabrier
        long(silent), $   ; produce output - or not
        /F_VALUE)



endif else begin
        print,'No stars, no luminosities (or magnitudes, in this case)!'
	return
endelse


;
;  Brant's newest version has different order
;
BolMag= NewStarLuminosities(0,*)
;UMag= NewStarLuminosities(1,*)
;BMag= NewStarLuminosities(2,*)
;VMag= NewStarLuminosities(3,*)
;RMag= NewStarLuminosities(4,*)
;IMag= NewStarLuminosities(5,*)
;JMag= NewStarLuminosities(6,*)
;HMag= NewStarLuminosities(7,*)
;KMag= NewStarLuminosities(8,*)
;usdssMag= NewStarLuminosities(9,*)
;gsdssMag= NewStarLuminosities(10,*)
;rsdssMag= NewStarLuminosities(11,*)
;isdssMag= NewStarLuminosities(12,*)
;zsdssMag= NewStarLuminosities(13,*)
;
;
;
usdssMag= NewStarLuminosities(1,*)
gsdssMag= NewStarLuminosities(2,*)
rsdssMag= NewStarLuminosities(3,*)
isdssMag= NewStarLuminosities(4,*)
zsdssMag= NewStarLuminosities(5,*)
UMag= NewStarLuminosities(6,*)
BMag= NewStarLuminosities(7,*)
VMag= NewStarLuminosities(8,*)
RMag= NewStarLuminosities(9,*)
IMag= NewStarLuminosities(10,*)
JMag= NewStarLuminosities(11,*)
HMag= NewStarLuminosities(12,*)
KMag= NewStarLuminosities(13,*)




end


;-------------------------------------------------------------------


;=====================================
; luminosity version of the above
;=====================================
pro load_all_stellar_luminosities, Nstars, Time, Mass_Stars, Age_Stars, Z_Stars, $
                                Lum_Bol, $
				Lum_U, $
				Lum_B, $
				Lum_V, $
				Lum_R, $
				Lum_I, $
				Lum_J, $
				Lum_H, $
				Lum_K, $
				Lum_sdss_u, $
				Lum_sdss_g, $
				Lum_sdss_r, $
				Lum_sdss_i, $
				Lum_sdss_z, $
				notsolar=notsolar



load_all_stellar_magnitudes, Nstars, Time, Mass_Stars, Age_Stars, Z_Stars, $
			BolMag, UMag, BMag, VMag, RMag, IMag, JMag, HMag, KMag, $
			usdssMag, gsdssMag, rsdssMag, isdssMag, zsdssMag



; -----------------------
;  Convert to Luminosity
; -----------------------
if keyword_set(notsolar) then begin
	Lum_Bol= (10^(-0.4*BolMag)) * Mass_Stars
	Lum_U= (10^(-0.4*UMag)) * Mass_Stars
	Lum_B= (10^(-0.4*BMag)) * Mass_Stars
	Lum_V= (10^(-0.4*VMag)) * Mass_Stars
	Lum_R= (10^(-0.4*RMag)) * Mass_Stars
	Lum_I= (10^(-0.4*IMag)) * Mass_Stars
	Lum_J= (10^(-0.4*JMag)) * Mass_Stars
	Lum_H= (10^(-0.4*HMag)) * Mass_Stars
	Lum_K= (10^(-0.4*KMag)) * Mass_Stars
	Lum_sdss_u= (10^(-0.4*usdssMag)) * Mass_Stars
	Lum_sdss_g= (10^(-0.4*gsdssMag)) * Mass_Stars
	Lum_sdss_r= (10^(-0.4*rsdssMag)) * Mass_Stars
	Lum_sdss_i= (10^(-0.4*isdssMag)) * Mass_Stars
	Lum_sdss_z= (10^(-0.4*zsdssMag)) * Mass_Stars
endif else begin
	;Lum_Bol= (10^(-0.4*(BolMag-4.75))) * Mass_Stars
	;Lum_U= (10^(-0.4*(UMag-5.60))) * Mass_Stars
	;Lum_B= (10^(-0.4*(BMag-5.51))) * Mass_Stars
	;Lum_V= (10^(-0.4*(VMag-4.84))) * Mass_Stars
	;Lum_R= (10^(-0.4*(RMag-4.48))) * Mass_Stars
	;Lum_I= (10^(-0.4*(IMag-4.13))) * Mass_Stars
	;Lum_J= (10^(-0.4*(JMag-3.70))) * Mass_Stars
	;Lum_H= (10^(-0.4*(HMag-3.37))) * Mass_Stars
	;Lum_K= (10^(-0.4*(KMag-3.33))) * Mass_Stars
	;Lum_sdss_u= (10^(-0.4*(usdssMag-6.28))) * Mass_Stars
	;Lum_sdss_g= (10^(-0.4*(gsdssMag-4.95))) * Mass_Stars
	;Lum_sdss_r= (10^(-0.4*(rsdssMag-4.45))) * Mass_Stars
	;Lum_sdss_i= (10^(-0.4*(isdssMag-4.35))) * Mass_Stars
	;Lum_sdss_z= (10^(-0.4*(zsdssMag-4.36))) * Mass_Stars

	; brant's new values - see below
	Lum_Bol= (10^(-0.4*(BolMag-4.74))) * Mass_Stars
	Lum_U= (10^(-0.4*(UMag-5.56))) * Mass_Stars
	Lum_B= (10^(-0.4*(BMag-5.45))) * Mass_Stars
	Lum_V= (10^(-0.4*(VMag-4.80))) * Mass_Stars
	Lum_R= (10^(-0.4*(RMag-4.46))) * Mass_Stars
	Lum_I= (10^(-0.4*(IMag-4.10))) * Mass_Stars
	Lum_J= (10^(-0.4*(JMag-3.66))) * Mass_Stars
	Lum_H= (10^(-0.4*(HMag-3.32))) * Mass_Stars
	Lum_K= (10^(-0.4*(KMag-3.28))) * Mass_Stars
	Lum_sdss_u= (10^(-0.4*(usdssMag-6.75))) * Mass_Stars
	Lum_sdss_g= (10^(-0.4*(gsdssMag-5.33))) * Mass_Stars
	Lum_sdss_r= (10^(-0.4*(rsdssMag-4.67))) * Mass_Stars
	Lum_sdss_i= (10^(-0.4*(isdssMag-4.48))) * Mass_Stars
	Lum_sdss_z= (10^(-0.4*(zsdssMag-4.42))) * Mass_Stars
endelse


end







;-------------------------------------------------------------------


;========================================
;   Test the luminosities
;==========================================
pro print_luminosity_info, frun, snapnum


;  Use brant's code to look up Bruzual & Charlot 2003
; stellar SED's to determine luminosities 
; of new (and eventually old) stellar particles
;
;
; this is mainly for test purposes
;

if not keyword_set(snapnum) then snapnum= 20

if not keyword_set(frun) then begin
   print, "  "
   print, "print_luminosity_info, frun, snapnum"
   print, "  "
   print, "  "
   print, "  "
   return
endif


;--------------------------------------
;  Set up Parameters
;--------------------------------------

loadedsnap= 0

    if not keyword_set(loadedsnap) then begin
      ;if (fload_snapshot(frun, snapnum)) then begin
      if (fload_snapshot_bh(frun, snapnum)) then begin
        print, "PROBLEM: opening file"
        return
      endif
    endif


Ttime= float(fload_time(1))
print, "Time= ", Ttime



; new stars
; ------------
Nstars= long(fload_npart(4))
Mass_Stars= 1.0e+10*fload_newstars_mass(1)   ; Brant's code needs mass in solar masses
Age_Stars= fload_newstars_age(1)
Z_Stars= fload_newstars_z(1)

; now convert to true age!
Age_Stars=float(Ttime-Age_Stars)


; get the luminosities
;  - in units of solar luminosities
print, "load luminosities"
load_all_stellar_luminosities, Nstars, Ttime, Mass_Stars, Age_Stars, Z_Stars, $
        Lum_Bol, Lum_U, Lum_B, Lum_V, Lum_R, Lum_I, Lum_J, Lum_H, Lum_K, $
        Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z   ;, /notsolar


NewStarBoloLum= total(Lum_Bol)
print, "New Stars Bolometric Luminosity= ", NewStarBoloLum



; all stars
; ------------
Nstars= long(fload_npart(2)+fload_npart(3)+fload_npart(4))
Mass_Stars= 1.0e+10*fload_allstars_mass(1)
Age_Stars= fload_allstars_age(1)
Z_Stars= fload_allstars_z(1)

; now convert to true age!
Age_Stars=float(Ttime-Age_Stars)


; get the luminosities
;  - in units of solar luminosities
print, "load luminosities"
load_all_stellar_luminosities, Nstars, Ttime, Mass_Stars, Age_Stars, Z_Stars, $
	Lum_Bol, Lum_U, Lum_B, Lum_V, Lum_R, Lum_I, Lum_J, Lum_H, Lum_K, $
	Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z   ;, /notsolar



AllStarBoloLum= total(Lum_Bol)
print, "All Stars Bolometric Luminosity= ", AllStarBoloLum

print, "percent from new stars= ", 100.0*NewStarBoloLum/AllStarBoloLum



; -------------------
;  Print some stuff
; -------------------

;print_age_u, Age_Stars, Lum_u, filename='ageu.eps', $
;		frun=frun, snapnum=snapnum, msg='Brant'


;print_z_age, Z_Stars, Age_Stars, filename='zage.eps'


; -------------
;  Done
; -------------


end




;-------------------------------------------------------------------




;========================================
;  This code got combined with Brant's
; into the routine above.  Hence this 
; routine is commented out and no longer
; used.
;========================================

;pro BC_lums, frun, snapnum


;;  Use Paul Martini's code to look up 
;; Bruzual & Charlot 2003 SED's for specific
;; age and metallicity stellar populations.
;;

;if not keyword_set(snapnum) then snapnum= 20


;if not keyword_set(frun) then begin
;   print, "  "
;   print, "BC_luminosities, frun, snapnum"
;   print, "  "
;   print, "  "
;   print, "  "
;   return
;endif





;;--------------------------------------
;;  Set up Parameters
;;--------------------------------------

;loadedsnap= 0

;    if not keyword_set(loadedsnap) then begin
;      ;if (fload_snapshot(frun, snapnum)) then begin
;      if (fload_snapshot_bh(frun, snapnum)) then begin
;        print, "PROBLEM: opening file"
;        return
;      endif
;    endif



;Nstars= 100L
;Nstars= fload_npart(4)
;Time= fload_time(1)
;Mass_Stars= 1.0e+10*fload_newstars_mass(1)   ; Brant's code needs mass in solar masses
;Age_Stars= fload_newstars_age(1)
;Z_Stars= fload_newstars_z(1)


;Time= 0.0  ; for now, we'll try this


;; ---------------------------------------
;; Actually determine the luminosities
;; ---------------------------------------
;if(Nstars gt 0) then begin

;        NewStarLuminosities= fltarr(5,Nstars)      ;luminosities

;        S = CALL_EXTERNAL('/home03/tcox/Tools/C-Routines_for_IDL/ComputeBCMags/getBCmags.so', $
;                'determine_magnitudes', $
;                Nstars, $
;                Time, $
;                Mass_Stars, $
;                Age_Stars, $
;                Z_Stars, $
;                NewStarLuminosities)
;endif else begin
;        print,'No stars, no luminosities!'
;endelse


;; returns magnitudes
;;  1. bolometric
;;  2. u
;;  3. b
;;  4. v
;;  5. k
;;


;; -----------------------
;;  Convert to Luminosity
;; -----------------------
;Lum_Bol= (10^(-0.4*NewStarLuminosities(0,*))) * Mass_Stars
;Lum_u= (10^(-0.4*(NewStarLuminosities(1,*)-5.60))) * Mass_Stars
;Lum_b= (10^(-0.4*(NewStarLuminosities(2,*)-5.51))) * Mass_Stars
;Lum_v= (10^(-0.4*(NewStarLuminosities(3,*)-4.84))) * Mass_Stars
;Lum_k= (10^(-0.4*(NewStarLuminosities(4,*)-3.33))) * Mass_Stars


;; -------------------
;;  Print some stuff
;; -------------------

;print_age_u, Age_Stars, Lum_u, filename='ageu.eps', $
;		frun=frun, snapnum=snapnum, msg='Paul'

;print_z_age, Z_Stars, Age_Stars, filename='zage.eps'

;; -------------
;;  Done
;; -------------



;end










;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Plot u-band luminosities for fun
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





pro print_age_u, age, u, filename=filename, $
		frun=frun, snapnum=snapnum, msg=msg


if not keyword_set(msg) then msg=''

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

xaxistitle = 'Age (Gyr)'
yaxistitle = 'L!Lu!N (L!D!9n!3!N)'

xmax = max(age)
xmin = 0
;ymax = max(u)
;ymin = min(u)
ymax = 6.0e+8
ymin = 5.0e+4



;---------------------------


!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
	/ylog, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


oplot, age, u, thick=3.0, psym=1, color= 150

xyouts, 0.7, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0
xyouts, 0.7, 0.85, snapnum, /normal, charthick=1, size=1.33, color=0

if msg NE '' then xyouts, 0.7, 0.75, msg, /normal, charthick=3.0, size=1.33, color=0


device, /close


end










;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Plot metallicity and age
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





pro print_z_age, z, age, msg=msg, filename=filename


if not keyword_set(msg) then msg=''

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename

yaxistitle = 'Age (Gyr)'
xaxistitle = 'Z'

ymax = max(age)
ymin = 0
xmax = max(z)
xmin = min(z)



;---------------------------


!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
	/xlog, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


oplot, z, age, thick=3.0, psym=1, color= 150

;xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0

if msg NE '' then xyouts, 0.25, 0.80, msg, /normal, charthick=3.0, size=1.33, color=0


device, /close


end











;======================================================================================
;
;     June 13th, 2006.
;
;     I will now update my routines to use Brant's newest version of his colors
;     code. Below is is sample script which, and below that is his email about
;     the details.
;
;
;
;======================================================================================

function colors,N,mass,age,metallicity,model


luminosities = fltarr(14,N)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;	colors.pro
;
;	TAKES:
;	
;	mass        in Solar Masses (no h)
;	age         in Gyr          (no h)
;	metallicity in Solar units  (Z_sun = 0.02)
;	model:
;		0 == Salpeter IMF, 0.1, 100.0
;		1 == Chabrier IMF, 0.1, 100.0
;		default = Salpeter
;
;
;	USES:
;		Padova 1994 stellar tracks
;
;
;	RETURNS:
;
;		Luminosities in Solar Luminosities (no h)
;		NOTE: returns solar luminosities in the chosen
;		band!
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;==========================
;
; original version had these flipped
;
;
;Luminosity index legend
; 0 = bolometric luminosity
; 1 = Johnsons U
; 2 = Johnsons B
; 3 = Johnsons V
; 4 = Johnsons R
; 5 = Johnsons I
; 6 = Cousins J
; 7 = Cousins H
; 8 = Cousins K
; 9 = Sloan u
;10 = Sloan g
;11 = Sloan r
;12 = Sloan i
;13 = Sloan z
;
;==========================
;
;  correct version
;
;
;Luminosity index legend
; 0 = bolometric luminosity
; 1 = Sloan u
; 2 = Sloan g
; 3 = Sloan r
; 4 = Sloan i
; 5 = Sloan z
; 6 = Johnsons U
; 7 = Johnsons B
; 8 = Johnsons V
; 9 = Johnsons R
;10 = Johnsons I
;11 = Cousins J
;12 = Cousins H
;13 = Cousins K


;VEGA system
;from www.ucolick.org/~cnaw/sun.html
;following Fukugita et al. 1995, PASP, 105, 945
s_UBVRIJHK = fltarr(8)
s_UBVRIJHK(0) = 5.56;  //U (BESSEL)
s_UBVRIJHK(1) = 5.45;  //B (BESSEL)
s_UBVRIJHK(2) = 4.80;  //V (BESSEL)
s_UBVRIJHK(3) = 4.46;  //R (KPNO)
s_UBVRIJHK(4) = 4.10;  //I (KPNO)
s_UBVRIJHK(5) = 3.66;  //J (BESSEL)
s_UBVRIJHK(6) = 3.32;  //H (BESSEL)
s_UBVRIJHK(7) = 3.28;  //K (BESSEL)

;AB  magnitudes -> SDSS unprimed
;from www.ucolick.org/~cnaw/sun.html
;following Fukugita et al. 1995, PASP, 105, 945
s_ugrizJHK = fltarr(5)
s_ugrizJHK(0) = 6.75; //u SDSS
s_ugrizJHK(1) = 5.33; //g SDSS
s_ugrizJHK(2) = 4.67; //r SDSS
s_ugrizJHK(3) = 4.48; //i SDSS
s_ugrizJHK(4) = 4.42; //z SDSS

solar_mags = fltarr(14)
solar_mags(0) = 4.74 ;bolometric from Allen's Astrophysical Quantities, p. 341
;solar_mags(1) = s_UBVRIJHK(0)    ; brant had these incorrectly flipped
;solar_mags(2) = s_UBVRIJHK(1)    ; in his original version
;solar_mags(3) = s_UBVRIJHK(2)
;solar_mags(4) = s_UBVRIJHK(3)
;solar_mags(5) = s_UBVRIJHK(4)
;solar_mags(6) = s_UBVRIJHK(5)
;solar_mags(7) = s_UBVRIJHK(6)
;solar_mags(8) = s_UBVRIJHK(7)
;solar_mags(9) = s_ugrizJHK(0)
;solar_mags(10) = s_ugrizJHK(1)
;solar_mags(11) = s_ugrizJHK(2)
;solar_mags(12) = s_ugrizJHK(3)
;solar_mags(13) = s_ugrizJHK(4)
solar_mags(1) = s_ugrizJHK(0)
solar_mags(2) = s_ugrizJHK(1)
solar_mags(3) = s_ugrizJHK(2)
solar_mags(4) = s_ugrizJHK(3)
solar_mags(5) = s_ugrizJHK(4)
solar_mags(6) = s_UBVRIJHK(0)
solar_mags(7) = s_UBVRIJHK(1)
solar_mags(8) = s_UBVRIJHK(2)
solar_mags(9) = s_UBVRIJHK(3)
solar_mags(10) = s_UBVRIJHK(4)
solar_mags(11) = s_UBVRIJHK(5)
solar_mags(12) = s_UBVRIJHK(6)
solar_mags(13) = s_UBVRIJHK(7)



S = CALL_EXTERNAL('/home03/brant/code/idl/colors/colors', $
	'main', $
	long(N), $
	age, $           ;age lin gyr
	metallicity, $   ;metallicity in solar units
	luminosities, $  ;luminosities in abs magnitudes
	long(model), $	 ;model==0 -> Salpeter, model==1 Chabrier
	/F_VALUE)

for i=0,13 do luminosities(i,*) = (10.0^(-0.4*(luminosities(i,*)-solar_mags(i))))*mass(*)

return,luminosities ;in Solar Luminosities in each band (no factors of h)

end






;===========================================================================================
;
;From brobertson@cfa.harvard.edu Tue Jun 13 18:49:43 2006
;Date: Mon, 5 Jun 2006 23:54:52 -0400
;From: Brant Robertson <brobertson@cfa.harvard.edu>
;To: Thomas J. Cox <tcox@cfa.harvard.edu>,
;     Philip Hopkins <phopkins@cfa.harvard.edu>,
;     Yuexing Li <yxli@cfa.harvard.edu>,
;     Elisabeth Krause <ekrause@cfa.harvard.edu>,
;     Sukanya Chakrabarti <schakrab@cfa.harvard.edu>
;Cc: Brant Robertson <brobertson@cfa.harvard.edu>,
;     Lars Hernquist <lars@cfa.harvard.edu>
;Subject: new colors code w/ Chabrier IMF option
;
;Hi Everyone,
;
;	By popular demand, I have revised the colors code.
;
;	Please see
;
;		/home/brant/idl/colors/colors.pro
;
;	This script is an IDL function; you supply arrays of mass, age, and  
;metallicity and it returns the bolometric luminosity (in solar  
;luminosities) + luminosities in UBVRIJHK (Bessel + KPNO, in solar  
;luminosities in that band, Vega system) and Sloan ugriz (in solar  
;luminosities in that band, AB system).  The stellar population  
;modeling is from Bruzual and Charlot 2003.  Please read the  
;colors.pro file before using it, as the calling sequence, etc. has  
;changed slightly.
;
;	You can now choose between Salpeter and Chabrier IMF's with a simple  
;flag, which should make it much more useful for everyone since most  
;serious work nowadays is with an IMF closer to Chabrier than  
;Salpeter.  See the script for details.  Also, I've updated and  
;documented the Vega and AB absolute magnitudes of the sun.
;
;	There's an IDL script that calls the colors.pro function at
;
;		/home/brant/idl/colors/testing/test_colors.pro
;
;	Unless someone finds an error or odd behavior, there likely will be  
;no need to further update the colors code until the next Bruzual and  
;Charlot models appear.  If there are problems or questions, please  
;let me know.
;
;	The update was a significant amount of work but hopefully it was  
;worthwhile.  This code is designed as a tool for communal use, but  
;please let me know if you plan on using it for a project since  
;advanced notice will help reduce duplicate work or conflicts  
;(especially once I'm in Chicago).  In that vein, this code was  
;developed in part for use with the Cloudy-based PAH/near-IR  
;calculation I'm doing and should not be used for similar calculations.
;If you use the code in a paper, you should cite Bruzual and Charlot  
;2003, mention that you're using the Padova 1994 model version, and  
;state which IMF and magnitude system you are using.
;
;If you have questions about any of the above, please let me know.
;
;Thanks,
;Brant
;===========================================================================================
;
;
;
;
;Date: Fri, 16 Jun 2006 15:19:41 -0400
;From: Brant Robertson <brobertson@cfa.harvard.edu>
;To: Elisabeth Krause <ekrause@cfa.harvard.edu>,
;Philip Hopkins <phopkins@cfa.harvard.edu>,
;Thomas J. Cox <tcox@cfa.harvard.edu>
;Cc: Brant Robertson <brobertson@cfa.harvard.edu>
;Subject: new colors routine correction : band order
;
;Hi,
;
;The ordering of the bands in the luminosity array was incorrect
;relative to the solar magnitudes in the new version of the colors
;routine.
;
;The proper order is now:
;Bolometric u g r i z U B V R I J H K
;
;If you used the new colors code for a calculation before today, 
;you may have been using the wrong solar magnitude.  This should be
;fixed now. 
;
;Sorry for any confusion.
;-B
;
;
;
;
;===========================================================================================
;
;  Email from Brant to Stijn, Wed. 16 Aug 2006
;
;  Again, BC03 provide the final answers to us.  They aren't clear as to
;what U filter they use, although most U filters differ in their Vega
;solar absolute magnitude by only ~0.02.  I might use the KPNO
;filters, which you can find at:
;
;http://www.noao.edu/kpno/mosaic/filters/filters.html
;
;You'll still have to download the BC03 SED library.
;
;
;
;
;
;
;
;

;===========================================================================================
