pro sfr_multi, junk, $
		filename=filename, $
		cumulative=cumulative, $
		h=h


if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_multi, junk, filename=filename, /h"
   print, "  "
   print, "  "
   return
endif


initialize_plotinfo, 1


filename='minorsfr.eps'

setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;   Load the runs to display
;-------------------------------------------
; cooling
;fruns= ["Z2m-u1",   $
;	"Z2m-u1i",   $
;	"Z2m-u1k"]

; sf law
;fruns= ["Z2m-u34",   $
;       "Z2m-u12",   $
;       "Z2m-u35"]

; resoltuion
;fruns= ["Z2m4x-u1",   $
;       "Z2m2x-u1",   $
;       "Z2m-u1",  $
;	"Z2mh-u2",  $
;	"Z2mq-u2",  $
;	"Z2me-u2"]

; feedback n
;fruns= ["Z2m-u7",   $
;       "Z2m-u1",   $
;       "Z2m-u12",   $
;       "Z2m-u6"]

; alternate n's
;fruns= ["Z2m-u5",   $
;       "Z2m-u10",   $
;       "Z2m-u8",  $
;       "Z2m-u27a",  $
;       "Z2m-u28a",  $
;       "Z2m-u29a"]


; MW orientations - all of em
;fruns= ["Z2m-u1","Z9m-u1","Z10m-u1","Z15m-u1","Z11m-u1","Z12m-u1","Z13m-u1","Z14m-u1","Z17m-u1","Z16m-u1"]
;       tmerg= [2.7,2.7,2.7,2.8,2.8,3.4,3.8,4.6,7.5,10.0]
;       lbls= ['Z2','Z9','Z10','Z15','Z11','Z12','Z13','Z14','Z17','Z16']
;       msg= '(b)'




; mh comparison
;fruns= ["mhmaj-u1","mhmaj-u2","mhmaj-u3","mhmaj-u4","mhmaj-u5","mhmaj-u6","mhmaj-u7"]
;fruns= ["mhmaj-u4",   $
;       "mhmaj-u5"]
;fruns= ["mhmaj-u6","mhmaj-u7"]


;fruns= ["mhmaj-u6","mhmaj-u7","mhmaj-u5"]
;lbls= ["isothermal","energy: geometric","entropy"]



; springel comparison
;fruns= ["A1m-u1","A1m-u2","A1m-u3"]
;lbls= ["entropy", "energy: arithmetic (S00)", "energy: geometric"]


;fruns= ["A1m-u1","A1m-u4","A1m-u7"]
;lbls= ["c!I*!N=0.004 (S00)","c!I*!N=0.012", "c!I*!N=0.03"]


; gas distribution comparison
;fruns= ["A1m-u4",   $
;       "A1m-u5"]

;fruns= ["B1m-u4",   $
;       "B1m-u5"]


; kennicutt comparison
;fruns= ["A1m-u1",   $
;       "A1m-u4"]



; Sc kennicutt comparison
;fruns= ["Sc11i-u3",   $
;	"Sc11i-u2",   $
;	"Sc11i-u1",   $
;	"Sc11i-u4"]
;fruns= ["Sc11i-u11",   $
;       "Sc11i-u12",   $
;       "Sc11i-u13",   $
;       "Sc11i-u14",   $
;	"Sc11i-u15"]
;fruns= ["Sc11i-u5",   $
;       "Sc11i-u6",   $
;       "Sc11i-u7"]
;fruns= ["Sc11i-u8",   $
;       "Sc11i-u9",   $
;       "Sc11i-u10"]
;fruns= ["Sc11i-u16",   $
;       "Sc11i-u17",   $
;       "Sc11i-u18",   $
;       "Sc11i-u19"]
;fruns= ["Sc11i-u17","Sc11i-u8", "Sc11i-u32","Sc11i-u9","Sc11i-u19","Sc11i-u18","Sc11i-u34","Sc11i-u36"]
;fruns= ["Sc11i-u31","Sc11i-u10","Sc11i-u33","Sc11i-u35"]
;fruns= ["Sc11i-u16","Sc11i-u5", "Sc11i-u28","Sc11i-u6","Sc11i-u29","Sc11i-u30"]
;fruns= ["Sc11i-u25","Sc11i-u7", "Sc11i-u26","Sc11i-u27"]
;fruns= ["Sc11i-u11","Sc11i-u12","Sc11i-u1", "Sc11i-u2","Sc11i-u23","Sc11i-u24","Sc11i-u3"]
;fruns= ["Sc11i-u13","Sc11i-u4", "Sc11i-u12","Sc11i-u20","Sc11i-u14","Sc11i-u21","Sc11i-u22","Sc11i-u15"]

;fruns= ["Sbc11i4-u20","Sbc11i4-u14", "Sbc11i4-u22","Sbc11i4-u12"]
;fruns= ["Sbc11i4-u12","Sbc11i4-u24","Sbc11i4-u1"]
;fruns= ["Sbc11i4-u28","Sbc11i4-u29", "Sbc11i4-u30","Sbc11i4-u5"]

;fruns= ["Sbc11i4-u6","Sbc11i4-u16"]
;fruns= ["Sbc11i4-u8","Sbc11i4-u17"]
;fruns= ["Sbc11i4-u4","Sbc11i4-u38"]
;fruns= ["Sbc11i4-u4","Sbc11i4-u43","Sbc11i4-u8","Sbc11i4-u42","Sbc11i4-u40","Sbc11i4-u41"]
;fruns= ["Sbc11i4-u40","Sbc11i4-u41"] & lbls=["n1low","n1high"]
;fruns= ["Sbc11i4-u8","Sbc11i4-u42"] & lbls=["n0low","n0high"]
;fruns= ["Sbc11i4-u4","Sbc11i4-u43"] & lbls=["n2low","n2high"]



;fruns= ["Sc11i-u4", "Sbc11i4-u4"]


;fruns= ["Sbc11i-u3", "Sbc11i-u2", "Sbc11i-u1", "Sbc11i-u4"]
;fruns= ["Sbc11i-u5", "Sbc11i-u6", "Sbc11i-u7"]
;fruns= ["Sbc11i-u8", "Sbc11i-u9", "Sbc11i-u10"]

;fruns= ["Sbc201a-u4", "Sc201-u4"]
;fruns= ["Sbc201a-u8","Sbc201a-u32", "Sbc201a-u45", "Sbc201a-u50"]

;fruns= ["fbi-u1", "fbp1-u1", "fbp2-u1", "fbp3-u1", "fbp4-u1", "fbps21-u1"]



; -----------------
; Sc Major Mergers
; -----------------
;fruns= ["Sc201-u4","Sc203-u4","Sc206-u4","Sc212-u4"]
;       tmerg= [1.7,1.8,2.1,1.5]
;	lbls= ['Sc201','Sc203','Sc206','Sc212']




; Sbc comparison
;fruns= ["Sbc201a-u4","Sbc203-u4","Sbc206-u4","Sbc212-u4"]
;       tmerg= [1.7,1.8,2.1,1.5]
;       lbls= ['Sbc201','Sbc203','Sbc206','Sbc212']
;	msg= ' '





; ------------------
; Sbc major mergers
; ------------------

; c_star changes
;fruns= ["Sbc201a-u4", "Sbc201a-u14", "Sbc201a-u44"]  ; some difference
;fruns= ["Sbc201a-u43", "Sbc201a-u47"]    ; not much difference
;fruns= ["Sbc201a-u40", "Sbc201a-u48"]

; methods paper
; -------------
;fruns= ["Sbc11i4-u40","Sbc11i4-u46"] & lbls=["n1low","n1high"]
;  fruns= ["Sbc11i4-u53","Sbc11i4-u52"] & lbls=["n1low","n1high"]
;fruns= ["Sbc11i4-u8","Sbc11i4-u45"] & lbls=["n0low","n0high"]
;  fruns= ["Sbc11i4-u50","Sbc11i4-u51"] & lbls=["n0low","n0high"]
;  fruns= ["Sbc11i4-u67","Sbc11i4-u50","Sbc11i4-u51"] & lbls=["n0toolow","n0low","n0high"]
;fruns= ["Sbc11i4-u4","Sbc11i4-u43"] & lbls=["n2low","n2high"]
; -------------
;fruns= ["Sbc201a-u4","Sbc201a-u43"] & lbls=["n2low","n2high"]   ; paper version
;fruns= ["Sbc201a-u40","Sbc201a-u46"] & lbls=["n1low","n1high"]
;  fruns= ["Sbc201a-u53","Sbc201a-u52"] & lbls=["n1low","n1high"]   ; paper version
;fruns= ["Sbc201a-u8","Sbc201a-u45"] & lbls=["n0low","n0high"]
;  fruns= ["Sbc201a-u50","Sbc201a-u51"] & lbls=["n0low","n0high"]   ; paper version
; -------------
;fruns= ["Sbc11i4-u60","Sbc11i4-u61"] & lbls=["n2low","n2high"]   ; low density threshold
;fruns= ["Sbc11i4-u65","Sbc11i4-u64"] & lbls=["n1low","n1high"]
;fruns= ["Sbc11i4-u62","Sbc11i4-u63"] & lbls=["n0low","n0high"]
; -------------
; all of em
;fruns= ["Sbc201a-u4","Sbc201a-u43","Sbc201a-u53","Sbc201a-u52","Sbc201a-u50","Sbc201a-u51"]
;	tmerg= [1.8,1.8,1.8,1.8,1.8,1.8]



; cstar variations
;fruns= ["Sbc201a-u74","Sbc201a-u64","Sbc201a-u4","Sbc201a-u44","Sbc201a-u14","Sbc201a-u54"]
;	lbls= ["c!I*!N=0.3","c!I*!N=0.06","n2low","c!I*!N=0.015","c!I*!N=0.01","c!I*!N=0.004"]
;fruns= ["Sbc201a-u74","Sbc201a-u64","Sbc201a-u4","Sbc201a-u44","Sbc201a-u54"]
;       lbls= ["c!I*!N=0.3","c!I*!N=0.06","c!I*!N=0.03 (n2low)","c!I*!N=0.015","c!I*!N=0.004"]
;        lbls= ["c!I* !N=0.3","  !I !N 0.06","  !I !N 0.03 (n2low)","  !I !N 0.015","  !I !N 0.004"]




	tmerg= [1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8]
	msg= ' '



; orientations - all of em
;fruns= ["Sbc201a-u4", "Sbc202-u4", "Sbc203-u4", "Sbc206-u4", "Sbc207-u4", "Sbc217-u4"]
;	tmerg= [1.8,1.8,1.8,1.8,1.8,1.8]
;	lbls= ['Sbc201','Sbc202','Sbc203','Sbc206','Sbc207','Sbc217']
;	msg= '(a)'

; p-p, p-r, r-r
;fruns= ["Sbc201a-u4", "Sbc202-u4", "Sbc203-u4"]  ; R_peri= 11.0
;fruns= ["Sbc204-u4", "Sbc208-u4", "Sbc209-u4"]   ; R_peri= 5.0
;fruns= ["Sbc205-u4", "Sbc211-u4"]               ; R_peri= 44.0

; orbits
;fruns= ["Sbc201a-u4", "Sbc204-u4", "Sbc205-u4", "Sbc215-u4", "Sbc202-u4", "Sbc208-u4", "Sbc203-u4", "Sbc209-u4", "Sbc211-u4"]
;fruns= ["Sbc201a-u4", "Sbc204-u4", "Sbc205-u4", "Sbc208-u4", "Sbc209-u4", "Sbc211-u4", "Sbc215-u4"]
;	tmerg= [1.8,1.4,3.9,1.5,1.6,4.2,12.0]
;	lbls= ['Sbc201','Sbc204','Sbc205','Sbc208','Sbc209','Sbc211','Sbc215']
;fruns= ["Sbc201a-u4", "Sbc204-u4", "Sbc205-u4", "Sbc215-u4"]    ; prograde-prograde
;fruns= ["Sbc202-u4", "Sbc208-u4"]    ; prograde-retrograde
;fruns= ["Sbc203-u4", "Sbc209-u4", "Sbc211-u4"]    ; retrograde-retrograde

; slow orbits
;fruns= ["Sbc201a-u4","Sbc212-u4","Sbc218-u4","Sbc213-u4","Sbc214-u4","Sbc216-u4","Sbc219-u4","Sbc210-u4"]
;	tmerg= [1.8, 1.6, 1.7, 2.0, 2.6, 6.0, 1.8, 1.7]
;	lbls= ['Sbc201','Sbc212','Sbc218','Sbc213','Sbc214','Sbc216','Sbc219','Sbc210']
;fruns= ["Sbc212-u4", "Sbc218-u4"]   ; same slow orbit, one p-p the other p-r
;fruns= ["Sbc201a-u4", "Sbc212-u4"] ; slow-fast comparison (R_peri= 11.0)
;fruns= ["Sbc205-u4", "Sbc214-u4"] ; slow-fast comparison (R_peri= 44.0)
;fruns= ["Sbc215-u4", "Sbc216-u4"] ; slow-fast comparison (R_peri= 100.0)

; really long mergers
;fruns= ["Sbc201a-u4","Sbc215-u4","Sbc216-u4","Sbc211-u4","Sbc205-u4"]
;       tmerg= [1.8, 12.0, 6.0,4.2,3.9]
;       lbls= ['Sbc201','Sbc215','Sbc216','Sbc211','Sbc205']
;	msg= '(b)'
;fruns= ["Sbc201a-u4","Sbc214-u4","Sbc213-u4","Sbc212-u4","Sbc218-u4"]
;       tmerg= [1.8, 2.6, 2.0,1.6,1.7]
;       lbls= ['Sbc201','Sbc214','Sbc213','Sbc212','Sbc218']
;	msg= '(c)'
;fruns= ["Sbc201a-u4","Sbc204-u4","Sbc208-u4","Sbc209-u4"]
;       tmerg= [1.8, 1.4, 1.5,1.6]
;       lbls= ['Sbc201','Sbc204','Sbc208','Sbc209']
;	msg= '(d)'
;fruns= ["Sbc201a-u4","Sbc210-u4","Sbc219-u4"]
;       tmerg= [1.8,1.7,1.8]
;       lbls= ['Sbc201','Sbc210','Sbc219']
;	msg= '(e)'



;fruns= ["Sbc201a-u4", "Sbc202-u4", "Sbc203-u4", "Sbc206-u4", "Sbc207-u4", "Sbc217-u4","Sbc212-u4","Sbc218-u4","Sbc213-u4","Sbc214-u4","Sbc216-u4","Sbc219-u4","Sbc210-u4"]
;lbls= [" "," "," "," "," "," "," "," "," "," "," "," "," "," "]
msg= ' '



; resolution
;fruns= ["Sbc201a-u4", "Sbc201a-n4", "Sbc201a2x-n4", "Sbc201a4x-n4", "Sbc201a10x-n4"]
;fruns= ["Sbc201a-u4", "Sbc201a2x-n4", "Sbc201a4x-n4", "Sbc201a10x-n4", "Sbc201a10x-u4"]
;lbls= ['1x (n2low)', '2x', '4x', '10x', '10x, h=0.047']
;fruns= ["Sbc201a-u4", "Sbc201a10x-n4", "Sbc201a10x-u4"]
;lbls= ['n2low', '10x','10x-h']
;fruns= ["Sbc201a10x-n4", "Sbc201a10x-u4"]
;lbls= ['h=0.100','h=0.047']

;msg= ' '
;tmerg=[1.8,1.8,1.8,1.8,1.8,1.8]


; sph version
;fruns= ["Sbc201a-u4", "Sbc201a-u4d", "Sbc201a-u4c"] & lbls= ["n2low","arithmetic","geometric"]


; zero angular momentum
;fruns= ["Sbc203-u4", "Sbc210-u4","Sbc219-u4"]


; plus/minus 1 sigma
;fruns= ["Sbc201a-u4","Sbcm201-u4","Sbcp201-u4"]
;	tmerg= [1.8,2.1,1.6]
;	lbls= ['Sbc201','Sbc-1','Sbc+1']
;        msg= '(f)'


; different cooling
;fruns= ["Sbc201a-u4", "Sbc201a-u4f", "Sbc201a-u4g", "Sbc201a-u4h", "Sbc201a-u4i"]
;	tmerg=[1.8,1.8,1.8,1.8,1.8]
;	lbls= ['KWH','mzero.cie','m-30.cie','m-05.cir','m-00.cie']
;	msg= ' '
;fruns= ["Sbc201a-u4", "Sbc201a-u4f", "Sbc201a-u4g", "Sbc201a-u4k", "Sbc201a-u4j", "Sbc201a-u4i"]
;       tmerg=[1.8,1.8,1.8,1.8,1.8,1.8]
;       lbls= ['KWH','mzero.cie','m-30.cie','m-20.cie','m-10.cie','m-00.cie']
;       msg= ' '
;fruns= ["Sbc201a-u4", "Sbc201a-u4f", "Sbc201a-u4g", "Sbc201a-u4k", "Sbc201a-u4l", "Sbc201a-u4j", "Sbc201a-u4i"]
;       tmerg=[1.8,1.8,1.8,1.8,1.8,1.8,1.8]
;       lbls= ['KWH','mzero.cie','m-30.cie','m-20.cie','m-15.cie','m-10.cie','m-00.cie']
;       msg= ' '
;fruns= ["Sbc201a-u4f", "Sbc201a-u4i"]
;       tmerg=[1.8,1.8]
;       lbls= ['mzero.cie','m-00.cie']
;       msg= ' '


; different bulge
;fruns= ["Sbc201a-u4", "Sbc201a-u4hbl"]
;       tmerg=[1.8,1.8]
;       lbls= ['Exponential','Hernquist']
;       msg= ' '

; radial cut off
;fruns= ["Sbc201a-u4", "Sbc201ac60-u4","Sbc201ac80-u4"]
;       tmerg=[1.8,1.8,1.8]
;       lbls= ['Fiducial','80 cutoff','60 curoff']
;       msg= ' '


; matt covington's energy detail
;fruns= ["Sbc201a-u4","Sbc201a-u4p"]
;tmerg=[1.8,1.8]
;msg=' '

; -------------------------------------

;fruns= ["V5m-u1","V7m-u1", "V4m-u1", "V2m-u1","V1m-u1", "V3m-u1", "V6m-u4"]
;	lbls= ['V5','V7','V4','V2','V1','V3','V6']
;        tmerg= [1.2,2.1,2.0,1.4,2.0,2.3,1.7]
;fruns= ["Y7mf-u1","Y1mf-u1","Y2mf-u1","Y3mf-u1","Y4mf-u1","Y5mf-u1","Y6mf-u1"]  
;	lbls= ['Y7','Y1','Y2','Y3','Y4','Y5','Y6']
;        tmerg= [2.1,2.3,2.6,3.0,3.6,4.3,5.5]
;fruns= ["D7mf-u1","D1mf-u1","D2mf-u1","D3mf-u1","D4mf-u1","D5mf-u1","D6mf-u1"]
;	lbls= ['D7','D1','D2','D3','D4','D5','D6']
;        tmerg= [2.0,2.0,2.2,2.7,3.2,4.0,5.2]
;fruns= ["V1m-u1", "V3m-u1"]

;fruns= ["V6m-u1","V6m-u2","V6m-u3"]
;	lbls= fruns

;fruns= ["X1i-u1","Y1i-u1","D3i-u1"]



; ------------------------------------
;  G Models

;fruns= ["G2i-u1", "G2in-u1", "G2im-u1"]

;fruns= ["G0i-u1", "G1i-u1", "G2im-u1", "G3il-u1"] & lbls= ['G0','G1','G2','G3']
fruns= ["G0i-u1a", "G1i-u1a", "G2im-u1a", "G3il-u1a"] & lbls= ['G0','G1','G2','G3']
oldstars= [	0.098 + 0.002, $
		0.47 + 0.03, $
		1.35 + 0.15, $
		4.11 + 0.89]

;fruns= ["G3G0-u1", "G3G0a-u1", "G3G0b-u1", "G3G0c-u3"]
;fruns= ["G3G3a-u1", "G3G3b-u1"]

; prograde
;fruns= ["G3G3b-u1", "G3G2-u3", "G3G1-u3", "G3G0e-u3"] & lbls= ['G3G3 (2.5)','G3G2 (2.9)','G3G1 (4.2)','G3G0 (6.0)']
;fruns= ["G2G2-u1", "G2G1-u3", "G2G0-u3"] & lbls= ['G2G2 (1.3)','G2G1 (1.5)','G2G0 (2.7)']
;fruns= ["G1G1a-u1", "G1G0-u3"] & lbls= ['G1G1 (1.3)','G1G0 (2.1)']
;fruns=  ["G0G0a-u1"] & lbls= ['G0G0 (1.4)']

; smaller majors with lower smoothing length
;fruns=  ["G0G0a-u1","G0G0a-u3"] & lbls= ['G0G0 (1.4)','G0G0 (h=0.035)']
;fruns=  ["G1G1a-u1","G1G1a-u3"] & lbls= ['G1G1 (1.3)','G1G1 (h=0.055)']

; retrograte
;fruns= ["G3G3r-u1","G3G2r-u3","G3G1r-u3"] & lbls= ['G3G3r','G3G2r','G3G1r']
;fruns= ["G2G2r-u1","G2G1r-u3","G2G0r-u3"] & lbls= ['G2G2r','G2G1r','G2G0r']
;fruns= ["G1G1r-u1","G1G0r-u3"] & lbls= ['G1G1r','G1G0r']
;fruns=  ["G0G0r-u1b"] & lbls= ['G0G0r']


;fruns= ["G3il-u1","G3bli-u1","G3bliv2-u1","G3bliv3-u1","G3bliv4-u1","G3bliv5-u1","G3bliv6-u1"]
;fruns= ["G3il-u1","G3bliv5-u1"]
;fruns= ["G3G3b-u1","G3blG3bl-u1"]
;fruns= ["G3G3b-u1","G3blv5G3blv5-u1"]
;fruns= ["G3blv5G3blv5-u1","G3blv5G2-u3","G3blv5G1-u3","G3blv5G0-u3"] & lbls= ['G3blG3bl','G3blG2','G3blG1','G3blG0'] & tmerg= [2.4, 2.9, 4.5, 6.0]


; just majors
;fruns= ["G3G3a-u1", "G3G3b-u1","G3G3r-u1"] & tmerg=[3.4,2.5,2.6] & lbls=['G3a','G3b','G3r']
;fruns= ["G2G2-u1","G2G2r-u1"] & tmerg=[1.3,1.3] & lbls=['G2','G2r']
;fruns= ["G1G1a-u1","G1G1r-u1"] & tmerg=[1.3,1.3] & lbls=['G1','G1r']
;fruns= ["G0G0a-u1","G0G0r-u1b"] & tmerg=[1.4,1.4] & lbls=['G0','G0r']
;fruns= ["G3blv5G3blv5-u1"] & tmerg= [2.4,2.4] & lblbs=['G3bl']

; ------------------------------------------------------

; vary gas fraction
;fruns= ["G3G3b-u1","G3gf1G3gf1b-u1","G3gf2G3gf2b-u1"]
;lbls= ['20%','50%','70%']
;fruns= ["G3gf4G3gf4","G3G3b-u1","G3gf1G3gf1b-u1","G3gf2G3gf2b-u1","G3gf3G3gf3"]
;lbls= ['11%','20%','42%','58%','75%']
msg= ' '

; vary bulge fraction
;fruns= ["G3G3b-u1","G3BT1G3BT1-u1","G3BT2G3BT2-u1","G3BT3G3BT3-u1"]
;fruns= ["G3G3b-u1","G3blv5G3blv5-u1","G3BT1G3BT1-u1","G3BT2G3BT2-u1","G3BT3G3BT3-u1"]
;fruns= ["G3BT1G3BT1-u1","G3BT1G2-u1","G3BT1G1-u1","G3BT1G0-u1"]
;fruns= ["G3BT2G3BT2-u1","G3BT2G2-u1","G3BT2G1-u1","G3BT2G0-u1"]
;fruns= ["G3BT3G3BT3-u1","G3BT3G2-u1","G3BT3G1-u1","G3BT3G0-u1"]
msg= ' '

;fruns= ["G3il-u1a","G3bliv5-u1","G3BT1i-u1","G3BT2i-u1","G3BT3i-u1"]
;msg= ' '

; minors with different n (feedback n)
;fruns= ["G3G3b-u1","G3G1-u3","G3G1-u4"]
;fruns= ["G3G3b-u1","G3G2-u3","G3G2-u4"]
;fruns= ["G3G3b-u1","G3G3b-u2","G3G0e-u3","G3G0e-u4"]

;fruns= ["G3G3b-u2", "G3G2-u4", "G3G1-u4", "G3G0e-u4"]
;fruns= ["G2G2-u2", "G2G1-u4", "G2G0-u4"]
;fruns= ["G1G1a-u2", "G1G0-u4"]
;fruns=  ["G0G0a-u2"]
;msg= ' '

;fruns= ["G2G2-u1", "G2G0-u3", "G2G0-u4"]

;fruns= ["G3gf1G3gf1b-u1","G3gf1G2-u1","G3gf1G1-u1","G3gf1G0-u1"]
;fruns= ["G3gf2G3gf2b-u1","G3gf2G2-u1","G3gf2G1-u1","G3gf2G0-u1"]
;msg= ' '

;fruns= ["G3G2a-u1","G3G2-u3","G3G2c-u1","G3G2b-u1","G3G2r-u3","G3G2d-u1"]
;fruns= ["G3G2e-u1","G3G2-u3","G3G2f-u1","G3G2g-u1","G3G2h-u1"]

;fruns= ["G3G1a-u1","G3G1-u3","G3G1c-u1","G3G1b-u1","G3G1r-u3","G3G1d-u1"]


; Mako & UpsAnd comparison
;fruns= ["G3G3b-u1", "G3blv5G3blv5-u1","G3bG3b"]
;fruns= ["G3BT1G3BT1-u1", "G3b1G3b1"]
;fruns= ["G3BT2G3BT2-u1", "G3b2G3b2"]
;fruns= ["G3BT3G3BT3-u1", "G3b3G3b3"]
msg= ' '

; Bulge-to-disk ratios
;    -- exp. bulge (upsand) ---
;fruns= ["G3blv5G3blv5-u1", "G3BT1G3BT1-u1","G3BT2G3BT2-u1","G3BT3G3BT3-u1"]
;    -- exp. bulge (mako) ---
;fruns= ["G3bG3b", "G3b1G3b1","G3b2G3b2","G3b3G3b3"]
;    -- hernquist bulge (mako) ---
;fruns= ["G3bG3b", "G3b6G3b6","G3b5G3b5","G3b4G3b4"]


; Bulge-to-disk ratios  & gas fraction
;fruns= ["G3b5_gf3_mm","G3b5_gf2_mm","G3b5_gf1_mm","G3b5G3b5", "G3b5_gf4_mm"]
;fruns= ["G3b5_gf3_mm","G3b5_gf2_mm","G3b5G3b5", "G3b5_gf4_mm"]






; ---------------------------------------------------

;fruns= ["I1i-u1", "I1i-u2", "I1i-u3", "I1i-u4"]
;fruns= ["I1i-u1", "I1i-u4", "I1i-u5"]
;fruns= ["I1i-u5", "I1i-u10", "I1i-u11","I1i-u12"]
;fruns= ["I1i-u5", "I1i-u13", "I1i-u14","I1i-u15"]
;fruns= ["I1gf1i-u1", "I1gf1i-u2", "I1gf1i-u3", "I1gf1i-u4"]
;fruns= ["I1gf1i-u1", "I1gf1i-u4", "I1gf1i-u5"]
;fruns= ["I1gf2i-u1", "I1gf2i-u2", "I1gf2i-u3", "I1gf2i-u4"]
;fruns= ["I1gf2i-u1", "I1gf2i-u4", "I1gf2i-u5"]


; ---------------------------------------------------



zcolor= 50
deltacolor= (250-zcolor)/n_elements(fruns)
;deltacolor= 0
lstyle= 0
zcolor_orig= zcolor
;thickstyle= 2.0
thickstyle= 5.0

;--------------------------------------
;  Print the Shit
;--------------------------------------


yaxistitle="!5SFR (M!D!9n!5!N yr!E-1!N)"
if keyword_set(cumulative) then yaxistitle="Mass (10!e10!n h!E-1!NM!D!9n!3!N)"
yaxistitle="!5SFR / Stellar Mass (0.1 Gyr!E-1!N) "
xaxistitle = "Time (Gyr/h)"
xaxistitle = "Time (Gyr)"

;xmax = 14.0
;xmax = 10.0
;xmax = 8.0
;xmax = 7.5
xmax = 6.0
;xmax = 5.0
;xmax = 4.25
;xmax = 4.0
;xmax = 3.0
;xmax = 2.8
;xmax = 2.0
;xmax = 1.5
;xmax = 1.3
;xmax = 1.0
xmin = 0

;ymax = 300
;ymax = 180
;ymax = 120
;ymax = 90
;ymax = 80
;ymax = 60.0
;ymax = 40.0
;ymax = 30.0
;ymax = 20.0
;ymax = 15
;ymax = 13
;ymax = 12
;ymax = 10.0
;ymax = 6.0
;ymax = 5.0
;ymax = 3.0
;ymax = 2.0
;ymax = 1.5
;ymax = 1.3
;ymax = 1.2
ymax = 1.0
;ymax = 0.75
;ymax = 0.6
;ymin = 0 
;ymin = 0.01
ymin = 0.001
;ymin = 0.00005


; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]
;!p.font= 0

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	ytickformat='exp_label', $
	/nodata

;stop


        ifinal= n_elements(fruns)-1
        for i=0,ifinal do begin
        ; get sfr data
	;-------------------------------------------
        minor_open_sfr_file, fruns[i], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass

        time= transpose(sfrtime)

        sfr= transpose(sfrsfr)
        sm= transpose(sfrmfs)
        gm= transpose(sfrgasmass)


	    ; physical units
	    if keyword_set(h) then sfrtime = sfrtime / h
	    ;if i eq 1 then sfrtime = sfrtime/(0.7)    ; did this for comparing h and noh


            if keyword_set(cumulative) then begin
                oplot, sfrtime, sfrmfs, psym=-3, color= zcolor, linestyle= lstyle, thick= 4.0
            endif else begin
                oplot, sfrtime, sfrsfr/(sfrmfs+oldstars[i]), psym=-3, color= zcolor, linestyle= lstyle, thick= thickstyle
            endelse
 
	    lstyle= lstyle+1
	    ;thickstyle= thickstyle+2.0
	    zcolor= zcolor+deltacolor
	endfor



; cooling
;xyouts, 0.7, 0.86, "Cooling", /normal, charthick=1, size=1.33, color=0
;oplot, [1.65,2.2],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.7, 0.80, "Explicit", /normal, charthick=1, size=1.33, color=50
;xyouts, 0.7, 0.75, "mzero.cie", /normal, charthick=1, size=1.33, color=100
;xyouts, 0.7, 0.70, "m-00.cie", /normal, charthick=1, size=1.33, color=150

; sfr law
;;xyouts, 0.7, 0.86, "SF Law", /normal, charthick=1, size=1.33, color=0
;xyouts, 0.74, 0.86, '!7s!3!N!D!20E!3!N', /normal, charthick=1, size=1.33, color=0
;oplot, [1.65,2.2],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.7, 0.80, "constant", /normal, charthick=1, size=1.33, color=50
;xyouts, 0.7, 0.75, "dynamical", /normal, charthick=1, size=1.33, color=100
;xyouts, 0.7, 0.70, '1/!7q!3!N', /normal, charthick=1, size=1.33, color=150

; resolution
;xyouts, 0.7, 0.86, "Resolution", /normal, charthick=1, size=1.33, color=0
;oplot, [1.65,2.2],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.7, 0.80, '4', /normal, charthick=1, size=1.33, color=50
;xyouts, 0.7, 0.76, '2', /normal, charthick=1, size=1.33, color=75
;xyouts, 0.7, 0.72, '1', /normal, charthick=1, size=1.33, color=100
;xyouts, 0.7, 0.68, '1/2', /normal, charthick=1, size=1.33, color=125
;xyouts, 0.7, 0.64, '1/4', /normal, charthick=1, size=1.33, color=150
;xyouts, 0.7, 0.60, '1/8', /normal, charthick=1, size=1.33, color=175

; feedback - within on n
;xyouts, 0.65, 0.86, 'n=2 Feedback', /normal, charthick=1, size=1.33, color=0
;;xyouts, 0.65, 0.86, '!7b!3 for n=2', /normal, charthick=1, size=1.33, color=0
;oplot, [1.50,2.35],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.65, 0.80, 'high', /normal, charthick=1, size=1.33, color=50
;xyouts, 0.65, 0.75, 'medium high', /normal, charthick=1, size=1.33, color=87
;xyouts, 0.65, 0.70, 'medium low', /normal, charthick=1, size=1.33, color=124
;xyouts, 0.65, 0.65, 'low', /normal, charthick=1, size=1.33, color=161

; alternate feedback n's
;xyouts, 0.55, 0.86, 'Alternate Feedback n', /normal, charthick=1, size=1.33, color=0
;oplot, [1.20,2.35],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.6, 0.80, 'n=1 high fb', /normal, charthick=1, size=1.33, color=50
;xyouts, 0.6, 0.76, 'n=1 medium fb', /normal, charthick=1, size=1.33, color=75
;xyouts, 0.6, 0.72, 'n=1 low fb', /normal, charthick=1, size=1.33, color=100
;xyouts, 0.6, 0.68, 'n=0 high fb', /normal, charthick=1, size=1.33, color=125
;xyouts, 0.6, 0.64, 'n=0 medium fb', /normal, charthick=1, size=1.33, color=150
;xyouts, 0.6, 0.60, 'n=0 low fb', /normal, charthick=1, size=1.33, color=175

; mh comparison
;xyouts, 0.6, 0.85, "Isothermal Gas", /normal, charthick=1, size=1.33, color=50
;xyouts, 0.6, 0.80, "Full Gas Treatment", /normal, charthick=1, size=1.33, color=150


; springel comparison
;xyouts, 0.3, 0.80, "Conservative Entropy", /normal, charthick=1, size=1.33, color=50
;xyouts, 0.3, 0.75, "Arithmetic", /normal, charthick=1, size=1.33, color=100
;xyouts, 0.3, 0.70, "Geometric", /normal, charthick=1, size=1.33, color=150


; minor paper - iso figure
; ------------------------------
xyouts, 0.80, 0.90, 'G3', /normal, charthick=3, size=1.33, color=zcolor_orig+deltacolor+deltacolor+deltacolor
xyouts, 0.80, 0.72, 'G2', /normal, charthick=3, size=1.33, color=zcolor_orig+deltacolor+deltacolor
xyouts, 0.80, 0.52, 'G1', /normal, charthick=3, size=1.33, color=zcolor_orig+deltacolor
xyouts, 0.80, 0.32, 'G0', /normal, charthick=3, size=1.33, color=zcolor_orig
toppos= 0.32
botpos= 0.01
;xarrow= 0.60
;yarrow= 0.71
;arrow, xarrow, toppos, yarrow, botpos, /data, COLOR=zcolor_orig, THICK=2.0, hthick=3.0




xyouts, 0.25, 0.85, msg, /normal, charthick=3.0, size=1.5, color= 0



; std
;std= 0
std= 1
if std eq 1 then begin
    ;x0= 0.75
    x0= 0.23
    y0= 0.92
    zcolor= zcolor_orig
    ;thickstyle= 2.0
    thickstyle= 5.0
    for i=0, n_elements(fruns)-1 do begin
	y0= y0-0.04
	if n_elements(lbls) gt 0 then begin
                ;xyouts, x0, y0, lbls[i], /normal, charthick=thickstyle, size=1.5, color=zcolor
                ;xyouts, x0, y0, lbls[i], /normal, charthick=thickstyle, size=1.333, color=zcolor
                    zcolor= zcolor+deltacolor
		    ;thickstyle= thickstyle+2.0
	endif else begin
		;xyouts, x0, y0, fruns[i], /normal, charthick=3, size=1.5, color=zcolor
	            zcolor= zcolor+deltacolor
	endelse
    endfor
endif










;--------------------------------------
;--------------------------------------

device, /close


end


