

;===========================================================================
;
; TJ adapted from code written by Patrik, Chris, and Desika
;
;==========================================================================


;broadbandfile='broadband_002.fits'
;broadbandfile='broadband_032.fits'
;broadbandfile='broadband_042.fits'
broadbandfile='broadband_092.fits'

; select which camera to look at
camera= 0   ; this has 0-4


; which bands do we want for the image?
; see the filters.txt file
; below, bands 3, 2, and 1 make up the image
;bands=[10,9,8]


; image processing (a la Lupton)
alpha= 0.2  ; 1.5
Q= 9.0      ; 20.0
m= 0.0      ; 0.01 also works


start_hdu = find_HDU(broadbandfile, "CAMERA0-BROADBAND")
camera= camera + start_hdu

image = mrdfits (broadbandfile, camera)

Sz = size (image)
xs = sz [1]
ys = sz [2]
if not keyword_set(scale) then scale=[1.,1.,1.]    ; scale each image
if not keyword_set (bands) then bands = [3,2,1]    ; which bands to use for image
r = image (*,*, bands [0])*scale[0]
g = image (*,*, bands [1])*scale[1]
b = image (*,*, bands [2])*scale[2]

;r (where (r < 0)) = 0
;g (where (g < 0)) = 0
;b (where (b < 0)) = 0

i = (r+g+b)/3+1e-20
print,max(i),max(r),max(g),max(b)
Rr = r*asinh(alpha*Q*(I-m))/(Q*i)
gg = g*asinh(alpha*Q*(I-m))/(Q*i)
bb = b*asinh(alpha*Q*(I-m))/(Q*i)
print,max(rr),max(gg),max(bb)

;rr(where(i leq 0))=0
;gg(where(i leq 0))=0
;bb(where(i leq 0))=0

; makes n x n x 3 array where each slice is an image
img = [[[rr] ], [[gg] ], [[bb] ]]

print, "Max(img)= ", max (img)
;stop
;if max(img) gt 1 then begin
;    for x = 0, xs- 1 do begin
;        for y = 0, ys- 1 do begin
;            mx =max (img (x, y,*))  
;            if (mx gt 1) then img (x, y,*) = img (x, y,*)/mx
;        end
;    end
;end

b=bytscl(img,min=0,max=1)
file_name = broadbandfile + "-" + string (camera, format = '(i3.3)')+".jpg"
;file_name = broadbandfile + "-" + string (camera, format = '(i3.3)')+"-bands2.jpg"
print, "Saving image " + file_name
write_jpeg, file_name, b, true = 3,qual=99

end

