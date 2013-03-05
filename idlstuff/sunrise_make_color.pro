;
;
; From chayward@cfa.harvard.edu Sat Oct 27 15:03:15 2007
; Date: Fri, 26 Oct 2007 15:42:51 -0400
; From: Chris Hayward <chayward@cfa.harvard.edu>
; To: Thomas J. Cox <tcox@cfa.harvard.edu>
; Subject: [Fwd: image script]
; 
; 
; 
; -------- Original Message --------
; Subject: image script
; Date: Thu, 4 Oct 2007 14:56:39 -0700
; From: Patrik Jonsson <patrik@governator.ucsc.edu>
; To: chayward@cfa.harvard.edu
; 
;
;


pro make_color,file,hdu,alpha,Q,m,scale=scale,bands=bands, imscale = imscale, image_path=image_path,one_each=one_each,noerase=noerase

sz=size(hdu)  
if(sz[0] gt 0) then begin
    start=hdu[0]
    stop=hdu[1]
end else begin
    start=hdu
    stop=hdu
end

if not keyword_set (image_path) then begin
    if not keyword_set(noerase) then $
      erase
    xyouts,4,8,file
end

for v=start,stop do begin
image = mrdfits (file, v)

Sz = size (image)
xs = sz [1]
ys = sz [2]
if not keyword_set(scale) then scale=[1.,1.,1.]
if not keyword_set (bands) then bands = [3,2,1]
r = image (*,*, bands [0])*scale[0]
g = image (*,*, bands [1])*scale[1]
b = image (*,*, bands [2])*scale[2]

;r (where (r < 0)) = 0
;g (where (g < 0)) = 0
;b (where (b < 0)) = 0

i = (r+g+b)/3+1e-20
print, "i: ", max(i), min(i)
print, "r: ", max(r), min(r)
print, "g: ", max(g), min(g)
print, "b: ", max(b), min(b)
;rr = r*asinh (alpha*Q*(I-m))/(Q*i)
;gg = g*asinh (alpha*Q*(I-m))/(Q*i)
;bb = b*asinh (alpha*Q*(I-m))/(Q*i)
rr = r
gg = g
bb = b
print, "rr: ", max(rr), min(rr)
print, "gg: ", max(gg), min(gg)
print, "bb: ", max(bb), min(bb)

;rr(where(i leq 0))=0
;gg(where(i leq 0))=0
;bb(where(i leq 0))=0

;--------

nclrs= 255

rr_scaled= nclrs * (rr - min(rr)) / (max(rr) -  min(rr))
gg_scaled= nclrs * (gg - min(gg)) / (max(gg) -  min(gg))
bb_scaled= nclrs * (bb - min(bb)) / (max(bb) -  min(bb))

pixels= sz[1]
img3d= fltarr(3,pixels,pixels)
img3d(0,*,*)= byte(rr_scaled)
img3d(1,*,*)= byte(gg_scaled)
img3d(2,*,*)= byte(bb_scaled)

;--------

img = [[[rr] ], [[gg] ], [[bb] ]]

;print, max (img)
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

if not keyword_set (imscale) then imscale = 1.0
imsize=3
margin=.5
lmargin=.5
bmargin=.5
imnum=v-start
;tv,b,v-start,true=3,xsize=1,/inches,/order
if not keyword_set (image_path) then begin
    if not keyword_set(one_each) then begin
        xpos=(lmargin+(imnum mod 2)*(imsize+margin))*imscale
        ypos=(bmargin+(imnum/2)*(imsize+margin))*imscale
    end else begin
        xpos=lmargin*imscale
        ypos=bmargin*imscale
    end
    tv,b,xpos,ypos,true=3,xsize=imsize,/inches,/order
end else begin
    ;file_name = image_path + '/'+ file + "-" + $
    file_name = file + "-" + $
                string (imnum, format = '(i3.3)')+".jpg"
    print, "Saving image " + file_name
    ;write_jpeg, file_name, b, true = 3,qual=99
    write_jpeg, file_name, img3d, true=1, quality=99
end

end
end


;-----------------------------------------------------------


pro make_color_all,base,view,alpha,Q,m,scale=scale,bands=bands, imscale = imscale, save_images= save_images,image_path=image_path,start_hdu=start_hdu,noerase=noerase

spawn,"ls "+base+"_???.fits",files

start_file=files [0]
if keyword_set(start_hdu) then $
  start_hdu = find_HDU(start_file, start_hdu) $
else $
  start_hdu = find_HDU(start_file, "CAMERA0-BROADBAND")

print, 'Generating color images starting with HDU ', start_HDU
if keyword_set (save_images) then begin
    if not keyword_set(image_path) then $
      image_path = "./color/"
    spawn, "mkdir "+image_path
end

for i =0,n_elements(files) - 1 do begin
    make_color,files[i],view+start_hdu,alpha,Q,m, scale = scale, $
               bands = bands, imscale =imscale, image_path = image_path, $
               noerase=noerase
end
end


