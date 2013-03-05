;+
; NAME:
;       ERROR_BARS.pro
;
; PURPOSE:
;	This procedure plots error bars over plotted data.
;
; CATEGORY:
;       Graphics
;
; CALLING SEQUENCE:
;       ERROR_BARS, X, Ymin, Ymax
;
; INPUTS:
;	X:  A vector containing the X-axis values of the points, of type
;	    integer or floating point.
;	Ymin:  A vector containing the minimum error limits of the
;	       points on the Y-axis, of type integer or floating point.
;	Ymax:  A vector containing the maximum error limits of the
;	       points on the Y-axis, of type integer or floating point.
;
; KEYWORD PARAMETERS:
;	COLOR:  The color index value of the error bars.  The default
;	        is the IDL default plot color.
;	THICK:  The thickness of the line plotting the error bars.  The
;	        default is a thickness of 1.
;	WIDTH:  The length of the horizontal end lines of the error
;	        bars.  The default is 0 (no end lines).
;	ROTATE:  If set, the procedure switches the X- and Y-axes,
;	         and so plots horizontal error bars on a vertical
;	         plot.  The default is to plot vertical error bars on a
;	         horizontal plot.
;
; USES:
;	VAR_TYPE.pro
;
; PROCEDURE:
;	This procedure plots lines forming vertical error bars over
;	a pre-existing plot.
;
; EXAMPLE:
;	Plot three points with error bars of width 1.
;	  plot, [1.,2.,3.], [1.1,1.2,1.3], psym=4
;	  error_bars, [1.,2.,3.], [1.1,1.2,1.3]-0.5, [1.1,1.2,1.3]+0.5
;
; MODIFICATION HISTORY:
;	Written by:	Daithi A. Stone, (stoned@uvic.ca) 2001-02-14.
;-

;***********************************************************************

PRO ERROR_BARS, X, $
                Ymin, Ymax, $
                COLOR=color, $
                THICK=thick, $
                WIDTH=width, $
                ROTATE=rotateopt

;***********************************************************************
;Default Settings

;Color
if var_type(color) eq 0 then color = !p.color

;Line thickness
if not(keyword_set(thick)) then thick = 1

;Length of end lines of error bar
if var_type(width) eq 0 then width = 0

;Input vectors
xv = [x]
yminv = [ymin]
ymaxv = [ymax]

;Length of vectors
nx = n_elements(xv)

;Rotation option
rotateopt = keyword_set(rotateopt)

;***********************************************************************
;Plot the Error Bars

for i=0,nx-1 do begin
  ;Main error bar
  if rotateopt then begin
    oplot, [yminv[i],ymaxv[i]], [xv[i],xv[i]], thick=thick, color=color
  endif else begin
    oplot, [xv[i],xv[i]], [yminv[i],ymaxv[i]], thick=thick, color=color
  endelse
  ;Horizontal end lines of the error bar
  if width ne 0 then begin
    ;Lower end line of error bar
    if rotateopt then begin
      oplot, [yminv[i],yminv[i]], [xv[i]-width/2.,xv[i]+width/2.], $
             thick=thick, color=color
    endif else begin
      oplot, [xv[i]-width/2.,xv[i]+width/2.], [yminv[i],yminv[i]], $
             thick=thick, color=color
    endelse
    ;Upper end line of error bar
    if rotateopt then begin
      oplot, [ymaxv[i],ymaxv[i]], [xv[i]-width/2.,xv[i]+width/2.], $
             thick=thick, color=color
    endif else begin
      oplot, [xv[i]-width/2.,xv[i]+width/2.], [ymaxv[i],ymaxv[i]], $
             thick=thick, color=color
    endelse
  endif
endfor

;***********************************************************************
;The End

return
END
