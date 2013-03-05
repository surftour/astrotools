pro lg_mollview, file_in, select_in, $
CHARSIZE=charsize, COLT = colt, COORD = coord, CROP = crop, $
FACTOR=factor, FLIP=flip, $
GAL_CUT=gal_cut, GIF = gif, GRATICULE = graticule, $
HELP = help, $
HIST_EQUAL = hist_equal, HXSIZE = hxsize, $
IGRATICULE=igraticule, $
LOG = log, $
MAX = max_set, MIN = min_set, $
NESTED = nested_online, NOBAR = nobar, NOLABELS = nolabels, NO_DIPOLE=no_dipole, NO_MONOPOLE=no_monopole, $
OFFSET = offset, ONLINE = online, OUTLINE=outline, $
PNG=png, POLARIZATION=polarization, PREVIEW = preview, PS = ps, PXSIZE = pxsize, $
QUADCUBE = quadcube, $
ROT = rot, $
SAVE = save, SUBTITLE = subtitle, $
TITLEPLOT = titleplot, $
UNITS = units, XPOS = xpos, YPOS = ypos

;+
; NAME:
;	MOLLVIEW, GNOMVIEW, CARTVIEW, ORTHVIEW
;
; PURPOSE:
; 	tools to view a Mollweide/gnomic/cartesian/orthographic projection of maps binned
;	in Healpix or COBE Quad-Cube pixelisation
;
; CALLING SEQUENCE:
; 	xxxxVIEW, File, [Select, ] $
;                       [CHARSIZE=, COLT=, COORD=, CROP=, $
;                       FACTOR=, FITS=, FLIP=, $
;                       GAL_CUT=, GIF=, GRATICULE=, $
;                       HALF_SKY =, HELP =, HIST_EQUAL=, HXSIZE=, $
;                       IGRATICULE=igraticule, $
; 	                LOG=, $
;                       MAX=, MIN=, $ 
; 	                NESTED=, NOBAR=, NOLABELS=, NOPOSITION =, $
;                       OFFSET =, ONLINE=, OUTLINE=, $
; 	                PNG=, POLARIZATION=, PREVIEW=,$
;                       PS=, PXSIZE=, PYSIZE=, $
;                       QUADCUBE= , $
;                       NO_DIPOLE=, NO_MONOPOLE=, $
;                       RESO_ARCMIN= , ROT=, $
;                       SAVE=, SUBTITLE=, $
;                       TITLEPLOT=, $
; 	                UNITS=, XPOS=, YPOS=]
;                        
;  all the arguments and parameters are identical for all the
;  routines, excepted stated otherwise.
;
;
; INPUTS:
; 	File = 
;          by default,           name of a FITS file containing 
;               the healpix map in an extension or in the image field
;          if Online is set :    name of a variable containing
;               the healpix map
;          if Save is set   :    name of an IDL saveset file containing
;               the healpix map stored under the variable  data
;
; OPTIONAL INPUTS:
;       Select =  if the file read is a multi column BIN table, Select indicates
;                 which column is to be plotted (the default is to plot the first one)
;               can be either a name : value given in TTYPEi of the FITS file
;                        NOT case sensitive and
;                        can be truncated, 
;                        (only letters, digits and underscore are valid)
;               or an integer        : number i of the column
;                            containing the data, starting with 1
;                   the default value is 1
;               (see the Examples below)
;
; OPTIONAL INPUT KEYWORDS:
;       CHARSIZE : overall multiplicative factor applied to the size of all
;               characters appearing on the plot
;                default = 1.0
;
; 	COLT : color table number, in [0,40]
;        	if unset the color table will be 33 (Blue-Red)
;
;       COORD : vector with 1 or 2 elements describing the coordinate system of the map 
;                either 'C' or 'Q' : Celestial2000 = eQuatorial,
;                       'E'        : Ecliptic,
;                       'G'        : Galactic 
;               if coord = ['x','y'] the map is rotated from system 'x' to system 'y'
;               if coord = ['y'] the map is rotated to coordinate system 'y' (with the
;               original system assumed to be Galactic unless indicated otherwise in the file)
;                  see also : Rot
;
;       CROP : if set the image file (gif, png) only contains the mollweide map and
;               not the title, color bar, ...
;               (see also : GIF, PNG)
;
;       FACTOR : multiplicative factor to be applied to the data (default = 1.0)
;               the data plotted is of the form FACTOR*(data + OFFSET)
;               see also : OFFSET, LOG
;
;       FITS : string containing the name of an output fits file with
;       the rectangular map in the primary image
;	      if set to 0 or not set : no .FITS done
;	      if set to 1            : output the plot in plot_gnomic.fits
;	      if set to a file name  : output the plot in that file 
;             ** gnomview only **
;
;       FLIP : if set the longitude increase to the right, whereas by
;               default (astro convention) it increases towards the left
;
;       GAL_CUT: (positive float) specifies the symmetric galactic cut in degree
;             outside of which the the monopole and/or dipole fitting is done
;             (see also: NO_DIPOLE, NO_MONOPOLE)
;             ** mollview only **
;
;	GIF : string containing the name of a .GIF output
;	      if set to 0 or not set : no .GIF done
;	      if set to 1            : output the plot in plot_mollweide.gif
;	      if set to a file name  : output the plot in that file 
;             (see also : CROP, PNG, PS and PREVIEW)
;
; 	GRATICULE : if set, puts a graticule with delta_long = delta_lat = default
;         if graticule is set to a scalar > gmin delta_long = delta_lat = graticule
;         if set to [x,y] with x,y > gmin then delta_long = x and delta_let = y
;         ** cartview : default =  5, gmin =  0 **
;         ** gnomview : default =  5, gmin =  0 **
;         ** mollview : default = 45, gmin = 10 **
;         ** orthview : default = 45, gmin = 10 **
;
;       HALF_SKY: if set, only shows only one half of the sky 
;          (centered on (0,0) or on the location parametrized by Rot) instead of the full sky
;             ** orthview only **
;        
;       HELP : if set, the routine header is printed (by doc_library)
;             and nothing else is done
;
;	HIST_EQUAL : if not set, uses linear color mapping and 
;                     		puts the level 0 in the middle
;                     		of the color scale (ie, green for Blue-Red)
;				unless MIN and MAX are not symmetric
;	       	      if set,     uses a histogram equalized color mapping
;			(useful for non gaussian data field)
;                     (see also : LOG)
;
; 	HXSIZE: horizontal dimension (in cm) of the Hardcopy plot : Only for postscript printout
;    		** mollview : default = 26 cm ~ 10 in **
;               ** mollview : default = 15 cm         **
;    		(useful for large color printer)
;               (see also : PXSIZE)
;
;       IGRATICULE: if set, puts a graticule in the input coordinates
;          if both graticule and igraticule are set, these ones will
;          be represented by dashes
;
; 	LOG: display the log of map (see also : HIST)
;         only applies to pixel with signal > 0.
;         (see OFFSET to offset signal)
;
; 	MAX : max value plotted, 
;		every data point larger than MAX takes the same color as MAX
;
; 	MIN : min value plotted, 
;		every data point smaller than MIN takes the same color as MIN
;
;	NESTED: specify that the online file is ordered in the nested scheme
;
; 	NOBAR : if set, no color bar is present
;
;	NOLABELS : if set, no color bar label (min and max) is present
;
;	NOPOSITION : if set, the astronomical location of the map
;	        central point is not indicated
;               ** gnomview only **
;
;       NO_DIPOLE: if set to 1 (and GAL_CUT is not set) 
;                the best fit monopole *and* dipole over all valid pixels are removed
;                * if GAL_CUT is set to b>0, the best monopole and dipole fit is done on all valid
;                pixels with |galactic latitude|>b (in deg) and is removed from all
;                pixels
;             can not be used together with NO_MONOPOLE 
;             (see: GAL_CUT, NO_MONOPOLE)
;               ** mollview only **
;
;       NO_MONOPOLE: if set to 1 (and GAL_CUT is not set) 
;                the best fit monopole over all valid pixels is removed
;                * if GAL_CUT is set to b>0, the best monopole fit is done on all valid
;                pixels with |galactic latitude|>b (in deg) and is removed from all
;                pixels
;             can not be used together with NO_DIPOLE 
;             (see: GAL_CUT, NO_DIPOLE)
;               ** mollview only **
;
;       OFFSET: additive offset to apply to data (default = 0)
;               the data plotted is of the form FACTOR*(data + OFFSET)
;               can be used together with LOG
;               see also : FACTOR, LOG
;               Note : does NOT apply to polarization direction or amplitude
;
; 	ONLINE: if set, you can put directly A HEALPIX VECTOR on File (and
;    		without header): useful when the vector is already
;    		available on line, and avoid to have to write it on disk
;    		just to be read by mollview
;		N.B. : the content of file_in is NOT altered in the
;		process
;               **  can not be used with /SAVE  **
;
;       OUTLINE : structure containing the coordinates of one (or several) outline(s) to
;               be overplot.
;           the structure, should contain for each outline, the following fields : 
;           - 'COORD' coordinate system (either, 'C', 'G', or 'E') of the contour
;           - 'RA'  or longitude coordinates (array)
;           - 'DEC' or lattitude coordinates (array of same size)
;           - 'LINESTYLE' : +2 : black dashes
;                           +1 : black dots
;                            0 : black solid [default]
;                           -1 : black dots on white background
;                           -2 : black dashes on white background
;
;	PNG : string containing the name of a .PNG output
;	      if set to 0 or not set : no .PNG done
;	      if set to 1            : output the plot in plot_gnomic.png
;	      if set to a file name  : output the plot in that file 
;             (see also : CROP, GIF, PS and PREVIEW)
;
;       POLARIZATION: 
;         if set to 0, no polarization information is plotted.
;
;         otherwise, and if the input data contains polarisation information
;             (ie, Stokes parameter Q and U for each pixel)
;
;         if set to 1 
;             the AMPLITUDE P = sqrt( U^2 + Q^2 ) of the polarisation is plotted
;
;         if set to 2 
;             the ANGLE phi = 0.5*ATAN(U/Q) of the polarisation is plotted
;             Note: the angles are color coded with a fixed color table (independant of Colt)
;
;         if set to 3
;             -the temperature is color coded (with a color table defined by Colt)
;             -and the polarisation is overplot as a headless vector
;
;	PREVIEW : if set, there is a 'ghostview' preview of the postscript file (see : PS)
;                    or a 'xv' preview of the gif or png file (see : GIF, PNG)
;
;	PS :  if set to 0 or not set : output to screen
;	      if set to 1            : output the plot in plot_mollweide.ps
;                                                         plot_gnomic.ps
;                                                         plot_cartesian.ps
;	      if set to a file name  : output the plot in that file 
;               (see : PREVIEW)
;
; 	PXSIZE: number of horizontal screen_pixels / postscript_color_dots of the plot
;    		** mollview : default = 800, gnomview and cartview : default = 500 **
;    		(useful for high definition color printer)
;
; 	PYSIZE: number of vertical screen_pixels or postscript_color_dots of the plot
;    		default = PXSIZE
;    		(useful for high definition color printer)
;                ** gnomview only **
;
;       RESO_ARCMIN: resolution of the gnomic map in arcmin
;       (default=1.5)
;                ** gnomview only **
;
; 	ROT :   vector with 1, 2 or 3 elements specifing the rotation angles in DEGREE
;               to apply to the map in the 'output' coordinate system (see coord)
;             = ( lon0, [lat0, rat0]) 
;               lon0 : longitude of the point to be put at the center of the plot
;		       the longitude increases Eastward, ie to the left of the plot 
;                      (unless flip is set)
; 		      =0 by default
;               lat0 : latitude of the point to be put at the center of the plot
; 		      =0 by default
;               rot0 : anti clockwise rotation to apply to the sky around the
;                     center (lon0, lat0) before projecting
;                     =0 by default
;
; 	SAVE: if set, tells that the file is in IDL saveset format, 
;    		the variable saved should be DATA 
;                 ** can not be used with /ONLINE **
;
; 	SUBTITLE : String containing the subtitle to the plot (see TITLEPLOT)
;
; 	TITLEPLOT : String containing the title of the plot, 
;     		if not set the title will be File (see SUBTITLE)
;
;	UNITS : character string containing the units, to be put on the right
;		side of the color bar (see : NOBAR)
;
;	XPOS : The X position on the screen of the lower left corner
;	        of the window, in device coordinate
;
;	YPOS : The Y position on the screen of the lower left corner 
;               of the window, in device coordinate
;
; NOTES
; 	this routine doesn't use the IDL map_set because it is precisely bugged for 
; 	the mollweide projection (at least in IDL 4.0)
;
; SIDE EFFECTS
; 	this routine uses ghostview when PS and PREVIEW are selected 
;	or xv when GIF or PNG and PREVIEW are selected
;
; EXAMPLES
;       ;to plot the signal of the COBE-DMR 4 year map at 53 GHz
;       read_fits_sb, 'dmr_skymap_53a_4yr.fits', dmr53a, /merge  ; read it only one time
;       mollview, dmr53a, /online, 'Sig', /quad
;
;       ;to plot it in Galactic coordinate instead of Ecliptic
;       mollview, drm53a, /online, 'Sig', /quad, coord='g'
;
; COMMONS USED : view_data
;
; PROCEDURES USED: 
;       in the Healpix package :
;	  index_word, read_fits_sb, vec2pix_ring, vec2pix_nest, euler_matrix
;         see  http://www.tac.dk/~healpix
;       it also requires the IDL astro library
;         http://idlastro.gsfc.nasa.gov/homepage.html
;       and the COBE analysis software
;         http://www.gsfc.nasa.gov/astro/cobe/cgis.html
;
; MODIFICATION HISTORY:
; 	October 1997, Eric Hivon, TAC
; 	Nov, 5, 1997,  correction of some bugs for postscript output
; 	13-Nov-97, H. Dole, IAS: save and log keywords
;  	4-Dec-97, H. Dole, IAS: online keyword
; 	16-Dec-97, E.H, TAC: added pxsize, hxsize, subtitle, nobar
;	17-Dec-97, split the loop for projection, added nolabels
;	March-98, E.H. added UNITS keyword
;	April-May-98 E.H. added NESTED_ONLINE, XPOS, YPOS, NOPREVIEW
;       March-99     E.H. Caltech, improved the GIF output
;              modified to deal with structures
;              added Select, COORD, ROT, QUADCUBE  suppressed LON0
;       April-99     E.H. Caltech, improved graticule
;       Nov-99         added flip
;       Feb-00   added rmmonopole and dipole, changed common
;       March-00   changed to no_monopole and no_dipole, changed common
;       Sept-00    added polarisation plotting (Polarization)
;       June-02  : EH, Caltech. Hacked G. Giardino's polview into cartview
;       June-02    partial consolidation of gnomview/mollview/cartview
;-

defsysv, '!healpix', exists = exists
if (exists ne 1) then init_healpix

@viewcom ; define common
data_plot = 0 ; empty common array

loadsky                         ; cgis package routine, define rotation matrices
projection = 'MOLLWEIDE'
routine = 'mollview'

uroutine = strupcase(routine)
if keyword_set(help) then begin
    doc_library,'mollview'
    return
endif

if keyword_set(gif) then begin
    message_gif, code=routine, error=error_gif
    if (error_gif) then return
endif

if (n_params() lt 1 or n_params() gt 2) then begin
    PRINT, 'Wrong number of arguments in '+uroutine
    print,'Syntax : '
    print, uroutine+', File, [Select, ]'
    print,'              [CHARSIZE=, COLT=, COORD=, CROP=, '
    print,'              FACTOR=, FLIP=, GAL_CUT=, GIF=, GRATICULE=, '
    print,'              HELP=, '
    print,'              HIST_EQUAL=, HXSIZE=,'
    print,'              IGRATICULE=,'
    print,'              LOG=, '
    print,'              MAX=, MIN=, NESTED=, NOBAR=, NOLABELS=, '
    print,'              NO_DIPOLE, NO_MONOPLE, '
    print,'              OFFSET=, ONLINE=, OUTLINE=,'
    print,'              PNG=,'
    print,'              POLARIZATION=, PREVIEW=, '
    print,'              PS=, PXSIZE=, PYSIZE=, QUADCUBE= ,'
    print,'              ROT=, SAVE=, '
    print,'              SUBTITLE=, TITLEPLOT=, '
    print,'              UNITS=, XPOS=, YPOS=]'
    print
    print,' Type '+uroutine+', /help '
    print,'   for an extended help'
    return
endif

IF (undefined(file_in)) then begin
    print,routine+': Undefined variable as 1st argument'
    return
endif
do_flip = keyword_set(flip)

if (!D.n_colors lt 4) then begin
    print,' : Sorry ... not enough colors ('+strtrim(string(!d.n_colors),2)+') available'
    return
endif

if (keyword_set(no_monopole) and keyword_set(no_dipole)) then begin
    print,routine+': choose either NO_MONOPOLE or NO_DIPOLE'
    print,'    (removal of best fit monopole only or best fit monopole+dipole)'
    return
endif

polar_type = 0
if keyword_set(polarization) then polar_type = polarization

loaddata, $
  file_in, select_in,$
  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, title_display, sunits, $
  SAVE=save, ONLINE=online, NESTED=nested_online, UNITS=units, COORD=coord, FLIP=flip, $
  ROT=rot, QUADCUBE=quadcube, LOG=log, ERROR=error, $
  POLARIZATION=polarization, FACTOR=factor, OFFSET=offset
if error NE 0 then return

data2moll, $
  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, $
  planmap, Tmax, Tmin, color_bar, planvec, vector_scale, $
  PXSIZE=pxsize, LOG=log, HIST_EQUAL=hist_equal, MAX=max_set, MIN=min_set, FLIP=flip,  $
  NO_DIPOLE=no_dipole, NO_MONOPOLE=no_monopole, UNITS=sunits, DATA_plot = data_plot, GAL_CUT=gal_cut, $
  POLARIZATION=polarization

lg_mollview_proj2out, $
  planmap, Tmax, Tmin, color_bar, 0., title_display, $
  sunits, coord_out, do_rot, eul_mat, planvec, vector_scale, $
  CHARSIZE=charsize, COLT=colt, CROP=crop, GIF = gif, GRATICULE = graticule, $
  HXSIZE=hxsize, NOBAR = nobar, NOLABELS = nolabels, PNG = png, PREVIEW = preview, PS=ps, PXSIZE=pxsize, $
  SUBTITLE = subtitle, TITLEPLOT = titleplot, XPOS = xpos, YPOS = ypos, $
  POLARIZATION=polarization, OUTLINE=outline, /MOLL, FLIP=flip, COORD_IN=coord_in, IGRATICULE=igraticule

w_num = !d.window

return
end

