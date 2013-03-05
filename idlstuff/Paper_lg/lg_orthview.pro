pro lg_orthview, file_in, select_in, $
CHARSIZE=charsize, COLT = colt, COORD = coord, CROP = crop, $
FACTOR=factor, FLIP=flip, $
GAL_CUT=gal_cut, GIF = gif, GRATICULE = graticule, $
HALF_SKY = half_sky, HELP = help, $
HIST_EQUAL = hist_equal, HXSIZE = hxsize, $
IGRATICULE=igraticule, $
LOG = log, $
MAX = max_set, MIN = min_set, $
NESTED = nested_online, NOBAR = nobar, NOLABELS = nolabels, NO_DIPOLE=no_dipole, NO_MONOPOLE=no_monopole, $
OFFSET=offset, ONLINE = online, OUTLINE=outline, $
PNG=png, POLARIZATION=polarization, PREVIEW = preview, PS = ps, PXSIZE = pxsize, $
QUADCUBE = quadcube, $
ROT = rot, $
SAVE = save, SUBTITLE = subtitle, $
TITLEPLOT = titleplot, $
UNITS = units, XPOS = xpos, YPOS = ypos

;+
; for extended description see mollview or the paper documentation
;-


defsysv, '!healpix', exists = exists
if (exists ne 1) then init_healpix

@viewcom ; define common
data_plot = 0 ; empty common array

loadsky                         ; cgis package routine, define rotation matrices
projection = 'ORTHOGRAPHIC'
routine = 'orthview'

uroutine = strupcase(routine)
if keyword_set(help) then begin
    doc_library,'orthview'
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
    print,'              FLIP=, GAL_CUT=, GIF=, GRATICULE=, '
    print,'              HALF_SKY=,     HELP=, '
    print,'              HIST_EQUAL=, HXSIZE=,  '
    print,'              IGRATICULE=,'
    print,'              LOG=, '
    print,'              MAX=, MIN=, NESTED=, NOBAR=, NOLABELS=, '
    print,'              OFFSET=, ONLINE=, OUTLINE=,'
    print,'              PNG=,'
    print,'              POLARIZATION=polarization, PREVIEW=, '
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

do_fullsky = 1 - keyword_set(half_sky)

loaddata, $
  file_in, select_in,$
  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, title_display, sunits, $
  SAVE=save, ONLINE=online, NESTED=nested_online, UNITS=units, COORD=coord, FLIP=flip, $
  ROT=rot, QUADCUBE=quadcube, LOG=log, ERROR=error, $
  POLARIZATION=polarization, FACTOR=factor, OFFSET=offset
if error NE 0 then return

data2orth, $
  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, $
  planmap, Tmax, Tmin, color_bar, planvec, vector_scale, $
  PXSIZE=pxsize, LOG=log, HIST_EQUAL=hist_equal, MAX=max_set, MIN=min_set, FLIP=flip,  $
  NO_DIPOLE=no_dipole, NO_MONOPOLE=no_monopole, UNITS=sunits, DATA_plot = data_plot, GAL_CUT=gal_cut, $
  POLARIZATION=polarization, HALF_SKY=half_sky

lg_mollview_proj2out, $
  planmap, Tmax, Tmin, color_bar, 0., title_display, $
  sunits, coord_out, do_rot, eul_mat, planvec, vector_scale, $
  CHARSIZE=charsize, COLT=colt, CROP=crop, GIF = gif, GRATICULE = graticule, $
  HXSIZE=hxsize, NOBAR = nobar, NOLABELS = nolabels, PNG = png, PREVIEW = preview, PS=ps, PXSIZE=pxsize, $
  SUBTITLE = subtitle, TITLEPLOT = titleplot, XPOS = xpos, YPOS = ypos, $
  POLARIZATION=polarization, OUTLINE=outline, /ORTH, FLIP=flip, HALF_SKY=half_sky, COORD_IN=coord_in, IGRATICULE=igraticule

w_num = !d.window

return
end

