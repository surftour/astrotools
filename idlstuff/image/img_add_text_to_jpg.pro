
pro img_add_text_to_jpg, inputfile, outputfile, texttoadd, x0, y0

font= "Courier-Bold"
;font= "arial"

ptsize= 30

color= "black"

cmd = "convert"
cmd = cmd + " -font "+font
cmd = cmd + " -antialias"
cmd = cmd + " -pointsize "+strcompress(string(ptsize),/remove_all)
cmd = cmd + " -fill "+color
cmd = cmd + " -draw 'text "
cmd = cmd + strcompress(string(x0),/remove_all) + "," + strcompress(string(y0),/remove_all) + " "
cmd = cmd + """
cmd = cmd + texttoadd
cmd = cmd + """
cmd = cmd + "'"
cmd = cmd + ' -quality 95 '
cmd = cmd + ' -colorspace CMY '
cmd = cmd + " "+inputfile+" "+outputfile+" "
print, cmd
spawn, cmd




end



