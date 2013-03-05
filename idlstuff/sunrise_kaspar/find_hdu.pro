function find_HDU, file, name
i = 0
while 1 do begin
    if (strlen (name) lt 8) then begin
        for j=1, 8-strlen (name) do name = name + ' '
    end
    h = headfits (file, exten = i)
    n =fxpar (h, 'extname')
    ;;print, '"',n,'"'
    if (strcmp (n, name,/fold_case)) then return, i else i = i+ 1
end
end


