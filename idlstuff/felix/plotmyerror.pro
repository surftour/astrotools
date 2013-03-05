PRO plotmyerror,xx,x,y,e,thick=thick,color=color,symsize=symsize,linestyle=linestyle

    yy = interpol(y,x,[xx,xx])
    yy = yy(0)
    ee = interpol(e,x,[xx,xx])
    ee = ee(0)
    ; plot a dot, a line and small end-ticks
    ;plotsym,0,/fill
    usersym,1.5*[-1,-1,1,1,-1],1.5*[-1,1,1,-1,-1],thick=4.0
    oplot,[xx,xx],[yy,yy],psym=8,symsize=symsize,color=color,thick=thick
    myoploterror,[xx,xx],[yy,yy],[ee,ee],thick=thick,errcolor=color,errstyle=linestyle,hatlength=2.0*!D.X_VSIZE / 100
    ;oplot,[xx,xx],[yy-ee,yy+ee],thick=thick,color=color,linestyle=linestyle


END
