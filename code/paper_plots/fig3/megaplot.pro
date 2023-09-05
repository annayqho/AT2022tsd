;+
; NAME:
;   megaplot
;
; PURPOSE: 
;   A wrapper for (2D) IDL plotting with potentially asymmetric error bars.
;   Better defaults and more flexible input than the general plot command. 
;
; CALLING SEQUENCE:
;   megaplot, x, y, [...]
;   (see below)
;       
; INPUTS:
;    x, y - Coordinates to plot
;
; INPUT KEYWORDS:
;    noerase - Keep drawing without erasing previous plot (used to add data)
;  Graphical options (can be vectors):
;    psym - Plotting symbol.   0:circle, 3:star, 4:up-triangle, 5:down-triangle,
;                              8:square, 9:diamond, 10:X, 11:4-pt-X-star, 12:left-triangle,
;                              13:right-triangle, 14:4-pt-+-star, 15:pentagon, 16:hexagon,
;                              17:plus, 18:sideways-hourglass, 19:hourglass
;    size         - Symbol size
;    fill         - Fill
;    color        - Color
;    thick        - Thickness of lines
;  Error bars and limit arrows:
;    xerr, yerr   - Error bar lengths (symmetric error bars in data units)
;    xmerr, ymerr - Lower error bar lengths
;    xperr, yperr - Upper error bar lenghts
;    limit        - Indicates y-axis limit: 1=upper, -1=lower
;    xlimit       - Indicates x-axis limit: 1=upper, -1=lower
;    hatlength    - Length of hats on error bars
;    errcolor     - Error bar color
;  More point-specific options:
;    skip      - Ignore this point
;    caption   - Write text next to this data point
;  Options largely just passed to plot:
;    /xlog, /ylog     - Flag to indicate logarithmic plot
;    /xrange, /yrange - Specify plot ranges explicity
;    /xstyle, /ystyle - Set axis behavior using the complicated IDL bit-coding logic
;    charsize         - Text size
;    xtickformat, ytickformat - Axis format
;    xtickv, ytickv           - Tick values
;    xticklen, yticklen       - Tick lengths
;    xminor, yminor           - Number of minor ticks between major ones
;  Secondary axis options:
;    /userightaxis - Draw a right-side y axis
;    rightyrange   - Range of right-side y axis
;    rightystyle   - ystyle code for right y axis
;    /rightylog    - Logarithmic axis
;    rightytitle   - Right axis title
;    rightynormpos - TBD
;    rytickv       - Right axis y tick values
;    ryminor       - Number of right axis y ticks
;  Other:
;    /backwards    - Draw the points starting at the end of the array
;    /randomize    - Draw the points in random order
;
; Almost all graphical options can be specified as vectors or scalars.
;
; Could use _EXTRA = pkey to allow further options, but not enabled yet.
;
;
; COMMENTS:
;    Used as a substitute for the IDL plot command.  Intended to draw logical axes/ticks
;    although doesn't perform well for logarithmic plots on small scales (IDL is extremely
;    finicky about how things are specified on log plots).  Could use some updating
;    of options, e.g to allow fill and lines to be different colors without having to
;    issue the command twice.
;
; Written by D. Perley.  Last modified 29 Jul 2021.  Documentation added 3 July 2022.
;-

; Plots too many points on the axis.  Need to work on a better system for
; increments that are not powers of 10.

pro megaplot, x, y, xerr=xerr, xmerr=xmerr, xperr=xperr, yerr=yerr, ymerr=ymerr, yperr=yperr, psym=psym, size=size, fill=fill, color=color, thick=thick, limit=limit, xlimit=xlimit, skip=skip, caption=caption, xlog=xlog, ylog=ylog, xrange=xrange, yrange=yrange, xstyle=xstyle, ystyle=ystyle, xtitle=xtitle, ytitle=ytitle, title=title, charsize=charsize, xtickformat=xtickformat,ytickformat=ytickformat, xtickv=xtickv, ytickv=ytickv, xticklen=xticklen, yticklen=yticklen, xminor=xminor, yminor=yminor, rightyrange=rightyrange, rightystyle=rightystyle, rightylog=rightylog, rightytitle=rightytitle, rightynormpos=rightynormpos, position=position, noerase=noerase, hatlength=hatlength, errcolor=errcolor, userightaxis=userightaxis, rytickv=rytickv, ryminor=ryminor, backwards=backwards, randomize=randomize

; these might get overwritten so store them now
if n_elements(color) gt 0 then storecolor = color
if n_elements(yperr) gt 0 then storeyperr = yperr
if n_elements(ymerr) gt 0 then storeymerr = ymerr
if n_elements(xperr) gt 0 then storexperr = xperr
if n_elements(xmerr) gt 0 then storexmerr = xmerr
if n_elements(psym) gt 0 then storepsym = psym

if n_elements(thick) gt 0 then storethick = !p.thick
if n_elements(thick) gt 0 then !p.thick = thick

; xticklen appears to be completely broken and it's not clear why.

if n_elements(xstyle) gt 0 then xs = xstyle else xs = 0
if n_elements(xstyle) eq 0 then xstyle = 0

initposition = !p.position

if total(!p.position) eq 0 and n_elements(rightyrange) gt 0 then !p.position = [0.09, 0.09, 0.92, 0.98]
 ; needed to know where to put the right axis (probably a way around this but I'm lazy)

if n_elements(size) gt 0 then symsize = size
if n_elements(position) gt 0 then !p.position = position
if n_elements(hatlength) gt 0 then hatlength = !d.x_vsize / 100 * hatlength
if keyword_set(noerase) eq 0 then begin
  ; set up the axes

  if n_elements(x) eq 0 then x = [0,1]
  if n_elements(y) eq 0 then y = [0,1]

  if n_elements(xerr) eq n_elements(x) then xbuff = xerr/4. else xbuff = x/20.
  if n_elements(yerr) eq n_elements(y) then ybuff = yerr/4. else ybuff = y/20.
    ;above could be refined to take into account x/ystyle and log/linear
  
  minx = n_elements(xrange) eq 2 ? xrange[0] : min(x-xbuff)
  maxx = n_elements(xrange) eq 2 ? xrange[1] : max(x+xbuff)
  miny = n_elements(yrange) eq 2 ? yrange[0] : min(y-ybuff)
  maxy = n_elements(yrange) eq 2 ? yrange[1] : max(y+ybuff)

  ; plot's axis labeling defaults suck, so let's make our own defaults.
  if n_elements(xtickformat) eq 0 then begin
    if keyword_set(xlog) then begin
       if maxx lt 1e+6 and minx gt -1e+6 then xtickformat =''
    endif else begin
       if maxx lt 1e+6 and minx gt 1 then xtickformat = '(I)'
    endelse
  endif
  if n_elements(ytickformat) eq 0 then begin
    if keyword_set(ylog) then begin
       if maxy lt 1e+6 and miny gt 1 then ytickformat ='(I)'
       if miny lt 1 and  miny gt 0.1 then ytickformat = '(F5.1)'
    endif else begin
       if maxy lt 1e+6 and miny gt 1 then ytickformat = '(I)'
       if abs(maxy-miny) le 3.5 and abs(maxy-miny) gt 0.1 then ytickformat = '(F5.1)'
    endelse
  endif
  if keyword_set(xlog) and n_elements(xtickv) eq 0 then begin
    if minx le 0 then begin
      xpos = where(x gt 0., ct)
      if ct gt 2 then minx = min(x[xpos]) else begin
        print, 'Cannot plot non-positive data logarithmically'
        return
      endelse
    endif
    xminlog = floor(alog10(minx < maxx))
    xmaxlog = ceil(alog10(maxx > minx))
    if xmaxlog - xminlog eq 1 then cnt = [1,2,3,4,5,6,7,8,9]
    if xmaxlog - xminlog eq 2 then cnt = [1,2,4,6] ;formerly a 1.5 and an 8
    ;if xmaxlog - xminlog eq 2 then cnt = [1,2,5] ; formerly 3
    if xmaxlog - xminlog ge 3 then cnt = [1]
    tickv = [0]
    for i = xminlog, xmaxlog+1 do begin
      tickv = [tickv, cnt*(10.^i)]
    endfor
    xtickv = tickv[1:*]
  endif
  if n_elements(xtickv) gt 0 then $
    xs = xs or 4 ;sets appropriate bit to suppress axis
  if n_elements(ytickv) gt 0 then begin
     userytickv=ytickv
     yts = n_elements(userytickv)
  endif
  if keyword_set(ylog) and n_elements(ytickv) eq 0 then begin
    if miny le 0 then begin
      ypos = where(y gt 0., ct)
      if ct gt 2 then miny = min(y[ypos]) else begin
        print, 'Cannot plot non-positive data logarithmically'
        return
      endelse
    endif
    yminlog = floor(alog10(miny))
    ymaxlog = ceil(alog10(maxy))
    if ymaxlog - yminlog eq 0 then cnt = [1,2,3,4,5,6,7,8,9]
    if ymaxlog - yminlog gt 0 and ymaxlog - yminlog lt 3 then cnt = [1,2,4,6,8]
    if ymaxlog - yminlog ge 3 then cnt = [1]  ;[1,2,5]
    ;if ymaxlog - yminlog ge 4 then cnt = [1]
    tickv = [0]
    for i = yminlog, ymaxlog+1 do begin
      tickv = [tickv, cnt*(10.^i)]
    endfor
    ytickv = tickv[1:*]
  endif

 if n_elements(rightyrange) gt 0 then begin
    if ystyle lt 8 then ys = ystyle + 8 else ys = ystyle
 endif else begin
    if n_elements(ystyle) gt 0 then ys = ystyle
 endelse
 if n_elements(userightaxis) eq 0 then userightaxis=0
 if userightaxis then ys = ys + 4
 ;if n_elements(xtickv) gt 0 then begin
   ;xts = 1 
   ;xtf = '(A1)'
 ;endif else begin
 ; if n_elements(xtickformat) gt 0 then xtf = xtickformat 
 ;endelse
 if n_elements(xtickformat) gt 0 then xtf = xtickformat 


 ;ys = ystyle
 if n_elements(ytickv) gt 1 then begin
   ;;;yts = 1
   ytf = '(A1)' ;new commentout                           ; no clue why this code block exists...?
   ;if n_elements(ytickformat) eq 0 then ytf = '' else ytf = ytickformat
 endif else begin
  if n_elements(ytickformat) gt 0 then ytf = ytickformat 
  if n_elements(ytitle) gt 0 then ytt = ytitle
 endelse

 if n_elements(xminor) eq 0 then begin
    xminor = 10
    if keyword_set(xlog) then begin
      if xmaxlog - xminlog ge 4 then xminor = 9 else xminor = 0
      ; decade intervals:  [1], 2,   3,   4,   5,   6,   7,   8,   9,   [10]
      ; others:            [1], 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, [2] NOT ANYMORE
    endif
 endif

 plot, [minx, maxx], [miny, maxy], xrange=xrange, yrange=yrange, $
  xstyle=xs, ystyle=ys, xtickformat=xtf, ytickformat=ytf, xlog=xlog, ylog=ylog,$
  title=title, xtitle=xtitle, ytitle=ytt, xticks=xts, yticks=yts, xticklen=xticklen, yticklen=yticklen, $
  /nodata, /ynozero, xminor=xminor, yminor=yminor, charsize=charsize, ytickv=userytickv

 if n_elements(xtickv) gt 1 then begin
   if ((xstyle and 4) eq 0) and ((xs and 4) ne 0) then begin 
     ;draw bottom x-axis (if axes not explicitly suppressed at top level and not already drawn)
     if keyword_set(xlog) then $  ; minor ticks
         axis, xaxis=0, xtickformat='(A1)', /noerase, xlog=xlog, xrange=xrange, xstyle=1, xminor=9 
     axis, xticks=n_elements(xtickv), xtickv=xtickv, xtickformat=xtickformat, xticklen=xticklen, $
           xaxis=0, xlog=xlog, /noerase, /save, xtitle=xtitle, xminor=xminor, charsize=charsize
   endif
   if (((xstyle and 4) eq 0) and ((xs and 4) ne 0)) and ((xstyle and 8) eq 0) then begin 
     ;draw top x-axis (if axes, incl. top axis, is not explicitly suppressed at top level)
     if keyword_set(xlog) then $   ; minor ticks
         axis, xaxis=1, xtickformat='(A1)', /noerase, xlog=xlog, xrange=xrange, xticklen=xticklen, xstyle=1, xminor=9
     axis, xticks=n_elements(xtickv), xtickv=xtickv, xtickformat='(A1)', xticklen=xticklen, $
           xaxis=1, xlog=xlog, /noerase, /save, xminor=xminor
   endif
 endif

 if n_elements(ytickv) gt 1 then begin
   if n_elements(ytickformat) gt 0 then begin 
      ytf = ytickformat
      ;if ytf eq '(A1)' and ytickformat ne '(A1)' then ytf = ''
   endif else begin
      ytf = ''
   endelse
   axis, yticks=n_elements(ytickv), ytickv=ytickv, ytickformat=ytf, $
         yaxis=0, ylog=ylog, /noerase, /save, ytitle=ytitle, charsize=charsize
 endif

 if n_elements(rightyrange) gt 1 then begin
   if n_elements(ryminor) eq 0 then ryminor = 9
   if n_elements(rytickv) eq 0 then begin
   if keyword_set(rightylog) then begin
     yminlog = floor(alog10(rightyrange[0]))
     ymaxlog = ceil(alog10(rightyrange[1]))
     if ymaxlog - yminlog eq 1 then cnt = [1,2,3,4,5,6,7,8,9]
     if ymaxlog - yminlog eq 2 then cnt = [1,2,4,6] ;8
     if ymaxlog - yminlog eq 3 then cnt = [1,2,5]
     if ymaxlog - yminlog ge 4 then cnt = [1]
     tickv = [0]
     for i = yminlog, ymaxlog+1 do begin
       tickv = [tickv, cnt*(10.^i)]
     endfor
     rytickv = tickv[1:*]
     rytickv = rytickv[where(rytickv ge min(rightyrange) and rytickv le max(rightyrange))]
   endif else begin
     rymin = floor(min(rightyrange))
     rymax = ceil(max(rightyrange))
     rrang = rymax - rymin
     inclog = 10.^floor(alog10(rrang / 5.))
     inc =  rrang/5.
     linc = floor(inc/inclog) 
     if linc ge 3 then linc = 5
     inc = linc*inclog
     firsttick = inc*(floor(min(rightyrange)/inc))
     ;rytickv = firsttick + inc*(2+findgen(ceil(rrang/inc)))
     rytickv = firsttick + inc*(1+findgen(ceil(rrang/inc)))
     rightystyle=1
     ryminor = linc
   endelse
   endif
   if userightaxis then begin
     axis, !p.position[0], !p.position[1], /norm,$
         yaxis=0, yrange=rightyrange, $
         ystyle=rightystyle, ylog=rightylog, $
         ytitle=rightytitle, ytickv=rytickv, $
         yminor=ryminor, charsize=charsize ;normally don't necessary want yminor = 9...!!!!!
     rightytitle=''
   endif
   axis, !p.position[2], !p.position[1], /norm,$ ;formerly rightynormpos
         yaxis=1, yrange=rightyrange, $
         ystyle=rightystyle, ylog=rightylog, $
         ytitle=rightytitle, ytickv=rytickv, yticks=n_elements(rytickv)-1, $
         yminor=ryminor ;normally don't necessary want yminor = 9...!!!!!
 endif
 if n_params() eq 0 then return
endif

n = min([n_elements(x), n_elements(y)])

if n_elements(xerr) eq 0 then xerr = replicate(0.,n)
if n_elements(yerr) eq 0 then yerr = replicate(0.,n)
if n_elements(xmerr) eq 0 then xmerr = xerr
if n_elements(xperr) eq 0 then xperr = xerr
if n_elements(ymerr) eq 0 then ymerr = yerr
if n_elements(yperr) eq 0 then yperr = yerr

; some defaults
if n_elements(psym) eq 0 then psym = (ymerr gt 0 and yperr gt 0 and xmerr gt 0 and xperr gt 0)* replicate(-1,n)$
  ; no error bars: default to circles.   error bars: default to no symbol (just bars)
  else if n_elements(psym) eq 1 then psym = replicate(psym, n) $
  else if n_elements(psym) lt n then psym =[psym,99+intarr(n-n_elements(psym))]
if n_elements(symsize) eq 0 then symsize = replicate(1,n) $
  else if n_elements(symsize) eq 1 then symsize = replicate(symsize, n) $
  else if n_elements(symsize) lt n then symsize = [symsize,1+fltarr(n-n_elements(symsize))]
if n_elements(fill) eq 0 then fill = replicate(0,n) $
  else if n_elements(fill) eq 1 then fill = replicate(fill, n) $
  else if n_elements(fill) lt n then fill = [fill,intarr(n-n_elements(fill))]
initcolor = !p.color
if n_elements(color) eq 0 then color = replicate(initcolor,n) $
  else if n_elements(color) eq 1 then color = replicate(color, n) $
  else if n_elements(color) lt n then color=[color,initcolor+intarr(n-n_elements(color))]
if n_elements(limit) eq 0 then limit = replicate(0,n) $
  else if n_elements(limit) eq 1 then limit = replicate(limit, n) $
  else if n_elements(limit) lt n then limit=[limit,intarr(n-n_elements(limit))]
if n_elements(xlimit) eq 0 then xlimit = replicate(0,n) $
  else if n_elements(xlimit) eq 1 then xlimit = replicate(xlimit, n) $
  else if n_elements(xlimit) lt n then xlimit=[xlimit,intarr(n-n_elements(xlimit))]

order = lindgen(n)
if keyword_set(backwards) then order = n-1 - indgen(n)
if keyword_set(randomize) then order = sort(randomn(randomize, n)) ;can seed by setting to a value

if n_elements(errcolor) eq 1 then errcol = intarr(n) + errcolor 
if n_elements(errcolor) ge 2 then errcol = errcolor

; plot the data, point by point
for ii = 0L, n-1 do begin
     i = order[ii]

     if n_elements(skip) gt i then begin
        if skip[i] gt 0 then continue 
     endif
     ;if psym[i] eq 0 then plotsym, 0, symsize[i], fill=fill[i]
     !p.color = color[i]

     if psym[i] ge 0 and psym[i] le 98 then begin
       symplot, x[i], y[i], sym=psym[i], size=symsize[i], fill=fill[i], thick=!p.thick
     endif

     if (ymerr[i] gt 0 or yperr[i] gt 0) and limit[i] eq 1 then begin
       ymerr[i] = 0
       yperr[i] = 0  ;not quite ideal as this can alter input value, but whatever
       if psym[i] lt 0 or psym[i] ge 99 then begin
         usersym, [0,0], [0, -1]*symsize[i]
          oplot, [x[i]], [y[i]], psym=8, thick=!p.thick ; connect arrow to horizontal error bar
       endif
     endif

     ;if psym[i] eq 9 then begin
     ;  usersym, [0,1,0,-1,0]*1.4*symsize[i], [-1,0,1,0,-1]*1.4*symsize[i], fill=fill[i] ; diamond
     ;   oplot, [x[i]], [y[i]], psym=8
     ; endif

     ;there is a bug here: if [x/y]merr = 0 but [x/y]perr > 0 then no error bar will display
     ; semi-fixed now for y only - should really look into simplifying this, but was there
     ; a reason I made it so complicated?

     if n_elements(hatlength) gt 0 then nohat = (hatlength eq 0) else nohat = 0

     if n_elements(errcolor) gt 0 then begin ;totally bizarre that I have to do this
       if xmerr[i] gt 0 and (ymerr[i] gt 0 or yperr[i] gt 0) then begin                 
          if ymerr[i] gt 0 then oploterror, x[i], y[i], xmerr[i], ymerr[i], /lobar, $
             hatlength=hatlength, errcolor=errcol[i], thick=!p.thick, errthick=!p.thick, nohat=nohat*1
          if yperr[i] gt 0 then oploterror, x[i], y[i], xperr[i], yperr[i], /hibar, $
             hatlength=hatlength, errcolor=errcol[i], thick=!p.thick, errthick=!p.thick, nohat=nohat*1
       endif
       if xmerr[i] le 0 and (ymerr[i] gt 0 or yperr[i] gt 0)  then begin
          if ymerr[i] gt 0 then oploterror, x[i], y[i], ymerr[i], /lobar, $
             hatlength=hatlength, errcolor=errcol[i], thick=!p.thick, errthick=!p.thick, nohat=nohat*1
          if yperr[i] gt 0 then oploterror, x[i], y[i], yperr[i], /hibar, $
             hatlength=hatlength, errcolor=errcol[i], thick=!p.thick, errthick=!p.thick, nohat=nohat*1
       endif
       if xmerr[i] gt 0 and (ymerr[i] le 0 and yperr[i] le 0) then begin
          oploterror, x[i], y[i], xmerr[i], replicate(0.,n), /lobar, $
             hatlength=hatlength, errcolor=errcol[i], thick=!p.thick, errthick=!p.thick, nohat=nohat*1
          oploterror, x[i], y[i], xperr[i], replicate(0.,n), /hibar, $
             hatlength=hatlength, errcolor=errcol[i], thick=!p.thick, errthick=!p.thick, nohat=nohat*1
       endif
     endif else begin
       if xmerr[i] gt 0 and (ymerr[i] gt 0 or yperr[i] gt 0)  then begin
          if ymerr[i] gt 0 then oploterror, x[i], y[i], xmerr[i], ymerr[i], /lobar, $
             hatlength=hatlength, thick=!p.thick, errthick=!p.thick, nohat=nohat*1
          if yperr[i] gt 0 then oploterror, x[i], y[i], xperr[i], yperr[i], /hibar, $
             hatlength=hatlength, thick=!p.thick, errthick=!p.thick, nohat=nohat*1
       endif
       if xmerr[i] le 0 and (ymerr[i] gt 0 or yperr[i] gt 0)  then begin
          if ymerr[i] gt 0 then oploterror, x[i], y[i], ymerr[i], /lobar, $
             hatlength=hatlength, thick=!p.thick, errthick=!p.thick, nohat=nohat*1
          if yperr[i] gt 0 then oploterror, x[i], y[i], yperr[i], /hibar, $
             hatlength=hatlength, thick=!p.thick, errthick=!p.thick, nohat=nohat*1
       endif
       if xmerr[i] gt 0 and (ymerr[i] le 0 and yperr[i] le 0)  then begin
          oploterror, x[i], y[i], xmerr[i], replicate(0.,n), /lobar, $
             hatlength=hatlength, thick=!p.thick, errthick=!p.thick, nohat=nohat*1
          oploterror, x[i], y[i], xperr[i], replicate(0.,n), /hibar, $
             hatlength=hatlength, thick=!p.thick, errthick=!p.thick, nohat=nohat*1
       endif
     endelse 
      
     if limit[i] then begin
       if limit[i] eq 1 then begin
         if psym[i] gt -1 then $
           usersym, [0,0,-1,0,1]*symsize[i], [-1,-4,-3,-4,-3]*symsize[i], thick=!p.thick $
         else $
           usersym, [0,0,-1,0,1]*symsize[i], [0, -4,-3,-4,-3]*symsize[i], thick=!p.thick
         oplot, [x[i]], [y[i]], psym=8 ; arrow
       endif
       if limit[i] eq -1 then begin
         if psym[i] gt -1 then $
           usersym, [0,0,-1,0,1]*symsize[i], [1,4,3,4,3]*symsize[i], thick=!p.thick $
         else $
           usersym, [0,0,-1,0,1]*symsize[i], [0,4,3,4,3]*symsize[i], thick=!p.thick
         oplot, [x[i]], [y[i]], psym=8 ; arrow
       endif
     endif
     if xlimit[i] then begin
       if xlimit[i] eq 1 then begin
         if psym[i] gt -1 then $
           usersym, [-1,-4,-3,-4,-3]*symsize[i], [0,0,-1,0,1]*symsize[i], thick=!p.thick  $
         else $
           usersym, [0, -4,-3,-4,-3]*symsize[i], [0,0,-1,0,1]*symsize[i], thick=!p.thick
         oplot, [x[i]], [y[i]], psym=8 ; arrow
       endif
       if xlimit[i] eq -1 or xlimit[i] eq 255 then begin  ; 255 is unsigned -1
         if psym[i] gt -1 then $
           usersym, [1,4,3,4,3]*symsize[i], [0,0,-1,0,1]*symsize[i], thick=!p.thick $
         else $
           usersym, [0,4,3,4,3]*symsize[i], [0,0,-1,0,1]*symsize[i], thick=!p.thick
         oplot, [x[i]], [y[i]], psym=8 ; arrow
       endif
     endif

     if n_elements(caption) gt i then begin

       ;data labels
 
       if keyword_set(xlog) eq 0 then begin
         if n_elements(xrange) ge 2 then xrangesize = xrange[1]-xrange[0] $
                                    else xrangesize = maxx-minx
         xcap = x[i] - xrangesize*0.03
       endif else begin
         if n_elements(xrange) ge 2 then xrangesize = xrange[1]/xrange[0] $
                                    else xrangesize = maxx/minx
         xcap = x[i] / (xrangesize^0.03)
       endelse
       if keyword_set(ylog) eq 0 then begin
         if n_elements(yrange) ge 2 then yrangesize = yrange[1]-yrange[0] $
                                    else yrangesize = maxy-miny
         ycap = y[i] + yrangesize*0.02
       endif else begin
         if n_elements(yrange) ge 2 then yrangesize = yrange[1]/yrange[0] $
                                    else yrangesize = maxy/miny
         ycap = y[i] * (yrangesize^0.02)
       endelse
       xyouts, /data, xcap, ycap, caption[i]
     endif

  ;oploterror, x[i], y[i], xerr[i], yerr[i], psym=psym[i], size=symsize[i];;;
endfor
!p.color = initcolor
!p.position = initposition
if n_elements(storecolor) gt 0 then color = storecolor ; else delvarx, color
if n_elements(storeyperr) gt 0 then yperr = storeyperr
if n_elements(storeymerr) gt 0 then ymerr = storeymerr
if n_elements(storexperr) gt 0 then xperr = storexperr
if n_elements(storexmerr) gt 0 then xmerr = storexmerr
if n_elements(storepsym) gt 0 then psym = storepsym
if n_elements(storethick) gt 0 then !p.thick = storethick

end
