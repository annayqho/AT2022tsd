;+
; NAME:
;   symplot
;
; PURPOSE: 
;   Specify (and then use) a variety of custom symbols for IDL plotting. 
;   Used by megaplot for custom symbols  
;
; CALLING SEQUENCE:
;   symplot, x, y, err=err, sym=sym, fill=fill, size=size, thick=thick, norm=norm
;       
; INPUTS:
;   x, y  - Coordinates
;   err   - Basic error bar lengths
;   sym   - Symbol(s), see codes below
;   /fill - Filled symbol
;   size  - Symbol size
;   thick - Line thickness
;   /norm - Interpret coordinates as normalized plot coordinates
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       
; OPTIONAL OUTPUT KEYWORD PARAMETERS:
;       
; OUTPUTS: 
;
; Written by D. Perley.  Last modified 12 Jul 2018.  Documentation added 3 July 2022.
;-

pro symplot, x, y, err=err, sym=sym, fill=fill, size=size, thick=thick, norm=norm
                               ;sym=-1 means use current usersym
  ; fill and size must be scalars.  sym can be a vector.

  ax = [x]
  ay = [y]  
  if (keyword_set(size) eq 0) then size = 1
  if n_elements(sym) eq 1 and n_elements(ax) gt 1 then $
    sym = replicate(sym,n_elements(ax))


  if n_elements(sym) eq 0 then begin
     ;oplot, ax, ay, psym=8
     if keyword_set(norm) then plots, ax, ay, psym=8, norm=norm, thick=thick $
                          else oplot, ax, ay, psym=8, thick=thick
  endif else begin
    for s = -1, 19 do begin
      ws = where(sym eq s, cw)
      if cw eq 0 then continue
      if s ge 0 and s le 8 then plotsym, s, size, fill=fill, thick=thick  ; except 4/5..?
          ; IDLASTRO plotsym:
          ;     0 - circle
          ;     1 - downward arrow (upper limit)
          ;     2 - upward arrow (lower limt)
          ;     3 - 5 pointed star
          ;     4 - triangle                [overridden below for better sizing]
          ;     5 - upside down triangle    ["]
          ;     6 - left pointing arrow (left limit)
          ;     7 - right pointing arrow (right limit)
          ;     8 - square
          ; New options:
          ;     9  - diamond
          ;     10 - blocky X
          ;     11 - pointy X
          ;     12 - left-pointing triangle
          ;     13 - right-pointing triangle
          ;     14 - four-pointed star
          ;     15 - pentagon
          ;     16 - hexagon
          ;     17 - plus symbol
          ;     18 - bowtie
          ;     19 - hourglass
      sc=1.3*size
      if s eq 4 then $  ; down triangle
         usersym, [-1,0,1,-1]*sc,[-0.577,1.154,-0.577,-0.577]*sc, fill=fill, thick=thick
      if s eq 5 then $  ; up triangle
         usersym, [-1,0,1,-1]*sc,[0.577,-1.154,0.577,0.577]*sc, fill=fill, thick=thick

      if s eq 9 then $  ;diamond
        usersym, [0,1,0,-1,0]*sc, [-1,0,1,0,-1]*sc, fill=fill, thick=thick
      if s eq 10 then $ ; right-angles X
        usersym, [0,0.5,1,0.5,1,0.5,0,-0.5,-1,-0.5,-1,-0.5,0]*sc,$
                 [0.5,1,0.5,0,-0.5,-1,-0.5,-1,-0.5,0,0.5,1,0.5]*sc,$
                 fill=fill, thick=thick
      if s eq 11 then $ ; pointy X
         usersym, [-1,0  ,1,0.4,1,0   ,-1,-0.4,-1]*sc, $
                  [ 1,0.4,1,0, -1,-0.4,-1,0,    1]*sc, fill=fill, thick=thick
      if s eq 12 then $ ; left triangle
         usersym, [0.577,-1.154,0.577,0.577]*sc, [-1,0,1,-1]*sc, fill=fill, thick=thick
        ;usersym, [-0.452,-0.452,0.872,-0.452]*sc, [-0.75,0.75,0,-0.75]*sc,fill=fill
      if s eq 13 then $ ; right triangle
         usersym, [-0.577,1.154,-0.577,-0.577]*sc, [-1,0,1,-1]*sc, fill=fill, thick=thick
        ;usersym, [0.4,0.4,-0.6,0.4]*sc, [-0.5,0.5,0,-0.5]*sc,fill=fill
      ;oplot, ax[ws], ay[ws], psym=8
      if s eq 14 then $ ; four-point star
         usersym, [-1.1,-0.3,0,0.3,1.1,0.2,0,-0.3,-1.1]*sc,$
                  [0,0.3,1.1,0.3,0,-0.3,-1.1,-0.3,0]*sc, fill=fill, thick=thick
      if s eq 15 then $ ; pentagon
         usersym, [0, 0.951, 0.588,-0.588,-0.951, 0]*sc, $
                  [1, 0.309,-0.809,-0.809, 0.309, 1]*sc, fill=fill, thick=thick
      if s eq 16 then $ ;hexagon
        usersym, [0,0.866,0.866,0,-0.866,-0.866,0]*sc, $
                 [1,0.5,-0.5,-1,-0.5,0.5,1]*sc,fill=fill, thick=thick
      if s eq 17 then $ ; plus sign
        usersym, [-0.3,-0.3,0.3,0.3,1  ,1,    0.3,0.3,-0.3,-0.3,-1,-1,-0.3]*sc, $
                 [ 0.3, 1,   1, 0.3,0.3,-0.3,-0.3,-1, -1, -0.3,-0.3,0.3,0.3]*sc,$
                  fill=fill, thick=thick
      if s eq 18 then $ ; bowtie
        usersym, [-1,0,  1, 1, 0, -1,-1]*0.8*sc, $
                 [ 1,0.35,1,-1,-0.35,-1,1]*sc,$
                  fill=fill, thick=thick
      if s eq 19 then $ ; hourglass
        usersym, [ 1,0.35,1,-1,-0.35,-1,1]*sc, $
                 [-1,0,  1, 1, 0, -1,-1]*0.8*sc,$
                  fill=fill, thick=thick

      if keyword_set(norm) then plots, ax[ws], ay[ws], psym=8, thick=thick, norm=norm $
                                else oplot, ax[ws], ay[ws], psym=8, thick=thick
    endfor
  endelse

  if (keyword_set(err)) then begin
    errplot, ax, ay-err, ay+err 
  endif

end
