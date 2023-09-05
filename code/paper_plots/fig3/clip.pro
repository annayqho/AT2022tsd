;+
; NAME:
;   clip
;
; PURPOSE: 
;   Compact expression for converting a number to string format and/or removing whitespace
;   (also works on vectors).
;
; CALLING SEQUENCE:
;   out = clip(in, nchar)
;       
; INPUTS:
;   in    - A number (int, float, etc) or string; or vector of numbers/strings
;   nchar - Number of characters in output string(s)
;       
; OUTPUTS: 
;   out   - Output string or string array.
;
; See also rclip() to clip a string from the right side.
;
; Last modified 30 May 2014.  Documentation added 3 July 2022.
;-

function clip, strin, maxlen
  
  if size(strin, /type) ne 1 then str = string(strin) else str = string([fix(strin)]) ; if it's a byte, weird stuff happens
  if n_params() eq 1 then return, strtrim(str,2)

  for n = 0, n_elements(strin)-1 do begin
     str[n] = strmid(strtrim(str[n],2),0,maxlen)
     curlen = strlen(str[n])
     if curlen lt maxlen then $
       for c = 0, (maxlen-curlen-1) do $
         str[n] = str[n] + ' '
  endfor

  return, str
end
