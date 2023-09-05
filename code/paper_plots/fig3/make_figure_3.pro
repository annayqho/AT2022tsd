pro openimplot, plotname, xsize=xsize, ysize=ysize

   psopen, plotname, /encaps, xsize=xsize, ysize=ysize, /inches
   device, /symbol, font_index = 4
   device, /helvetica, font_index = 17
   device, /helvetica, /bold, font_index = 18
   device, /isolatin, font_index = 19
   device, /color
   colors = indgen(256)##[1,1,1] ;, [255,0,0], [0,255,0], [0,0,255]]
   colors[*,251] = [219,67,37]  ;[225,95,69]  ;[255,90,90] ;[219,67,37]
   colors[*,252] = [87,196,164] ;[0,255,0]
   colors[*,253] = [0,0,255]
   colors[*,254] = [255,255,0]
   colors[*,255] = [255,255,255]
   tvlct, transpose(colors)

   !p.font = 0
   !x.thick = 2
   !y.thick = 2
   !p.thick = 2

end

function abmag2ujy, m
  return, 3631e6 * 10.^(-m/2.5)
end

pro addmagaxis, xmax, tickscale, charsize=charsize, minmag=minmag, maxmag=maxmag
  if n_elements(minmag) eq 0 then minmag = 20
  if n_elements(maxmag) eq 0 then maxmag = 23
  axis, yaxis=1, yticklen=1e-19, ytitle='', ytickformat='(A1)'
  for m = 23.,17.,-.1 do oplot, [-0.002,0]*tickscale+xmax, abmag2ujy(m)+[0,0]
  for m = 25, 17, -1  do oplot, [-0.004,0]*tickscale+xmax, abmag2ujy(m)+[0,0]
  for m = maxmag, minmag, -1  do xyouts, xmax+0.0005, abmag2ujy(m)-0.45, clip(m), charsize=charsize
end

pro addimage, im, h, rac, decc, boxsize, minval=minval, maxval=maxval, invert=invert
   s = size(im)
   nx = s[1]
   ny = s[2]

   cd11 = sxpar(h,'CD1_1')
   cd12 = sxpar(h,'CD1_2')
   cd21 = sxpar(h,'CD2_1')
   cd22 = sxpar(h,'CD2_2')
   cdm = [[cd11,cd12],[cd21,cd22]]
   dx_corner = [-1,-1,1,1]*1.01/2.
   dy_corner = [-1,1,1,-1]*1.01/2.
   dgx_corner = fltarr(4)
   dgy_corner = fltarr(4)
   for i = 0, n_elements(dx_corner)-1 do begin
     dg = cdm # [dx_corner[i], dy_corner[i]]
     dgx_corner[i] = dg[0]*3600./boxsize
     dgy_corner[i] = dg[1]*3600./boxsize
   endfor

   plot, [0], [0], xrange=[0,1], yrange=[0,1], xsty=5, ysty=5, xtickformat='(A1)', ytickformat='(A1)'
   
   if keyword_set(invert) eq 0 then begin
     u = 0    ; no signal
     v = 250  ; saturated
     d = +1
   endif else begin
     u = 250
     v = 0
     d = -1
   endelse
   
   polyfill, [0,1,1,0,0],[0,0,1,1,0], color=u
   
   xscale = 3600. * sqrt(cd11^2 + cd21^2) ; size of a pixel in arcsec
   yscale = 3600. * sqrt(cd12^2 + cd22^2)
   pixscale = sqrt(xscale*yscale)

   dx = xscale/boxsize ; size of a pixel in graphical coordinates
   dy = yscale/boxsize
   
   adxy, h, rac, decc, xc, yc  ; center of display (image pixel)
   x0 = floor(xc - 1.44*(boxsize/xscale))
   x1 =  ceil(xc + 1.44*(boxsize/xscale))
   y0 = floor(yc - 1.44*(boxsize/yscale))
   y1 =  ceil(yc + 1.44*(boxsize/yscale))
   x0 = (x0 > 0)
   y0 = (y0 > 0)
   x1 = (x1 < nx-1)
   y1 = (y1 < ny-1)
   if x0 gt nx-1 then print, 'Outside image'
   if y0 gt ny-1 then print, 'Outside image'
   if x1 lt 0 then print, 'Outside image'
   if y1 lt 0 then print, 'Outside image'
    
   imcut = im[x0:x1,y0:y1] ; cutout around the target
   if n_elements(minval) eq 0 then minval = min(imcut)
   if n_elements(maxval) eq 0 then maxval = max(imcut)

   imsc = (imcut-minval)/(maxval-minval) ; normalized
   imcol = (fix(imsc*250) > 1) < 250
         
   cosdec = cos(!pi*decc/180.)
   
   npix = s[4]
   xx = indgen(x1-x0+1,y1-y0+1)
   yy = indgen(x1-x0+1,y1-y0+1)
   for x = x0, x1 do $
     xx[x-x0,*] = x
   for y = y0, y1 do $
     yy[*,y-y0] = y
   xyad, h, xx, yy, rapix, decpix
   xyad, h, xx+0.51, yy+0.51, rapix1, decpix1
   xyad, h, xx+0.51, yy-0.51, rapix2, decpix2
   xyad, h, xx-0.51, yy-0.51, rapix3, decpix3
   xyad, h, xx-0.51, yy+0.51, rapix4, decpix4
   
   for x = x0, x1 do begin  ; x and y are in original image system
   for y = y0, y1 do begin         
     ra = rapix[x-x0,y-y0]
     dec = decpix[x-x0,y-y0]
     dra_arcsec = (ra - rac)*3600.*cosdec
     ddec_arcsec = (dec - decc)*3600.
     g_xc = 0.5 - dra_arcsec*1./boxsize   ; position of center of pixel in graphical coordinate
     g_yc = 0.5 + ddec_arcsec*1./boxsize
     if g_xc lt -0.1 or g_xc gt 1.1 then continue
     if g_yc lt -0.1 or g_yc gt 1.1 then continue
     racorners = [rapix1[x-x0,y-y0], rapix2[x-x0,y-y0], rapix3[x-x0,y-y0], rapix4[x-x0,y-y0]]
     deccorners = [decpix1[x-x0,y-y0], decpix2[x-x0,y-y0], decpix3[x-x0,y-y0], decpix4[x-x0,y-y0]]
     dra_corners = (racorners - rac)*3600.*cosdec
     ddec_corners = (deccorners - decc)*3600.
     xcorners = 0.5 - dra_corners*1./boxsize
     ycorners = 0.5 + ddec_corners*1./boxsize
     if 0 then begin
      w = where(xcorners lt 0 or xcorners gt 1, ctx, complement=wx)
      w = where(ycorners lt 0 or ycorners gt 1, cty, complement=wy)
      if ctx gt 0 or cty gt 0 then print, ctx, cty
      if ctx eq 0 and cty eq 1 then begin
        xcorners = xcorners[wy]
        ycorners = ycorners[wy]
        print, ycorners
      endif
      if ctx eq 1 and cty eq 0 then begin
        xcorners = xcorners[wx]
        ycorners = ycorners[wx]     
      endif
     endif
     if max(xcorners) le 0  or min(xcorners) gt 1 then continue
     if max(ycorners) le 0 or  min(ycorners) gt 1 then continue
     polyfill, [xcorners, xcorners[0]], [ycorners, ycorners[0]], color=u+d*imcol[x-x0,y-y0]
   endfor
   endfor

   m = 1.5*xscale/boxsize
   ; lazy way of dealing with subpixels extending beyond edges: just paint over them.
   polyfill, [-m,0,0,-m],[-m,-m,1+m,1+m], color=255
   polyfill, [-m,-m,1+m,1+m],[-m,0,0,-m], color=255
   polyfill, [1+m,1,1,1+m],[-m,-m,1+m,1+m], color=255
   polyfill, [-m,-m,1+m,1+m],[1+m,1,1,1+m], color=255

end

pro make_figure_3

  readcol, 'ultraspeclc.txt', mjd, filt, ee, f, pm, ef, x, exp, filename, $
       format='d,a,a,f,a,f,a,f,a', /silent, comment='#'

  seq1 = 135+indgen(11)*2
  seq2 = max(seq1)+(1+indgen(9))*3
  seq = [seq1,seq2]

  xsize = 7.5
  ysize = 3.5
  
  fontsize = 6.8
  cs = fontsize/12. * 1.06  ; = 0.6 for 6.8 font
    
  plotname = 'lcim_ULTRASPEC_r.eps'
   
  rac = 50.0452358d  ; center of display (WCS)
  decc = 8.7488722d
  boxsize = 16. ; arcsec
  
  openimplot, plotname, xsize=xsize, ysize=ysize
  
  !p.multi = [0,20,3]
   
  !p.position = [0.05, 0.1, 0.95, 0.76+0.10]
   
  t0 = floor(min(mjd))
  t = mjd - t0   
  xrange = [0.63, 0.7]

  t0str = strtrim(long(t0),2)
  t0str = strmid(t0str,0,2)+','+strmid(t0str,2,3)
   
  plot, [0], [0], xrange=xrange, yrange=[-5,36], xsty=5, ysty=9, charsize=cs*2, $
     ytitle='Flux density (!4m!17Jy)', xtitle='MJD-'+t0str
   
  axis, xaxis=0, xrange=xrange*24.*60, /xsty, charsize=cs*2, $ ; xtitle='Minutes after MJD '+t0str,
      xminor=10, xtickv=900+10*findgen(12), xticks=12
  xyouts, 0.424, 0.02, /norm, 'Minutes after MJD '+t0str, charsize=cs
  axis, xaxis=1, xrange=xrange*24*60., /xsty, xtickformat='(A1)', xtitle='', charsize=cs*2, $
      xminor=10, xtickv=900+10*findgen(12), xticks=12
 
  addmagaxis, xrange[1], 0.2, charsize=cs, minmag=21
  xyouts, xrange[1]+0.0005, abmag2ujy(20)-0.8, clip(20), charsize=cs
  oplot, t, f, col=200
  megaplot, t, f, psym=0, yerr=ef, errcol=251, $
       size=0.5, hatlen=0.25,  /noerase, color=251
  megaplot, t[seq], f[seq], size=1.45, /noerase, /fill, color=251
  for i = 0, n_elements(seq)-1 do begin
     xyouts, t[seq[i]]-0.00055, f[seq[i]]-0.48, '!18'+rclip(i+1,2)+'!17', color=255, charsize=0.57
  endfor
  oplot, xrange, [0,0], color=6
  xyouts, xrange[1]+0.0026, 21, 'AB magnitude', orient=-90, charsize=cs

  xyouts, xrange[0]+0.002, 32., 'TNT/ULTRASPEC (r-band)', charsize=cs

  mx = 0.05
  pdx = 0.08-0.04
  pdy = pdx * xsize/ysize
  sdx = 0.09-0.045
  sdy = sdx * xsize/ysize

  allimages = repstr(filename,'.subc.','.')  ; original images  
  images = allimages[seq]
  subimages = repstr(images,'.fits','.subc.fits')

  nrow = 20
  ncol = 3

  for i = 0, n_elements(subimages)-1 do begin
      if i gt 20 then break
      px = i mod nrow
      py = i / nrow
      !p.position = [mx+sdx*px,    0.96-py*sdy-pdy,$
                     mx+sdx*px+pdx,0.96-py*sdy]

      im = readfits(subimages[i], h, /silent)
      sky = median(im)
      addimage, im, h, rac, decc, boxsize, min=sky-60, max=sky+400, /invert

      xyouts, 0.3, 1.07, '!18'+rclip(i+1,2)+'!17', color=251, charsize=cs
      
   endfor
   
   xyouts, mx-0.01, 0.97-sdy*0.7-0.01, /norm, 'sub', orient=90, charsize=cs

  psclose
end


