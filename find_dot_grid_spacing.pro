FUNCTION find_dot_grid_spacing, grid_image_input, start_pos_dot=start_pos_dot, step_size=step_size, $
             bootstrap=bootstrap, sobel_cutoff=sobel_cutoff, verbose=verbose, data_mask=data_mask, $
             dot_pos_map=dot_pos_map, rotation_grid=rotation_grid,fft_spacing=fft_spacing,$
            fix_radius=fix_radius,radius_guess=radius_guess,arcsec_step=arcsec_step, $
            correlation_refine=correlation_refine, num_steps=num_steps

;+
; NAME:
;        find_dot_grid_spacing
; PURPOSE:
;        locate positions of dots from dot grid pattern and use those to determine
;        the spatial scale and rotation of the image
; CATEGORY:
; CALLING SEQUENCE:
;        image_scale_xy = find_dot_grid_spacing(grid_image, start_pos_dot=1st_dot_pos, step_size=step_size_input)
; INPUTS:
;        grid_image_input = image of dot grid. It should be rotated to within ~0.2 degrees of having
;                               the dots aligned with the array axes
;        start_pos_dot    = [optional] pixel coordinates of dot in lower left corner; if it is not a two-element
;                               array, user will be prompted to click on lower-left dot.; updated with
;                               dot position determined by grid fitting process
;        step_size        = [optional] initial estimates of the pixel separation between dots; updated with
;                               step size determined from fitting process.
; OUTPUTS:
;        spatial_scale = calculated spatial scale in 
; KEYWORDS (Input):
;        bootstrap     = instructs the user to click on several dots on the image in order to
;                            to make an initial guess at the dot spacing
;        fft_spacing   = use an FFT technique to make initial guess at grid spacing;
;                            overrides any user input value for step_size.
;        data_mask     = a mask, the same size as the input image, to apply to the input image in
;                            order to eliminate points (e.g. the edges of the field) that might confuse
;                            the fitting algorithm
;        fix_radius    = use a fixed radius for the circle fitting routine
;        radius_guess  = initial guess for radius of dots
;        arcsec_step   = value to use for grid dot separation in arcseconds 
;                            [default = 1.88" suitable for 0.5 mm grid at DST]
;        verbose       = print more detailed information during and at end of fitting process
;        sobel_cutoff  = the input threshold to be applied to the image to define the 
;                            circle points around each dot to which to apply the fit
;        correlation_refine = [default=yes] use cross-correlation of the individual dot cutouts with
;                                 with the average dot image in order to better refine the dot positions.
;        num_steps     = two-element array used to manually set the number of dots that are fitted across
;                            the image. Can be used if the automatic determination of the number of dots 
;                            is too low (dots at the edge of the field are not fit) or 
;                            too high (dots do not fully fill the field of view).
;                                 
; KEYWORDS (Output):
;        dot_pos_map   = the map of [x,y,radius] parameters for all the fitted dots
;        rotation_grid = the calculated rotation in x and y of the dot grid
;
; COMMON BLOCKS:
;        None.
; SIDE EFFECTS:
; NOTES: 
;       Originally required tvcirc.pro from JHU/APL IDL library. That has been replaced to use Tvcircle.pro
;       from the more common IDL AstroLib.
; RESTRICTIONS:
;       Currently doesn't deal well with rotations of dot grid (with respect to array access) that are greater
;           than approximately 0.5 degrees (which causes an offset of about 1/2 a dot spacing distance over a 
;           1000 pixel image (or 50 dots).
; PROCEDURE:
; MODIFICATION HISTORY:
;        KPR  - April 2017 - original implementation
;-


IF N_ELEMENTS(bootstrap) NE 1           THEN bootstrap = 0
IF N_ELEMENTS(verbose) LT 1             THEN verbose = 0
IF N_ELEMENTS(sobel_hist_limit) LT 1    THEN  sobel_hist_limit = 0.9
IF N_ELEMENTS(correlation_refine) LT 1  THEN correlation_refine = 1
IF N_ELEMENTS(fix_radius)   LT 1        THEN fix_radius   = 1

do_sobel_trend_correction = 0

; set to 1 to prompt user to click on starting dot, or set to zero to use a predefined starting point.
IF N_ELEMENTS(start_pos_dot) NE 2 THEN BEGIN
    click_point = 1
ENDIF ELSE BEGIN
    click_point = 0
ENDELSE

fit_type = 'fit_circle'
;fit_type = 'mpfitellipse'
num_iter  = 5

dst_prime_focus_scale = 3.76  ; arcsec / mm
IF N_ELEMENTS(arcsec_step)  LT 1 THEN  arcsec_step = 0.5 * dst_prime_focus_scale

; this ratio of the dot spacing to the radius of the individual dots 
; It is appropriate for the dot grid used at the Dunn Solar Telescope but may be different for other dot grids.
dot_radius_to_spacing_ratio = 0.25

grid_image_size = SIZE(grid_image_input)
im_size         = [grid_image_size[1],grid_image_size[2]] 

grid_image      = grid_image_input / MEDIAN(grid_image_input)

IF click_point THEN BEGIN
    TVSCL,grid_image
    PRINT,'Click on center of grid dot in lower left corner of image'
    TVCRS,50,50,/Device
    cursor,clickx, clicky, 3, /Device, /Wait

    start_pos = [30,30]
    IF (clickx GE 5) AND (clickx LE 100) THEN start_pos[0] = clickx
    IF (clicky GE 5) AND (clicky LE 100) THEN start_pos[1] = clicky
    IF verbose GE 1 THEN PRINT, 'User selected starting position [ ' + STRTRIM( start_pos[0],2) + ', ' + STRTRIM( start_pos[1],2) + ' ]'
ENDIF ELSE BEGIN
    start_pos = start_pos_dot
    IF verbose GE 1 THEN PRINT, 'Input starting position [ ' + STRTRIM( start_pos[0],2) + ', ' + STRTRIM( start_pos[1],2) + ' ]'
ENDELSE
Print,start_pos

grid_im_sobel = SOBEL(grid_image)
;grid_im_sobel = madmax(grid_im)

IF do_sobel_trend_correction THEN BEGIN
    ; try to remove horizontal trends in the Sobel values (perhaps due to spatially dependent
    ;     blurring or scattered light
    grid_im_sobel_x_ave = REBIN(grid_im_sobel,grid_image_size[1],1)
    pp = poly_fit(findgen(grid_image_size[1]),grid_im_sobel_x_ave,1,yfit=yfit)
    grid_im_sobel_x_cor = REBIN(yfit,grid_image_size[1],grid_image_size[2])
    grid_im_sobel_x_cor /= MEAN(grid_im_sobel_x_cor)
    grid_im_sobel = grid_im_sobel / grid_im_sobel_x_cor
    
    grid_im_sobel_y_ave = REBIN(grid_im_sobel,1,grid_image_size[2])
    pp = poly_fit(findgen(grid_image_size[2]),grid_im_sobel_y_ave,1,yfit=yfit)
    grid_im_sobel_y_cor = REBIN(yfit,grid_image_size[1],grid_image_size[2])
    grid_im_sobel_y_cor /= MEAN(grid_im_sobel_y_cor)
    grid_im_sobel = grid_im_sobel / grid_im_sobel_y_cor
ENDIF

TVSCL,grid_im_sobel

; automatically determine the optimal cutoff value for the Sobel values to find 
; edge of circles for fitting procedure
IF N_ELEMENTS(sobel_cutoff) NE 1 THEN BEGIN
    num_hist_bins             = 200
    grid_im_sobel_hist        = histogram(grid_im_sobel,nbins=num_hist_bins,loc=xscl)
    grid_im_sobel_hist_sum    = FLTARR(num_hist_bins)
    grid_im_sobel_hist_sum[0] = grid_im_sobel_hist[0]
    for nn=1,num_hist_bins-1 do grid_im_sobel_hist_sum[nn] = grid_im_sobel_hist_sum[nn-1] + grid_im_sobel_hist[nn]
    sobel_cutoff              = xscl(MIN(WHERE(grid_im_sobel_hist_sum GE MAX(grid_im_sobel_hist_sum)*sobel_hist_limit)))
    IF verbose GE 1 THEN PRINT,'Determined sobel cutoff - ',STRTRIM(sobel_cutoff,2)
ENDIF 

grid_im_sobel_mask = grid_im_sobel GE sobel_cutoff

; If a data mask was provided then apply mask to avoid masked areas 
;    (e.g. edges where dots my be cut off or too close to edge)
IF Keyword_Set(data_mask) THEN BEGIN
    grid_im_sobel_mask = grid_im_sobel_mask * data_mask
ENDIF

selected_points = TOTAL(grid_im_sobel_mask)/FLOAT(N_ELEMENTS(grid_im_sobel_mask))
IF verbose GE 1 THEN PRINT,'Found ' + STRING(selected_points * 100,FORMAT='(F5.1)') + '% of points above sobel cutoff'

IF N_ELEMENTS(step_size) LT 2 THEN step_size_default = [19.6, 19.6] ELSE step_size_default = step_size

IF Keyword_Set(fft_spacing) THEN BEGIN
    grid_im_fft_x = FLTARR(grid_image_size[1])
    apod_win      = apod(grid_image_size[1], 1, 0.04, 0.0, 2)
    FOR nn=grid_image_size[2]*0.1,grid_image_size[2]*0.9 DO BEGIN
       grid_im_fft_x += ABS(FFT((grid_image[*,nn] - MEDIAN(grid_image[*,nn])) * apod_win,-1))
    ENDFOR
    grid_im_fft_x_strt    = FIX(grid_image_size[1] * 0.005) + 1
    grid_im_fft_x_max     = MAX(grid_im_fft_x[grid_im_fft_x_strt:grid_image_size[1]/2.], grid_im_fft_x_maxpos)
    grid_im_fft_x_maxpos += grid_im_fft_x_strt
    grid_im_fft_x_lc      = lc_find(-grid_im_fft_x,grid_im_fft_x_maxpos-10,grid_im_fft_x_maxpos+10,3)
    grid_im_fft_xstep     = grid_image_size[1]/grid_im_fft_x_lc[0]

    grid_im_fft_y = FLTARR(grid_image_size[2])
    apod_win      = apod(grid_image_size[2], 1, 0.04, 0.0, 2)
    FOR nn=grid_image_size[1]*0.1,grid_image_size[1]*0.9 DO BEGIN
       grid_im_fft_y += ABS(FFT((grid_image[nn,*] - MEDIAN(grid_image[nn,*])) * apod_win,-1))
    ENDFOR
    grid_im_fft_y_strt    = FIX(grid_image_size[2] * 0.005) + 1
    grid_im_fft_y_max     = MAX(grid_im_fft_y[grid_im_fft_y_strt:grid_image_size[1]/2.], grid_im_fft_y_maxpos)
    grid_im_fft_y_maxpos += grid_im_fft_y_strt
    grid_im_fft_y_lc      = lc_find(-grid_im_fft_y,grid_im_fft_y_maxpos-10,grid_im_fft_y_maxpos+10,3)
    grid_im_fft_ystep     = grid_image_size[2]/grid_im_fft_y_lc[0]
    
    step_size_default = [grid_im_fft_xstep, grid_im_fft_ystep]
ENDIF

IF Keyword_Set(bootstrap) THEN BEGIN
    stride = 19.5
    TVLCT,rrct,ggct,bbct,/GET
    LOADCT,5
        
    ;tvcirc,start_pos[0] + step_size_default[0] * 2, start_pos[1] + step_size_default[1] * 2, 18,col=100
    Tvcircle, 18, start_pos[0] + step_size_default[0] * 2, start_pos[1] + step_size_default[1] * 2, 100, /Device
    TVCRS,start_pos[0] + step_size_default[0] * 2, start_pos[1] + step_size_default[1] * 2,/Device
    PRINT,'Click on grid dot two steps away in x- and y- directions;'
    PRINT,'Dot should be approximately in center of circle. Click on center of appropriate dot'
    cursor, clickx_2, clicky_2, /Down, /Device
    wait,0.3
 
    stride_x = (clickx_2 - start_pos[0]) / 2.
    stride_y = (clicky_2 - start_pos[1]) / 2.
        
    FOR stp=1,10 do plots,start_pos[0] + stride_x * stp, start_pos[1] + stride_y * stp,psym=1,col=50,th=2,SymSize=2,/Device
    
    ;tvcirc,start_pos[0] + stride_x * 10, start_pos[1] + stride_y * 10, 14,col=100
    Tvcircle, 14, start_pos[0] + stride_x * 10, start_pos[1] + stride_y * 10, 100, /Device
    Tvcrs, start_pos[0] + stride_x * 10, start_pos[1] + stride_y * 10, /Device
    PRINT,'Click on grid dot ten steps away in x- and y- directions;'
    PRINT,'Appropriate dot should be approximately in center of circle. Click on center of dot'
    cursor,clickx_3, clicky_3, /Down, /Device
    wait,0.3

    stride_x = (clickx_3 - start_pos[0]) / 10.
    stride_y = (clicky_3 - start_pos[1]) / 10.
    
    FOR stp=11,30 do plots,start_pos[0] + stride_x * stp, start_pos[1] + stride_y * stp,psym=1,col=50,th=2,/Device
    
    ;tvcirc,start_pos[0] + stride_x * 30, start_pos[1] + stride_y * 30, 8,col=100
    Tvcircle, 8, start_pos[0] + stride_x * 30, start_pos[1] + stride_y * 30, 100, /Device
    TVCRS, start_pos[0] + stride_x * 30, start_pos[1] + stride_y * 30, /Device
    PRINT,'Click on grid dot twenty steps away in x- and y- directions;'
    PRINT,'Appropriate dot should be approximately in center of circle. Click on center of dot'
    cursor,clickx_4, clicky_4, /Down, /Device
    wait,0.3

    stride_x = (clickx_4 - start_pos[0]) / 30.
    stride_y = (clicky_4 - start_pos[1]) / 30.
    
    step_size = [stride_x, stride_y]
    TVLCT,rrct,ggct,bbct
    dot_spacing_type = 'bootstrapped'
ENDIF ELSE IF Keyword_Set(fft_spacing) THEN BEGIN
    step_size = step_size_default
    dot_spacing_type = 'FFT-determined'
ENDIF ELSE IF N_ELEMENTS(step_size) EQ 2 THEN BEGIN
    step_size = step_size
    dot_spacing_type = 'user-input'
ENDIF ELSE BEGIN 
    step_size = step_size_default
    dot_spacing_type = 'default'
ENDELSE

IF verbose GE 1 THEN PRINT,'Using ' + dot_spacing_type + ' [x,y] dot spacing of ' + STRING(step_size,FORMAT='("[", F7.4, ", ", F7.4, "].")')

IF N_ELEMENTS(radius_guess) LT 1        THEN radius_guess = MEAN(step_size) * dot_radius_to_spacing_ratio

; find the half size of the box size used to fit around the dot, the largest even number smaller than half of the dot spacing
box_half   = FIX(FIX(MEAN(step_size)/2.) / 2) * 2
IF N_ELEMENTS(num_steps) LE 1 THEN BEGIN
    ;box_size  = [21,21]
    box_size  = step_size
    num_steps = FIX( (im_size - start_pos - box_size) / step_size) + 1
ENDIF ELSE BEGIN
    num_steps = num_steps[0:1]
ENDELSE
IF verbose GE 1 THEN PRINT,'Size of dot grid: x = ' + STRTRIM(num_steps[0],2) + ' ; y = ' + STRTRIM(num_steps[1],2)

dot_pos = FLTARR(num_steps[0],num_steps[1],3)

trend_x = FLTARR( num_steps[1] )
trend_y = FLTARR( num_steps[0] )
radius_ave = radius_guess

FOR reps = 0,num_iter - 1 DO BEGIN
    TVSCL,grid_im_sobel_mask

    FOR xxp = 0,num_steps[0]-1 DO BEGIN
        FOR yyp = 0,num_steps[1]-1 DO BEGIN
            posx     = (start_pos[0] + xxp * step_size[0])>0<(im_size[0]-1)
            posy     = (start_pos[1] + yyp * step_size[1])>0<(im_size[1]-1)
            posx     = posx + trend_x[yyp]
            posy     = posy + trend_y[xxp]
            plots,posx,posy,psym=1,col=20,/dev
            dot_cutx = ([-box_half,box_half] + posx) >0<(im_size[0]-1)
            dot_cuty = ([-box_half,box_half] + posy) >0<(im_size[1]-1)
            dot_im   = grid_im_sobel_mask[dot_cutx[0]:dot_cutx[1],dot_cuty[0]:dot_cuty[1]]
            ;tvscl,SCALE(dot_im,5,5)
            dot_im_sz = [ N_ELEMENTS(dot_im[*,0]), N_ELEMENTS(dot_im[*,1]) ]
            dot_im_edge = WHERE(dot_im, dot_edge_count)
            dot_im_edge_x = dot_im_edge MOD dot_im_sz[0]
            dot_im_edge_y = FIX(dot_im_edge / dot_im_sz[0])
            
            ;tvcirc,dot_im_sz[0]/2.+dot_cutx[0],dot_im_sz[1]/2+dot_cuty[0],5.176,col=5e4
            
            IF dot_edge_count GE 5 THEN BEGIN
                CASE fit_type OF
                    'fit_circle' : BEGIN
                        circle_guess = [dot_im_sz[0]/2.,dot_im_sz[1]/2.,radius_ave]
                        circle_coord = fit_circle(dot_im_edge_x, dot_im_edge_y,circle_guess,radius_fix=fix_radius)
                        circle_coord[0:1] += [dot_cutx[0],dot_cuty[0]]
                    END
                    'mpfitellipse' : BEGIN
                        circle_guess = [radius_ave,radius_ave,dot_im_sz[0]/2.,dot_im_sz[1]/2.,0]
                        mpfit_params = mpfitellipse(dot_im_edge_x, dot_im_edge_y,circle_guess,/Circular,/Quiet)
                        circle_coord = [mpfit_params[2]>(-10)<10, mpfit_params[3]>(-10)<10, MEAN(mpfit_params[0:1]) ]
                        circle_coord[0:1] += [dot_cutx[0],dot_cuty[0]]
                    END
                ENDCASE
            ENDIF ELSE BEGIN
                circle_coord = [posx, posy, radius_ave]
            ENDELSE
            
            dot_pos(xxp, yyp, *) = circle_coord
            ;tvcirc,circle_coord[0],circle_coord[1],circle_coord[2],col=200
            Tvcircle, circle_coord[2], circle_coord[0], circle_coord[1], 200, /Device

        ENDFOR
    ENDFOR
    
    radius_ave = MEDIAN(dot_pos[*,*,2])

    trend_x  = REFORM(REBIN(dot_pos(*,*,0),1,num_steps[1],1))
    trend_x -= MEAN(trend_x)
    trend_y  = REFORM(REBIN(dot_pos(*,*,1),num_steps[0],1,1))
    trend_y -= MEAN(trend_y)

ENDFOR

IF correlation_refine THEN BEGIN

    dot_box_size       = ROUND(step_size * 1.6)
    dot_box_size      += (dot_box_size + 1) MOD 2
    box_half           = FIX(dot_box_size/2.)
    dot_im_ave         = FLTARR(dot_box_size[0],dot_box_size[1])  
    dot_box_corsz      = ROUND(dot_box_size * 0.8 )
    dot_pos_ave_ccor   = FLTARR(num_steps[0],num_steps[1],2)
    
    for xx=1,num_steps[0]-2 do begin
        for yy=1,num_steps[1]-2 do begin
            pos_guess = ROUND(REFORM(dot_pos[xx,yy,0:1]))
            dot_im_ave += grid_image[pos_guess[0]-box_half[0]:pos_guess[0]+box_half[0],$
                                        pos_guess[1]-box_half[1]:pos_guess[1]+box_half[1]]
      endfor
    endfor
    
    for xx=0,num_steps[0]-1 do begin
        for yy=0,num_steps[1]-1 do begin
            pos_guess = ROUND(REFORM(dot_pos[xx,yy,0:1]))
            pos_guess[0] = pos_guess[0]>box_half[0]<(im_size[0]-1-box_half[0])
            pos_guess[1] = pos_guess[1]>box_half[1]<(im_size[1]-1-box_half[1])
            dotim = grid_image[pos_guess[0]-box_half[0]:pos_guess[0]+box_half[0],pos_guess[1]-box_half[1]:pos_guess[1]+box_half[1]]
            dot_shift = xyoff(dot_im_ave,dotim,dot_box_corsz[0],dot_box_corsz[1],/Quiet)
            dot_pos_ave_ccor[xx,yy,*] = pos_guess - dot_shift
        endfor
    endfor

    dot_pos[*,*,0:1] = dot_pos_ave_ccor

ENDIF

step_x_diff_xaxis = (dot_pos(1:num_steps[0]-1,*,0) - dot_pos(0:num_steps[0]-2,*,0))
step_x_xaxis_ave  = MEAN( step_x_diff_xaxis )
step_y_diff_xaxis = (dot_pos(1:num_steps[0]-1,*,1) - dot_pos(0:num_steps[0]-2,*,1))
step_y_xaxis_ave  = MEAN( step_y_diff_xaxis )
spatial_scale_x   = arcsec_step/SQRT(step_x_xaxis_ave^2 + step_y_xaxis_ave^2)

step_x_diff_yaxis = (dot_pos(*,1:num_steps[1]-1,0) - dot_pos(*,0:num_steps[1]-2,0))
step_x_yaxis_ave  = MEAN( step_x_diff_yaxis )
step_y_diff_yaxis = (dot_pos(*,1:num_steps[1]-1,1) - dot_pos(*,0:num_steps[1]-2,1))
step_y_yaxis_ave  = MEAN( step_y_diff_yaxis )
spatial_scale_y   = arcsec_step/SQRT(step_x_yaxis_ave^2 + step_y_yaxis_ave^2)

trend_y_fit       = LinFit(rebin(dot_pos(*,*,0),num_steps[0],1,1),trend_y)
trend_x_fit       = LinFit(rebin(dot_pos(*,*,1),1,num_steps[1],1),trend_x)
yangle            = - ATAN(trend_x_fit[1]) * !RADEG
xangle            =   ATAN(trend_y_fit[1]) * !RADEG
rotation_grid     = [xangle, xangle]

step_size         = [step_x_xaxis_ave, step_y_yaxis_ave]
dot_pos_map       = dot_pos
start_pos_dot     = REFORM(dot_pos[0,0,0:1])

IF verbose GE 1 THEN BEGIN
    PRINT, 'Average Dot Spacing:'
    PRINT, step_x_xaxis_ave, FORMAT='("    x-axis: ", F7.3, " pixel / dot")'
    PRINT, step_y_yaxis_ave, FORMAT='("    y-axis: ", F7.3, " pixel / dot")'
    
    PRINT, 'Calculated Plate Scales:'
    PRINT, spatial_scale_x, FORMAT='("    x-axis: ", F7.5, " arcsec/pixel")'
    PRINT, spatial_scale_y, FORMAT='("    y-axis: ", F7.5, " arcsec/pixel")'
     
    PRINT, yangle, xangle,FORMAT='( "Estimated Grid Rotation - ", F6.2, ", ", F6.2, " deg")'
ENDIF

TVSCL,grid_image
FOR xx=0,num_steps[0]-1 DO BEGIN
    FOR yy=0,num_steps[1]-1 DO BEGIN
        ;tvcirc,dot_pos[xx,yy,0],dot_pos[xx,yy,1],dot_pos[xx,yy,2],col=200
        Tvcircle, dot_pos[xx,yy,2], dot_pos[xx,yy,0], dot_pos[xx,yy,1], 'red3', /Device
    ENDFOR
ENDFOR

RETURN, [spatial_scale_x, spatial_scale_y]

END


