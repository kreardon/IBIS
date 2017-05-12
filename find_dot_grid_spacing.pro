FUNCTION find_dot_grid_spacing, grid_image_input, start_pos_dot=start_pos_dot, step_size=step_size, $
             bootstrap=bootstrap, sobel_cutoff=sobel_cutoff, verbose=verbose, data_mask=data_mask, $
             dot_pos_map=dot_pos_map, rotation_grid=rotation_grid,$
            fix_radius=fix_radius,radius_guess=radius_guess,arcsec_step=arcsec_step

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
;        start_pos_dot    = pixel coordinates of dot in lower left corner; if it is not a two-element
;                               array, user will be prompted to click on lower-left dot.
;        step_size        = initial estimates of the pixel separation between dots
; OUTPUTS:
;        spatial_scale = calculated spatial scale in 
; KEYWORDS (Input):
;        bootstrap     = instructs the user to click on several dots on the image in order to
;                            to make an initial guess at the dot spacing
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
; KEYWORDS (Output):
;        dot_pos_map   = the map of [x,y,radius] parameters for all the fitted dots
;        rotation_grid = the calculated rotation in x and y of the dot grid

; COMMON BLOCKS:
;        None.
; SIDE EFFECTS:
; NOTES: 
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;        KPR
;-


IF N_ELEMENTS(bootstrap) NE 1        THEN bootstrap = 0
IF N_ELEMENTS(verbose) LT 1          THEN verbose = 0
IF N_ELEMENTS(sobel_hist_limit) LT 1 THEN  sobel_hist_limit = 0.9


; set to 1 to prompt user to click on starting dot, or set to zero to use a predefined starting point.
IF N_ELEMENTS(start_pos_dot) NE 2 THEN BEGIN
    click_point = 1
ENDIF ELSE BEGIN
    click_point = 0
ENDELSE

fit_type = 'fit_circle'
;fit_type = 'mpfitellipse'

IF N_ELEMENTS(radius_guess) LT 1 THEN radius_guess = 5.15
IF N_ELEMENTS(fix_radius)   LT 1 THEN fix_radius   = 1

dst_prime_focus_scale = 3.76  ; arcsec / mm
IF N_ELEMENTS(arcsec_step)  LT 1 THEN  arcsec_step = 0.5 * dst_prime_focus_scale


grid_image_size = SIZE(grid_image_input)

step_size_default = [19.6, 19.6]

grid_image = grid_image_input / MEDIAN(grid_image_input)

grid_im_sobel = SOBEL(grid_image)
;grid_im_sobel = madmax(grid_im)
grid_im_sobel_x_ave = REBIN(grid_im_sobel,grid_image_size[1],1)
pp = poly_fit(findgen(grid_image_size[1]),grid_im_sobel_x_ave,1,yfit=yfit)
grid_im_sobel_x_cor = REBIN(yfit,grid_image_size[1],grid_image_size[2])
grid_im_sobel_x_cor /= MEAN(grid_im_sobel_x_cor)
grid_im_sobel = grid_im_sobel / grid_im_sobel_x_cor
TVSCL,grid_im_sobel

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

; mask off edges where there may be extraneous structures
IF Keyword_Set(data_mask) THEN BEGIN
    grid_im_sobel_mask = grid_im_sobel_mask * data_mask
ENDIF

selected_points = TOTAL(grid_im_sobel_mask)/FLOAT(N_ELEMENTS(grid_im_sobel_mask))
IF verbose GE 1 THEN PRINT,'Found ' + STRING(selected_points * 100,FORMAT='(F5.1)') + '% of points above sobel cutoff'

IF click_point THEN BEGIN
    PRINT,'Click on center of grid dot in lower left corner of image'
    crs, clickx, clicky, /Device, /Quiet
    start_pos = [30,30]
    IF (clickx GE 10) AND (clickx LE 80) THEN start_pos[0] = clickx
    IF (clicky GE 10) AND (clicky LE 80) THEN start_pos[1] = clicky
    IF verbose GE 1 THEN PRINT, 'User selected starting position [ ' + STRTRIM( start_pos[0],2) + ', ' + STRTRIM( start_pos[1],2) + ' ]'
ENDIF ELSE BEGIN
    start_pos = start_pos_dot
    IF verbose GE 1 THEN PRINT, 'Input starting position [ ' + STRTRIM( start_pos[0],2) + ', ' + STRTRIM( start_pos[1],2) + ' ]'
ENDELSE

IF Keyword_Set(bootstrap) THEN BEGIN
    stride = 19.5
    TVLCT,rrct,ggct,bbct,/GET
    LOADCT,5
    
    tvcirc,start_pos[0] + step_size_default[0] * 2, start_pos[1] + step_size_default[1] * 2, 18,col=100
    PRINT,'Choosing grid dot two steps away in x- and y- directions;'
    PRINT,'Dot should be approximately in center of circle. Click on center of appropriate dot'
    crs, clickx_2, clicky_2, /Device, /Quiet
    stride_x = (clickx_2 - start_pos[0]) / 2.
    stride_y = (clicky_2 - start_pos[1]) / 2.
    
    FOR stp=1,10 do plots,start_pos[0] + stride_x * stp, start_pos[1] + stride_y * stp,psym=1,col=50,th=2,/Device

    tvcirc,start_pos[0] + stride_x * 10, start_pos[1] + stride_y * 10, 14,col=100
    PRINT,'Choosing grid dot ten steps away in x- and y- directions;'
    PRINT,'Appropriate dot should be approximately in center of circle. Click on center of dot'
    crs, clickx_3, clicky_3, /Device, /Quiet
    stride_x = (clickx_3 - start_pos[0]) / 10.
    stride_y = (clicky_3 - start_pos[1]) / 10.

    FOR stp=11,30 do plots,start_pos[0] + stride_x * stp, start_pos[1] + stride_y * stp,psym=1,col=50,th=2,/Device

    tvcirc,start_pos[0] + stride_x * 30, start_pos[1] + stride_y * 30, 8,col=100
    PRINT,'Choosing grid dot twenty steps away in x- and y- directions;'
    PRINT,'Appropriate dot should be approximately in center of circle. Click on center of dot'
    crs, clickx_4, clicky_4, /Device, /Quiet
    stride_x = (clickx_4 - start_pos[0]) / 30.
    stride_y = (clicky_4 - start_pos[1]) / 30.

    step_size = [stride_x, stride_y]
    TVLCT,rrct,ggct,bbct
    IF verbose GE 1 THEN PRINT,'Using bootstrapped [x,y] dot spacing of ' + STRING(step_size,FORMAT='("[", F7.4, ", ", F7.4, "].")')

ENDIF ELSE BEGIN
    IF N_ELEMENTS(step_size) LT 2 THEN step_size = step_size_default
    IF verbose GE 1 THEN PRINT,'Using input [x,y] dot spacing of ' + STRING(step_size,FORMAT='("[", F7.4, ", ", F7.4, "].")')
ENDELSE

;box_size  = [21,21]
box_size  = [17,17]
im_size   = [grid_image_size[1],grid_image_size[2]]
num_steps = FIX( (im_size - box_size) / step_size)

dot_pos = FLTARR(num_steps[0],num_steps[1],3)

trend_x = FLTARR( num_steps[1] )
trend_y = FLTARR( num_steps[0] )
radius_ave = radius_guess

FOR reps = 0,6 DO BEGIN
    TVSCL,grid_im_sobel_mask

    FOR xxp = 0,num_steps[0]-1 DO BEGIN
        FOR yyp = 0,num_steps[1]-1 DO BEGIN
            posx     = (start_pos[0] + xxp * step_size[0])>0<(im_size[0]-1)
            posy     = (start_pos[1] + yyp * step_size[1])>0<(im_size[1]-1)
            posx     = posx + trend_x[yyp]
            posy     = posy + trend_y[xxp]
            plots,posx,posy,psym=1,col=20,/dev
            dot_cutx = ([-10,10] + posx) >0<(im_size[0]-1)
            dot_cuty = ([-10,10] + posy) >0<(im_size[1]-1)
            dot_im   = grid_im_sobel_mask[dot_cutx[0]:dot_cutx[1],dot_cuty[0]:dot_cuty[1]]
            ;tvscl,SCALE(dot_im,5,5)
            dot_im_sz = [ N_ELEMENTS(dot_im[*,0]), N_ELEMENTS(dot_im[*,1]) ]
            dot_im_edge = WHERE(dot_im)
            dot_im_edge_x = dot_im_edge MOD dot_im_sz[0]
            dot_im_edge_y = FIX(dot_im_edge / dot_im_sz[0])
            
            ;tvcirc,dot_im_sz[0]/2.+dot_cutx[0],dot_im_sz[1]/2+dot_cuty[0],5.176,col=5e4

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
            dot_pos(xxp, yyp, *) = circle_coord
            tvcirc,circle_coord[0],circle_coord[1],circle_coord[2],col=200
        ENDFOR
    ENDFOR
    
    radius_ave = MEDIAN(dot_pos[*,*,2])

    trend_x  = REFORM(REBIN(dot_pos(*,*,0),1,num_steps[1],1))
    trend_x -= MEAN(trend_x)
    trend_y  = REFORM(REBIN(dot_pos(*,*,1),num_steps[0],1,1))
    trend_y -= MEAN(trend_y)

ENDFOR

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


RETURN, [spatial_scale_x, spatial_scale_y]

END


