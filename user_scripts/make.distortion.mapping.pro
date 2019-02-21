; This program computes the distortion maps that are needed to remap the IBIS images onto a 
; regular grid, with a desired spatial scale.
; The program starts from the map of the positions of the dots from the dot grid calibration, 
; interpolates those onto a grid with integer spacing and the desired dot separation (= spatial scale)

load_data       = 1
run_grid_check  = 1
do_save_results = 1

basedir = '/Users/kreardon/IBIS/ALMA-2016/Reardon/23Apr2017/'
date_str = '23Apr2017'
;channel  = 'whitelight'
channel  = 'spectral'
nb_wavelength = '6173'

; it might be preferred to specify the number of grid points to optimize the coverage (or avoid extra points at the edges).
num_steps = [49,49]
;num_steps = -1

IF load_data THEN BEGIN
    CASE channel OF
    'whitelight' : BEGIN
        RESTORE,verbose=0,basedir + 'DarkCalibration.whitelight.combineall.20170423_211809.series.ave.sav'
        wl_dark_ave = Series_Ave
        restore,verbose=0,'GridImages.whitelight.combineall.20170423_202240.series.ave.sav'
        wl_grid_info = Images_info
        wl_grid_ave = Series_Ave
        
        wl_cal_params = load_alma_calibration_info(date_str, 'ibis_wl')
        ibis_grid_use       = wl_grid_ave
        ibis_grid_use      -= wl_dark_ave
        ibis_grid_use       = ROTATE(ibis_grid_use,wl_cal_params.transpose)

        rotation_value = wl_cal_params.rot_to_grid
        step_size_output = [19.8987,19.4152]

        label = 'wl'
   
      END
    'spectral' : BEGIN
        RESTORE,verbose=0,basedir + 'DarkCalibration.spectral.combineall.20170423_211809.series.ave.sav'
        nb_dark_ave = Series_Ave
        restore,verbose=0,basedir + 'GridImages.spectral.byfilter.20170423_202240.series.ave.sav'
        nb_grid_info = Images_info
        nb_grid_input = Series_Ave_Input
        nb_grid_ave = Series_Ave

        nb_cal_params = load_alma_calibration_info(date_str, 'ibis_nb')
        filter_idx_select = WHERE(nb_cal_params.filter_ids EQ nb_wavelength)
        series_idx_select = WHERE(FIX(nb_grid_input[2,*,0])  EQ FIX(filter_idx_select[0]))
        PRINT,filter_idx_select,series_idx_select
        ibis_grid_use       = nb_grid_ave[*,*,series_idx_select]
        ibis_grid_use      -= nb_dark_ave
        ibis_grid_use       = ROTATE(ibis_grid_use,nb_cal_params.transpose)
        
        rotation_value = nb_cal_params.rot_to_grid
        step_size_output = [19.71,19.28]
        prefix = 'nb'
      END
    ENDCASE

    grid_im_rot =  ROT(SHIFT(ibis_grid_use,0,0) ,rotation_value,CUBIC=-0.5)
ENDIF

IF run_grid_check THEN BEGIN
    start_pos_input  = -1
    hh = find_dot_grid_spacing(grid_im_rot,start_pos_dot=start_pos_input,step_size=step_size_output,verbose=1,bootstrap=0,dot_pos_map=dot_pos_map,num_steps=num_steps)
ENDIF
$ls -l 
dot_pos_map_sz = SIZE(dot_pos_map)
np_half        = FIX(dot_pos_map_sz[1:2] / 2.)

; reorganize the maps of dot positions to match the scheme needed by the "doreg" destretching algorithm 
dot_pos_disp = FLTARR(2,dot_pos_map_sz[1],dot_pos_map_sz[2])
dot_pos_disp[0,*,*] = dot_pos_map[*,*,0]
dot_pos_disp[1,*,*] = dot_pos_map[*,*,1]

; Calculate expected locations of dots given a uniform plate scale with a specified spacing 
dot_pos_rdisp = dot_pos_disp
;step_size_even = MEAN(step_size_output)
step_size_even = 19.598    & wl_grid_scl        = ['0.096','arcsec/pixel']
;step_size_even = 18.865   & wl_grid_scl        = ['0.10','arcsec/pixel']
dot_pos_rdisp[0,*,*]  = REBIN((FINDGEN(dot_pos_map_sz[1],1) - np_half[0]) * step_size_even ,dot_pos_map_sz[1], dot_pos_map_sz[2])
dot_pos_rdisp[0,*,*] += dot_pos_disp[0,np_half[0],np_half[1]]    
dot_pos_rdisp[1,*,*]  = REBIN((FINDGEN(1,dot_pos_map_sz[2]) - np_half[0]) * step_size_even ,dot_pos_map_sz[1], dot_pos_map_sz[2]) 
dot_pos_rdisp[1,*,*] += dot_pos_disp[1,np_half[0],np_half[1]]    

; The "doreg" algortithm seems to introduce artifacts if the array of the control points
; does not have an integer spacing.
; So we calculate a grid of points with a integer spacing of grid points
step_size_integer = FIX(step_size_output) + 1

xran_min      = MIN(dot_pos_rdisp[0,*,np_half[1]])
xran_max      = MAX(dot_pos_rdisp[0,*,np_half[1]])
xran_interval = FIX((xran_max - xran_min) / step_size_integer[0]) + 1
xran_steps    = xran_interval + 1
xstart        = MEAN(dot_pos_disp(0,0,*))

yran_min      = MIN(dot_pos_rdisp[1,np_half[0],*])
yran_max      = MAX(dot_pos_rdisp[1,np_half[0],*])
yran_interval = FIX((yran_max - yran_min) / step_size_integer[1]) + 1
yran_steps    = yran_interval + 1
ystart        = MEAN(dot_pos_disp(1,*,0))

dot_pos_rdisp_even = FLTARR(2, xran_steps, yran_steps)
dot_pos_disp_even  = FLTARR(2, xran_steps, yran_steps)

dot_pos_rdisp_even[0,*,*] = REBIN(FINDGEN(xran_steps, 1) * step_size_integer[0], xran_steps, yran_steps) + ROUND(xstart)
dot_pos_rdisp_even[1,*,*] = REBIN(FINDGEN(1, yran_steps) * step_size_integer[1], xran_steps, yran_steps) + ROUND(ystart)

; Now we interpolate the dot coordinates determined from "find_dot_grid_spacing" onto the new, integer-spaced
; grid of control points. 
FOR nn=0,xran_steps-1 DO BEGIN
    dot_pos_disp_even(1,nn,*) = INTERPOL(dot_pos_disp(1,nn,*), dot_pos_rdisp(1,nn,*), dot_pos_rdisp_even(1,nn,*),/QUAD)
ENDFOR

FOR nn=0,yran_steps-1 DO BEGIN
    dot_pos_disp_even(0,*,nn) = INTERPOL(dot_pos_disp(0,*,nn), dot_pos_rdisp(0,*,nn), dot_pos_rdisp_even(0,*,nn),/QUAD)
ENDFOR

IF do_save_results THEN BEGIN
    ; save the control point (rdisp) and distortion (disp) maps
    variables = prefix + '_grid_even_rdisp, ' + prefix + '_grid_even_disp, ' + prefix + '_grid_scl, ' + prefix + '_grid_rot'
    res = EXECUTE(prefix + '_grid_even_rdisp = dot_pos_rdisp_even')
    res = EXECUTE(prefix + '_grid_even_disp  = dot_pos_disp_even')
    res = EXECUTE(prefix + '_grid_scl        = wl_grid_scl')
    res = EXECUTE(prefix + '_grid_rot        = -0.62')
    output_filename = prefix+'.' + nb_wavelength + '.' + date_str + '.destr.vect.to.even.scale.sav'
    res = EXECUTE('SAVE,Filename=output_filename,' + variables)
ENDIF

; remap the original dot grid image using the new distorion map
IF run_grid_check  THEN BEGIN
    grid_im_even_reg = doreg(grid_im_rot, dot_pos_rdisp_even, dot_pos_disp_even) 
    hh2 = find_dot_grid_spacing(grid_im_even_reg,start_pos_dot=-1,$
              step_size=[step_size_even,step_size_even],verbose=1,bootstrap=0,dot_pos_map=dot_pos_map,num_steps=num_steps)

ENDIF

END

;    spatscale_5896 = find_dot_grid_spacing(grid_im_even_reg_5896,start_pos_dot=start_pos_5896,$
;        step_size=[step_size_even,step_size_even],verbose=0,bootstrap=0,num_steps=num_steps)
;    spatscale_6173 = find_dot_grid_spacing(grid_im_even_reg_6173,start_pos_dot=start_pos_6173,$
;        step_size=[step_size_even,step_size_even],verbose=0,bootstrap=0,num_steps=num_steps)
;    spatscale_6563 = find_dot_grid_spacing(grid_im_even_reg_6563,start_pos_dot=start_pos_6563,$
;        step_size=[step_size_even,step_size_even],verbose=0,bootstrap=0,num_steps=num_steps)
;    spatscale_8542 = find_dot_grid_spacing(grid_im_even_reg_8542,start_pos_dot=start_pos_8542,$
;        step_size=[step_size_even,step_size_even],verbose=0,bootstrap=0,num_steps=num_steps)
