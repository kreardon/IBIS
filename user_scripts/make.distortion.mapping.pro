; This program computes the distortion maps that are needed to remap the IBIS images onto a 
; regular grid, with a desired spatial scale.
; The program starts from the map of the positions of the dots from the dot grid calibration, 
; interpolates those onto a grid with integer spacing and the desired dot separation (= spatial scale)

load_data       = 1
run_grid_check  = 1
do_save_results = 1

;basedir       = '/SMdata1/kreardon/IBIS/VAULT/30Sep2014/ibis/'
basedir       = '/Users/kreardon/IBIS/VAULT/IBIS/'
ave_ser_dir   = basedir + 'averaged_series'
date_str      = '30Sep2014'
day_id        = '20140930'
;channel       = 'whitelight'
channel       = 'spectral'
nb_wavelength = '7090'
dot_grid_series = '20140930_201441'

; it might be preferred to specify the number of grid points to optimize the coverage (or avoid extra points at the edges).
num_steps_wl = [48,49]
num_steps_nb = [50,51]
;num_steps = -1

IF load_data THEN BEGIN
    CASE channel OF
    'whitelight' : BEGIN
        dark_wl_files      = file_search(ave_ser_dir,'DarkCalibration.whitelight.combineall.' + day_id + '*.series.ave.sav', count=num_darks)     
        wl_darks_all       = FLTARR(1000,1000,num_darks) 
        wl_darks_timerange = DBLARR(2,num_darks)
        FOR nn=0,num_darks-1 DO BEGIN
            RESTORE,Verbose=0,dark_wl_files[nn]
            wl_darks_all[*,*,nn]      = Series_Ave
            validp                    = (WHERE(strlen(Series_Ave_Input[5,*,*]) GE 1))
            wl_darks_timerange[*,nn]  = [MIN(fits_date_convert((SERIES_AVE_INPUT[5,*,*])[validp]),max=maxval),maxval]
        ENDFOR

        flat_wl_files      = file_search(ave_ser_dir,'FlatFieldCalibration.whitelight.combineall.' + day_id + '*.series.ave.sav', count=num_flats)     
        wl_flats_all       = FLTARR(1000,1000,num_flats) 
        wl_flats_timerange = DBLARR(2,num_flats)
        FOR nn=0,num_flats-1 DO BEGIN
            RESTORE,Verbose=0,flat_wl_files[nn]
            wl_flats_all[*,*,nn]      = Series_Ave
            validp                    = (WHERE(strlen(Series_Ave_Input[5,*,*]) GE 1))
            wl_flats_timerange[*,nn]  = [MIN(fits_date_convert((SERIES_AVE_INPUT[5,*,*])[validp]),max=maxval),maxval]
            dark_match                = get_closest(REBIN(wl_darks_timerange,1,num_darks),MEAN(wl_flats_timerange[*,nn]))
            wl_flats_all[*,nn]       -= wl_darks_all[*,*,dark_match]
            wl_flats_all[*,nn]       /= MEDIAN(wl_flats_all[*,nn])
        ENDFOR
        
        dot_grid_file = file_search(ave_ser_dir,'GridImages.whitelight.combineall*' + dot_grid_series + '*.series.ave.sav', count=num_darks) 
        restore,verbose=0,dot_grid_file[0]

        wl_grid_info = Images_info
        wl_grid_ave = Series_Ave

        validp                    = (WHERE(strlen(Series_Ave_Input[5,0,*]) GE 1))
        wl_grid_timerange         = [MIN(fits_date_convert((SERIES_AVE_INPUT[5,*,*])[validp]),max=maxval),maxval]
        dark_match                = get_closest(REBIN(wl_darks_timerange,1,num_darks),MEAN(wl_grid_timerange))
        flat_match                = get_closest(REBIN(wl_flats_timerange,1,num_flats),MEAN(wl_grid_timerange))
        
        wl_cal_params = load_calibration_info(date_str, 'ibis_wl')
        ibis_grid_use       = wl_grid_ave
        ibis_grid_use      -= wl_darks_all[*,*,dark_match]
        ibis_grid_use      /= wl_flats_all[*,*,flat_match]
        ibis_grid_use       = ROTATE(ibis_grid_use,wl_cal_params.transpose)

        rotation_value = wl_cal_params.rot_to_grid
        step_size_output = 1/wl_cal_params.plate_scale * 1.88 ; 1.88" / 0.5 mm = dot spacing
        ;step_size_output = [19.8987,19.4152]
        num_steps = num_steps_wl

        prefix = 'wl'
        output_filename = prefix+'.' +  date_str + '.destr.vect.to.even.scale.sav'
   
      end
    'spectral' : begin
        dark_nb_files      = file_search(ave_ser_dir,'DarkCalibration.spectral.combineall.' + day_id + '*.series.ave.sav', count=num_darks)     
        nb_darks_all       = fltarr(1000,1000,num_darks) 
        nb_darks_timerange = dblarr(2,num_darks)
        for nn=0,num_darks-1 do begin
            restore,verbose=0,dark_nb_files[nn]
            nb_darks_all[*,*,nn]      = series_ave
            validp                    = (where(strlen(series_ave_input[5,*,*]) ge 1))
            nb_darks_timerange[*,nn]  = [min(fits_date_convert((series_ave_input[5,*,*])[validp]),max=maxval),maxval]
        endfor

        flat_nb_files      = file_search(ave_ser_dir,'flatfieldcalibration.spectral.byfilter+wave.' + day_id + '*.series.ave.sav', count=num_flats)     
        
        dot_grid_file = file_search(ave_ser_dir,'GridImages.spectral.byfilter*' + dot_grid_series + '*.series.ave.sav', count=num_darks) 
        restore,verbose=0,dot_grid_file[0]

        nb_grid_info  = images_info
        nb_grid_input = series_ave_input
        nb_grid_ave   = series_ave

        nb_cal_params = load_calibration_info(date_str, 'ibis_nb')
        filter_idx_select = where(nb_cal_params.filter_ids eq nb_wavelength)
        series_idx_select = where(fix(nb_grid_input[2,*,0])  eq fix(filter_idx_select[0]))
        print,filter_idx_select,series_idx_select

        validp                    = (where(strlen(nb_grid_input[5,0,*]) ge 1))
        nb_grid_timerange         = [min(fits_date_convert((nb_grid_input[5,*,*])[validp]),max=maxval),maxval]

        dark_match                = get_closest(rebin(nb_darks_timerange,1,num_darks),mean(nb_grid_timerange))

        ibis_grid_use       = nb_grid_ave[*,*,series_idx_select]
        ibis_grid_use      -= nb_darks_all[*,*,dark_match]
        ibis_grid_use       = rotate(ibis_grid_use,nb_cal_params.transpose)
        
        rotation_value = nb_cal_params.rot_to_grid
        ;step_size_output = [19.70,19.30]
        
        step_size_output = 1/nb_cal_params.plate_scale[0:1,filter_idx_select] * 1.88 ; 1.88" / 0.5 mm = dot spacing
        num_steps = num_steps_nb

        prefix = 'nb'
        output_filename = prefix+'.' + nb_wavelength + '.' + date_str + '.destr.vect.to.even.scale.sav'

      end
    endcase

    grid_im_rot =  ROT(SHIFT(ibis_grid_use,0,0) ,rotation_value,CUBIC=-0.5)
ENDIF

IF run_grid_check THEN BEGIN
    start_pos_input  = -1
    hh = find_dot_grid_spacing(grid_im_rot,start_pos_dot=start_pos_input,step_size=step_size_output,verbose=1,bootstrap=0,dot_pos_map=dot_pos_map,num_steps=num_steps)
ENDIF

dot_pos_map_sz = SIZE(dot_pos_map)
np_half        = FIX(dot_pos_map_sz[1:2] / 2.)

; reorganize the maps of dot positions to match the scheme needed by the "doreg" destretching algorithm 
dot_pos_disp = FLTARR(2,dot_pos_map_sz[1],dot_pos_map_sz[2])
dot_pos_disp[0,*,*] = dot_pos_map[*,*,0]
dot_pos_disp[1,*,*] = dot_pos_map[*,*,1]

; Calculate expected locations of dots given a uniform plate scale with a specified spacing 
dot_pos_rdisp = dot_pos_disp
;step_size_even = MEAN(step_size_output)
step_size_even = 19.5833   & even_grid_scl        = ['0.096','arcsec/pixel']
;step_size_even = 19.598    & even_grid_scl        = ['0.096','arcsec/pixel']
;step_size_even = 18.865   & even_grid_scl        = ['0.10','arcsec/pixel']
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
    res = EXECUTE(prefix + '_grid_scl        = even_grid_scl')
    res = EXECUTE(prefix + '_grid_rot        = grid_im_rot')
    res = EXECUTE('SAVE,Filename=output_filename,' + variables)
ENDIF

; remap the original dot grid image using the new distorion map
IF run_grid_check  THEN BEGIN
    grid_im_even_reg = doreg(grid_im_rot, dot_pos_rdisp_even, dot_pos_disp_even) 
    hh2 = find_dot_grid_spacing(grid_im_even_reg,start_pos_dot=-1,$
              step_size=[step_size_even,step_size_even],verbose=1,bootstrap=1,dot_pos_map=dot_pos_map,num_steps=num_steps-[1,1])

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
