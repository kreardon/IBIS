FUNCTION load_calibration_info, date_str, instrument_channel, cal_directory=cal_directory

; sample usage:
; nb_21apr_info = load_alma_calibration_info('21Apr2017','ibis_nb') 

; these orientations should be pretty standard for IBIS 
; they shouldn't need to be changed for most observations.
ibis_wl_transpose    =  3
ibis_nb_transpose    =  7
; these values for ROSA may change depending on how the feed optics were setup
;rosa_gband_transpose =  6
;rosa_cak_transpose   =  1

; this is a "standard" value which was adopted but could be changed
target_plate_scale   = 0.096

; look for a "calibration_files" directory as a subdirectory of the directory
; where this program is found
IF NOT KEYWORD_SET(cal_directory) THEN BEGIN
    repository_location = File_Dirname(Routine_Filepath(/Either),/Mark)
    calibration_location = repository_location + 'calibration_files/'
ENDIF ELSE BEGIN
    calibration_location = cal_directory
ENDELSE

    ; name of IDL save file containing calibration frames and
    ; name of variable in save file for desired calibration type
    ; the user may save dark and gain files with whatever names desired and update information here
    wl_dark_cal_file      = 'whitelight.dark.ave.25Apr2019.sav'
    wl_dark_name          = 'wl_dark_ave'

    wl_gain_cal_file      = 'whitelight.gain.ave.25Apr2019.sav'
    wl_gain_name          = 'wl_gain'

    nb_dark_cal_file      = 'nb.dark.ave.25Apr2019.sav'
    nb_dark_name          = 'nb_dark_ave'
    
    nb_7090_gain_cal_file = 'Fe7090.gain.info.20190425.ave.sav'
    nb_7090_gain_name     = ['nb_gain_7090_ave', '7090', 'nb_gain_info_7090']

    nb_7699_gain_cal_file = 'K7699.gain.info.20190425.ave.sav'
    nb_7699_gain_name     = ['nb_gain_k7699_ave','7699', 'nb_gain_k7699_info_ave']  

    nb_5434_gain_cal_file  = 'Fe5434.gain.info.20190425.ave.sav'
    nb_5434_gain_name      = ['nb_gain_5434_ave','5434', 'nb_gain_5434_info']

    nb_gain_cal_file      = [[nb_7090_gain_cal_file], [nb_7699_gain_cal_file], [nb_5434_gain_cal_file]]
    nb_gain_name          = [[nb_7090_gain_name], [nb_7699_gain_name], [nb_5434_gain_name]]

    nb_bad_pixel_file      = 'narrowband.bad.pixel.index.txt'
    wl_bad_pixel_file      = 'whitelight.bad.pixel.index.txt'

    nb_7090_blue_map_file  = calibration_location + 'blueshift.map.ave.7090.25Apr2019.sav'
    nb_7090_blue_map_name  = ['blueshift_map_7090', '7090']
    nb_7699_blue_map_file  = calibration_location + 'blueshift.map.ave.7699.25Apr2019.sav'
    nb_7699_blue_map_name  = ['blueshift_map_7699', '7699']
    nb_5434_blue_map_file  = calibration_location + 'blueshift.map.ave.5434.25Apr2019.sav'
    nb_5434_blue_map_name  = ['blueshift_map_5434', '5434']
    nb_blue_map_file       = [[nb_7090_blue_map_file], [nb_7699_blue_map_file], [nb_5434_blue_map_file]]
    nb_blue_map_name       = [[nb_7090_blue_map_name], [nb_7699_blue_map_name], [nb_5434_blue_map_name]]
    
    ; pick a reference time for some time-dependent values
    ; typically pick the time of midnight before the start of the observations
    time_ref                = (fits_date_convert('2019-04-25T00:00:00'))[0]

    ; this is how a fixed value for rotation would be defined
    ;rot_wl_to_sol_north     = [0.32]    ; determined from fitting to HMI
    ; this is a quadratic fit to the observed rotation
    ; where the independent variable is fractional days since the reference time
    rot_wl_to_sol_north     = [-14.502663, 39.010681, -26.255281] ; determined from fitting to HMI

    ; reference time is midnight in UT
    rot_wl_to_sol_reftime   = time_ref

    ; IBIS narrowband plate scales - determined from dot grid using find_dot_grid_spacing.pro
    ;scale_ibis_wl          = [ 0.0975688, 0.0977281]  ; determined from dot grid

    ; IBIS narrowband plate scales - determined from fitting to HMI using stx_findlocation.pro
    ; this generally requires running alignment over many whitelight vs. HMI combinations to beat down the noise
    ; but in this case we found a plate scale that was about 2% larger than from the dot grid
    ; Don't know which is correct in an absolute sense, but the key is to match the whitelight image
    ; to the HMI plate scale, so we'll go with HMI-derived value
    scale_ibis_wl           = [ 0.0976, 0.0976 ]   ; mean value, determined from fitting to HMI, about 0.2% larger
    scale_ibis_wl           = [ 0.0975204, 0.0976796] ; asymmetric values, combining HMI average with ratio derived from dot grid

    ; the rotation angle of the whitelight dot grid determined from find_dot_grid_spacing.pro
    ; multiple dot grid images are processed to reduce the noise
    rot_ibis_wl             = [ 0.08 ]

    ; IBIS narrowband plate scales - determined from dot grid using find_dot_grid_spacing.pro
    scale_ibis_7090         = [0.09748, 0.09752, 7090]
    scale_ibis_7699         = [0.09748, 0.09753, 7699]
    scale_ibis_5434         = [0.09707, 0.09711, 5434]
    scale_ibis_nb           = [[scale_ibis_7090],[scale_ibis_7699],[scale_ibis_5434]]
    ; the rotation angle of the whitelight dot grid determined from find_dot_grid_spacing.pro
    ; multiple dot grid images are processed to reduce the noise
    ; the rotation angle is assumed to be the same for all filters/wavelengths
    rot_ibis_nb             = [0.28]

    ; now we can take the difference of the grid angles to determine the relative rotation between the
    ; narrowband and whitelight channels. 
    ; this value will be applied to the science observations as well
    rot_nb_to_wl            = [[time_ref],[rot_ibis_nb - rot_ibis_wl]]
    
    ; values of any integer shifts that should be applied to the narrowband or whitelight data
    ; not used by load_calibrated_image.pro (as of July, 2019)
    shift_ibis_wl_even      = [0, 0]
    shift_ibis_nb_even      = [0, 0]

    ; optical shifts between images through different filters
    ; these are determined by computing the cross correlation between the Air Force target images taken through different filters
    ; one filter is choses as the "reference wavelength," or some average position could be chosen
    ; these shifts are due to optical effects from the (tilted) prefilter or the automatic repositioning of the camera 
    ; to adjust focus for each filter.

    ; original values
    ; shift_target_7090       = [ -1.0,   1.6, 7090 ] 
    ; shift_target_7699       = [  1.0,   0.5, 7699 ] 
    ; shift_target_5434       = [ -4.1,   0.1, 5434 ] 
    ; new revised values, 03 September, 2019
    shift_target_7090       = [ -1.0 + 0.07,   2.15 - 0.1, 7090 ] 
    shift_target_7699       = [  1.38 - 0.06,   0.6 + 0.07, 7699 ] 
    shift_target_5434       = [ -4.1 + 0.14,   0.1 + 0.15, 5434 ] 

    shift_target_nb         = [ [shift_target_7090], [shift_target_7699], [shift_target_5434] ]

    ; optical shifts between narrowband reference wavelength and whitelight channel
    ; this value is also calculated from the cross-correlation between the whitelight target image and the 
    ; target image at the reference wavelength
    shift_wl_to_nb          = [ 4.47, 9.18 ]
    
    ; load_destretch vectors
    IF FILE_TEST(calibration_location + '/destr.components.nb2wl.txt',/Read) THEN BEGIN
        vects_in    = read_ascii(calibration_location + '/destr.components.nb2wl.txt',type='float',record_start=0,data_start=1)
        rdisp_nb2wl = reform(vects_in.field1[0:1,*],2,49,49)
        disp_nb2wl  = reform(vects_in.field1[2:3,*],2,49,49)
    ENDIF ELSE BEGIN
        rdisp_nb2wl = -1
        disp_nb2wl  = -1
    ENDELSE
    
    ; define which filters were in the wheel for a given observing day, and which filters were used
    filter_ids              = ['', '8542', '7090', '6563', '', '7699', '5896', '5434']
    ;filters_used           = [2,5,7]   ; filter wheel order
    filters_used            = [2,5,7]   ; wavelength sampling order
    num_filters             = N_ELEMENTS(filters_used)

    ; this is a placeholder for the future definition of the time-dependent offsets between wavelengths due to
    ; atmospheric dispersion or time-dependent changes in the offset between whitelight and narrowband channels
    ;wl_to_nb_drift          = DBLARR(4,num_timesteps,num_filters) + 1
    wl_to_nb_drift          = [  -1795.28,    2979.13,   14061.74,  -21798.14,  -43832.50,   63662.14,  $
                                 67972.67,  -92785.72,  -52450.36,   67499.80,   16115.69,  -19611.18]
    wl_to_nb_drift          = REFORM(wl_to_nb_drift, 2, 6)

    found_atm_disp_file = FILE_TEST(FILE_SEARCH(calibration_location, 'atmospheric.dispersion.calc.25Apr2019.sav'),/Read)
    IF STRLOWCASE(instrument_channel) EQ 'ibis_nb' THEN BEGIN
        atm_dispersion_file = (FILE_SEARCH(calibration_location, 'atmospheric.dispersion.calc.25Apr2019.sav'))[0]
	IF FILE_TEST(atm_dispersion_file,/Read) THEN BEGIN
	    RESTORE,Verbose=0,atm_dispersion_file
	    ; provides ATM_DISP_CALC and SEQUENCE_TIMES
	    
	    ; pull out refraction values decomposed into x- and y-shifts (in solar heliocentric coordinates)
	    dispcalc_sfts       = [REFORM(atm_disp_calc.SFTS_HELIOCENT_EW,1,1048,4),REFORM(atm_disp_calc.SFTS_HELIOCENT_NS,1,1048,4)]
	    dispcalc_sfts_size  = SIZE(dispcalc_sfts)
	    num_timesteps       = dispcalc_sfts_size[2]
	    ; ATM_DISP_CALC.times_jd also stores the times for the calculated refraction
	    atm_dispersion_times = sequence_times

	    ; take difference between the refraction at the wavelengths of the whitelight image
	    ; and the refraction at each narrowband wavelength to get theoretical dispersion offsets 
	    ; between each of the channels/filters
	    ; the atmospheric dispersion values are stored in increasing wavelength order - [5434,7090,7200,7699]
	    ; 7200 (index 2) is the wavelength of the whitelight images, so all the shifts 
	    ; are differenced with respect to that wavelength
	    atm_dispersion_nb            = FLTARR(3,num_timesteps,num_filters)
	    atm_dispersion_nb[2,*,0]     = atm_disp_calc.wavelengths[1]*10
	    atm_dispersion_nb[2,*,1]     = atm_disp_calc.wavelengths[3]*10
	    atm_dispersion_nb[2,*,2]     = atm_disp_calc.wavelengths[0]*10
	    atm_dispersion_nb[0:1,*,0,0] = dispcalc_sfts[*,*,2] - dispcalc_sfts[*,*,1]
	    atm_dispersion_nb[0:1,*,1,0] = dispcalc_sfts[*,*,2] - dispcalc_sfts[*,*,3]
	    atm_dispersion_nb[0:1,*,2,0] = dispcalc_sfts[*,*,2] - dispcalc_sfts[*,*,0]
	ENDIF ELSE BEGIN
	    num_timesteps              = 300
	    atm_dispersion_nb          = FLTARR(3,2,num_filters)
	    FOR filtn = 0,num_filters-1 DO BEGIN & atm_dispersion_nb[2,*,filtn] = FLOAT(filter_ids[filters_used[filtn]]) & ENDFOR
	    atm_dispersion_times       = [time_ref, time_ref+1]
	ENDELSE
	wl_optical_drifts       = FLTARR(2,num_timesteps)
    ENDIF
    
    ;Precalculated hmi_pixel cutouts for the 
    hmi_pixel_cutout        = FLTARR(1,5)
    hmi_pixel_cutout[0,*]   = [543, 1542, 482, 1481, time_ref] 
    ;hmi_pixel_cutout[0,*]    = [87, 246, 77, 236, time_ref]
    
    ; ---------- ROSA Corrections --------
    ; values are saved here for future use but are not currently returned
    ; -- plate scale and rotation of ROSA images
    scale_rosa_gband  = [0.0, 0.0]
    rot_rosa_gband    = [-0.42]
    scale_rosa_cak    = [0.0, 0.0]
    rot_rosa_cak      = [+0.41]
    scale_rosa_gband  = [0.07571, 0.07748]
    rot_rosa_gband    = [-0.42]
    scale_rosa_cak    = [0.1543,  0.1542]
    rot_rosa_cak      = [+0.41]


; package all the above information into a structure to be returned to the user

rot_to_solnorth_combo = rot_wl_to_sol_north

CASE STRLOWCASE(instrument_channel) OF
    'ibis_wl' : BEGIN
        output_info = CREATE_STRUCT('transpose',          ibis_wl_transpose, $
                                    'rot_to_solnorth',    rot_to_solnorth_combo, $
                                    'rot_to_sol_reftime', rot_wl_to_sol_reftime, $
                                    'rot_to_grid',        rot_ibis_wl, $
                                    'plate_scale',        scale_ibis_wl, $
                                    'uniform_plate_scale', target_plate_scale, $
                                    'shift_even_scale',   shift_ibis_wl_even, $
                                    'optical_shift',      shift_wl_to_nb, $
                                    'hmi_pixel_cutout',   hmi_pixel_cutout, $
                                    'dark_file',          wl_dark_cal_file, $
                                    'dark_name',          wl_dark_name, $
                                    'gain_file',          wl_gain_cal_file, $
                                    'gain_name',          wl_gain_name, $
                                    'bad_pixel_file',     wl_bad_pixel_file)
    END
    'ibis_nb' : BEGIN
        ; add the relative rotation between the narrowband and whitelight channels to 
        ; the rotation-to-solar-north value
        rot_to_solnorth_combo[0]    += rot_nb_to_wl[1]
        output_info = CREATE_STRUCT('transpose',          ibis_nb_transpose, $
                                    'rot_to_solnorth',    rot_to_solnorth_combo, $
                                    'rot_to_sol_reftime', rot_wl_to_sol_reftime, $
                                    'rot_to_grid',        rot_ibis_nb, $
                                    'plate_scale',        scale_ibis_nb, $
                                    'uniform_plate_scale', target_plate_scale, $
                                    'shift_even_scale',   shift_ibis_nb_even,$
                                    'optical_shift',      shift_target_nb, $
                                    'hmi_pixel_cutout',   hmi_pixel_cutout, $
                                    'nb_to_wl_destr_ref', rdisp_nb2wl, $
                                    'nb_to_wl_destr_sft', disp_nb2wl, $
                                    'wl_to_nb_drift',     wl_to_nb_drift, $
                                    'wl_drifts',          wl_optical_drifts, $
                                    'atm_dispersion',     atm_dispersion_nb, $
                                    'atm_dispersion_times',atm_dispersion_times, $
                                    'filter_ids',         filter_ids, $
                                    'dark_file',          nb_dark_cal_file, $
                                    'dark_name',          nb_dark_name, $
                                    'gain_file',          nb_gain_cal_file, $
                                    'gain_name',          nb_gain_name, $
                                    'bluemap_file',       nb_blue_map_file, $
                                    'bluemap_name',       nb_blue_map_name, $
                                    'bad_pixel_file',     nb_bad_pixel_file)
    END
ENDCASE

RETURN,output_info

END
