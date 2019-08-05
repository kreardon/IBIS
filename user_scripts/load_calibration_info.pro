FUNCTION load_calibration_info, date_str, instrument_channel

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
repository_location = File_Dirname(Routine_Filepath(/Either),/Mark)
calibration_location = repository_location + 'calibration_files/'

    ; name of IDL save file containing calibration frames and
    ; name of variable in save file for desired calibration type
    ; the user may save dark and gain files with whatever names desired and update information here
    wl_dark_cal_file      = 'whitelight.darks.30Sep2014.sav'
    wl_dark_name          = 'wl_dark_ave'
    wl_gain_cal_file      = 'whitelight.gains.30Sep2014.sav'
    wl_gain_name          = 'wl_gain'

    nb_dark_cal_file      = 'nb.dark.30Sep2014.ave.sav'
    nb_dark_name          = 'nb_darks'
    nb_7090_gain_cal_file = ['nb.gain.7090.30Sep2014.ave.sav', '7090']
    nb_7090_gain_name     = ['nb_gain_7090_ave', '7090']
    nb_5876_gain_cal_file = ['nb.gain.5876.30Sep2014.ave.sav', '5876']
    nb_5876_gain_name     = ['nb_gain_5876_ave', '5876']
    nb_6563_gain_cal_file = ['nb.gain.6563.30Sep2014.ave.sav', '6563']
    nb_6563_gain_name     = ['nb_gain_6563_ave', '6563']
    nb_8542_gain_cal_file = ['nb.gain.8542.30Sep2014.ave.sav', '8542']
    nb_8542_gain_name     = ['nb_gain_8542_ave', '8542']

    nb_gain_cal_file      = [[nb_5876_gain_cal_file], [nb_6563_gain_cal_file], [nb_7090_gain_cal_file], [nb_8542_gain_cal_file]]
    nb_gain_name          = [[nb_5876_gain_name], [nb_6563_gain_name], [nb_7090_gain_name], [nb_8542_gain_name]]

    ; pick a reference time for some time-dependent values
    ; typically pick the time of midnight before the start of the observations
    time_ref                = (fits_date_convert('2014-09-30T00:00:00'))[0]

    ; this is how a fixed value for rotation would be defined
    ;rot_wl_to_sol_north     = [0.32]    ; determined from fitting to HMI
    ; this is how a linear trend to rotation would be defined
    ; where the independent variable is fractional days since the reference time
    rot_wl_to_sol_north     = [0.098793436,0.29216892]    ; determined from fitting to HMI
    ; reference time is midnight in UT
    rot_wl_to_sol_reftime   = time_ref

    ; IBIS narrowband plate scales - determined from dot grid using find_dot_grid_spacing.pro
    ;scale_ibis_wl          = [ 0.09074, 0.09295 ]  ; determined from dot grid
    ; IBIS narrowband plate scales - determined from fitting to HMI using stx_findlocation.pro
    ; this generally requires running alignment over many whitelight vs. HMI combinations to beat down the noise
    ; but in this case if found a plate scale that was about 0.2% larger
    scale_ibis_wl           = [ 0.09099, 0.09312 ]   ; determined from fitting to HMI, about 0.2% larger
    ; the rotation angle of the whitelight dot grid determined from find_dot_grid_spacing.pro
    ; multiple dot grid images are processed to reduce the noise
    rot_ibis_wl             = [ 0.57 ]

    ; IBIS narrowband plate scales - determined from dot grid using find_dot_grid_spacing.pro
    scale_ibis_8542         = [0.09539, 0.09749, 8542]
    scale_ibis_7090         = [0.09544, 0.09751, 7090]
    scale_ibis_6563         = [0.09533, 0.09742, 6563]
    scale_ibis_5876         = [0.09494, 0.09707, 5876]
    scale_ibis_nb           = [[scale_ibis_5876],[scale_ibis_6563],[scale_ibis_7090],[scale_ibis_8542]]
    ; the rotation angle of the whitelight dot grid determined from find_dot_grid_spacing.pro
    ; multiple dot grid images are processed to reduce the noise
    ; the rotation angle is assumed to be the same for all filters/wavelengths
    rot_ibis_nb             = [0.17]

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
    ; in this case, the H-alpha filter was chosen, since it's position was most "central"
    ; (i.e. the other filters had shifts on either side of the H-alpha position
    ; these shifts are due to optical effects from the (tilted) prefilter or the automatic repositioning of the camera 
    ; to adjust focus for each filter.
    shift_target_8542       = [ -2.8,   3.2, 8542 ] 
    shift_target_7090       = [ -1.1,   3.8, 7090 ] 
    shift_target_6563       = [ -0.0,  -0.0, 6563 ] 
    shift_target_5876       = [  4.3,  -3.6, 5876 ] 
    shift_target_nb         = [ [shift_target_5876], [shift_target_6563], [shift_target_7090], [shift_target_8542] ]

    ; optical shifts between narrowband reference wavelength and whitelight channel
    ; this value is also calculated from the cross-correlation between the whitelight target image and the 
    ; target image at the reference wavelength
    shift_wl_to_nb          = [19.6, -17.6]
    
    ; define which filters were in the wheel for a given observing day, and which filters were used
    filter_ids              = ['', '8542', '7773', '6563', '', '7699', '7090', '5876']
    ;filters_used           = [1,3,6,7]   ; filter wheel order
    filters_used            = [6,3,1,7]   ; wavelength sampling order
    num_filters             = N_ELEMENTS(filters_used)

    ; this is a placeholder for the future definition of the time-dependent offsets between wavelengths due to
    ; atmospheric dispersion or time-dependent changes in the offset between whitelight and narrowband channels
    num_timesteps           = 300
    wl_to_nb_drift          = DBLARR(4,num_timesteps,num_filters) + 1
    wl_optical_drifts       = FLTARR(2,num_timesteps) + 1
    atm_dispersion_nb       = FLTARR(2,num_timesteps,num_filters) + 1

    ;Precalculated hmi_pixel cutouts for the 
    hmi_pixel_cutout        = FLTARR(1,5)
    hmi_pixel_cutout[0,*]   = [346,1345,359,1358, time_ref] 

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
                                    'shift_even_scale',   shift_ibis_wl_even, $
                                    'optical_shift',      shift_wl_to_nb, $
                                    'hmi_pixel_cutout',   hmi_pixel_cutout, $
                                    'dark_file',          wl_dark_cal_file, $
                                    'dark_name',          wl_dark_name, $
                                    'gain_file',          wl_gain_cal_file, $
                                    'gain_name',          wl_gain_name)
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
                                    'shift_even_scale',   shift_ibis_nb_even,$
                                    'optical_shift',      shift_target_nb, $
                                    'wl_to_nb_drift',     wl_to_nb_drift, $
                                    'wl_drifts',          wl_optical_drifts, $
                                    'atm_dispersion',     atm_dispersion_nb, $
                                    'filter_ids',         filter_ids, $
                                    'dark_file',          nb_dark_cal_file, $
                                    'dark_name',          nb_dark_name, $
                                    'gain_file',          nb_gain_cal_file, $
                                    'gain_name',          nb_gain_name)
    END
ENDCASE

RETURN,output_info

END
