FUNCTION load_calibration_info, date_str, instrument_channel

; sample usage:
; nb_21apr_info = load_alma_calibration_info('21Apr2017','ibis_nb') 

ibis_wl_transpose    =  3
ibis_nb_transpose    =  7
;rosa_gband_transpose =  6
;rosa_cak_transpose   =  1

target_plate_scale   = 0.096

repository_location = File_Dirname(Routine_Filepath(/Either),/Mark)
calibration_location = repository_location + 'calibration_files/'

   
;'21apr2017' 

    wl_dark_cal_file = 'whitelight.darks.30Sep2014.sav'
    wl_dark_name     = 'wl_dark_ave'
    wl_gain_cal_file = 'whitelight.gains.30Sep2014.sav'
    wl_gain_name     = 'wl_gain'

    nb_dark_cal_file = 'nb.dark.30Sep2014.ave.sav'
    nb_dark_name     = 'nb_darks'
    nb_7090_gain_cal_file = ['nb.gain.7090.30Sep2014.ave.sav', 7090]
    nb_7090_gain_name     = ['nb_gain_7090_ave', 7090]
    nb_5876_gain_cal_file = ['nb.gain.5876.30Sep2014.ave.sav', 5876]
    nb_5876_gain_name     = ['nb_gain_5876_ave', 5876]
    nb_6563_gain_cal_file = ['nb.gain.6563.30Sep2014.ave.sav', 6563]
    nb_6563_gain_name     = ['nb_gain_6563_ave', 6563]
    nb_8542_gain_cal_file = ['nb.gain.8542.30Sep2014.ave.sav', 8542]
    nb_8542_gain_name     = ['nb_gain_8542_ave', 8542]

    nb_gain_cal_file  = [[nb_5876_gain_cal_file], [nb_6563_gain_cal_file], [nb_7090_gain_cal_file], [nb_8542_gain_cal_file]]
    nb_gain_name      = [[nb_5876_gain_name], [nb_6563_gain_name], [nb_7090_gain_name], [nb_8542_gain_name]]

    num_minutes   = 300
    time_start    = 14.5
    time_axis_min = (DINDGEN(num_minutes)/60.) + time_start
    time_ref      = (fits_date_convert('2014-09-30T18:00:00'))[0]

    rot_grid_to_sol_north = [[time_ref],[-0.249]]
    rot_wl_to_sol_north   = [[time_ref],[0.321]]    ; determined from fitting to HMI

    ;scale_ibis_wl     = [ 0.09668, 0.09680 ]
    ;scale_ibis_wl     = [ 0.09074, 0.09295 ]  ; determined from dot grid
    scale_ibis_wl     = [ 0.09099, 0.09312 ]   ; determined from fitting to HMI, about 0.2% larger
    rot_ibis_wl       = [ 0.57 ]

    scale_ibis_8542   = [0.09539, 0.09749, 8542]
    scale_ibis_7090   = [0.09544, 0.09751, 7090]
    scale_ibis_6563   = [0.09533, 0.09742, 6563]
    scale_ibis_5876   = [0.09494, 0.09707, 5876]
    scale_ibis_nb     = [[scale_ibis_5876],[scale_ibis_6563],[scale_ibis_7090],[scale_ibis_8542]]
    rot_ibis_nb       = [ 0.18 ]
    
    shift_ibis_wl_even = [0, 0]
    shift_ibis_nb_even = [0, 0]  

    shift_target_8542 = [ -0.0,  -0.0, 8542 ] 
    shift_target_7090 = [ -0.0,   0.0, 7090 ] 
    shift_target_6563 = [ -0.0,  -0.0, 6563 ] 
    shift_target_5876 = [ -0.0,  -0.0, 5876 ] 
    shift_target_nb   = [ [shift_target_5876], [shift_target_6563], [shift_target_7090], [shift_target_8542] ]
    
    filter_ids        = ['', '8542', '7773', '6563', '', '7699', '7090', '5876']
     ;filters_used     = [1,3,6,7]   ; filter wheel order
    filters_used      = [6,3,1,7]   ; wavelength sampling order
    num_filters       = N_ELEMENTS(filters_used)

    wl_to_nb_drift    =  DBLARR(4,num_minutes,num_filters)
    wl_optical_drifts = FLTARR(2,num_minutes)
    atm_dispersion_nb = FLTARR(2,num_minutes,num_filters)

    ;Precalculated hmi_pixel cutouts for the 
    hmi_pixel_cutout_target_1 = [346,1345,359,1358] 
    hmi_pixel_cutout_target_2 = [375,1374,324,1323] 
    hmi_pixel_cutout          = FLTARR(2,4)
    hmi_pixel_cutout[0,*]     = hmi_pixel_cutout_target_1
    hmi_pixel_cutout[1,*]     = hmi_pixel_cutout_target_2

    ; ---------- ROSA Corrections --------
    ; -- plate scale and rotation of ROSA images
    scale_rosa_gband  = [0.0, 0.0]
    rot_rosa_gband    = [-0.42]
    scale_rosa_cak    = [0.0, 0.0]
    rot_rosa_cak      = [+0.41]
    scale_rosa_gband  = [0.07571, 0.07748]
    rot_rosa_gband    = [-0.42]
    scale_rosa_cak    = [0.1543,  0.1542]
    rot_rosa_cak      = [+0.41]


rot_to_solnorth_combo = rot_grid_to_sol_north

CASE STRLOWCASE(instrument_channel) OF
    'ibis_wl' : BEGIN
        rot_to_solnorth_combo[*,1] += rot_ibis_wl[0]
        output_info = CREATE_STRUCT('transpose',        ibis_wl_transpose, $
                                    'rot_to_solnorth',  rot_to_solnorth_combo, $
                                    'rot_to_grid',      rot_ibis_wl, $
                                    'plate_scale',      scale_ibis_wl, $
                                    'shift_even_scale', shift_ibis_wl_even, $
                                    'hmi_pixel_cutout', hmi_pixel_cutout, $
                                    'dark_file',        wl_dark_cal_file, $
                                    'dark_name',        wl_dark_name, $
                                    'gain_file',        wl_gain_cal_file, $
                                    'gain_name',        wl_gain_name)
    END
    'ibis_nb' : BEGIN
        rot_to_solnorth_combo[*,1] += rot_ibis_nb[0]
        output_info = CREATE_STRUCT('transpose',       ibis_nb_transpose, $
                                    'rot_to_solnorth', rot_to_solnorth_combo, $
                                    'rot_to_grid',     rot_ibis_nb, $
                                    'plate_scale',     scale_ibis_nb, $
                                    'shift_even_scale',shift_ibis_nb_even,$
                                    'optical_shift',   shift_target_nb, $
                                    'wl_to_nb_drift',  wl_to_nb_drift, $
                                    'wl_drifts',       wl_optical_drifts, $
                                    'atm_dispersion',  atm_dispersion_nb, $
                                    'filter_ids',      filter_ids, $
                                    'dark_file',       nb_dark_cal_file, $
                                    'dark_name',       nb_dark_name, $
                                    'gain_file',       nb_gain_cal_file, $
                                    'gain_name',       nb_gain_name)
    END
ENDCASE

RETURN,output_info

END
