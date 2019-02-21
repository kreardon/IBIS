FUNCTION load_alma_calibration_info, date_str, instrument_channel

; nb_21apr_info = load_alma_calibration_info('21Apr2017','ibis_nb') 

ibis_wl_transpose    =  3
ibis_nb_transpose    =  7
rosa_gband_transpose =  6
rosa_cak_transpose   =  1

target_plate_scale   = 0.096

repository_location = File_Dirname(Routine_Filepath(/Either),/Mark)
calibration_location = repository_location + 'calibration_files/'

CASE STRLOWCASE(date_str) OF
'23apr2017' : BEGIN

    num_min = 310
    time_start    = 14.
    time_axis_min = (DINDGEN(num_min)/60.) + time_start
    time_ref      = (fits_date_convert('2017-04-23T00:00:00'))[0]

    
    ; -- angle between data reference direction (defined by grid) and solar north (for alignment to HMI).
    ;rot_grid_to_sol_north = -0.48
    ; -- need to include time rangges on rotation angle to account for time varying rotation changes 
    ;rot_grid_to_sol_north = [ 0.31, fits_date_convert('2017-04-23 14:00')]
    ;rot_grid_to_sol_north = [[rot_grid_to_sol_north],[ 0.31, fits_date_convert('2017-04-23 15:05')]]
    ;rot_grid_to_sol_north = [[rot_grid_to_sol_north],[-0.48, fits_date_convert('2017-04-23 15:10')]]
    ;rot_grid_to_sol_north = [[rot_grid_to_sol_north],[-0.48, fits_date_convert('2017-04-23 19:00')]]
    
    ; -- from a first iteration of the destretching from WL to HMI, we have the full measurements
    ;    of the rotation trends with time, which we are able to fit.
    ; -- the trends were fitted using the fractional hours after midnight as the free parameter 
    ;    in order to reducethe effects of numerical errors
;    rot_target1_param      = [2.71, -4.40]
;    rotation_trend_target1 = [(DINDGEN(70)/60. + 14)/24.]
;    rotation_trend_target1 = [[rotation_trend_target1],[rot_target1_param[0] + rotation_trend_target1 * rot_target1_param[1]]]

; -- because the earlier destretching assumed a fixed rotation (different for target 1 and target 2),
    ;    those were not included in the fitted rotation angles, so  we need to correct for the 
    ;    previously applied rotation angles
    ;rotation_hlptrend_target1[*,1]  = -0.2 - rotation_trend_target1[*,1]
    ;rotation_trend_target2[*,1]  = -0.48 - rotation_trend_target2[*,1]

    ; -- latest fit to measured rotation
    ; -- an arbitrary shift was applied to the target 1 rotations (which will be backed out below) 
    ;    so that a single fit will suffice for target 1 and target 2 (for simplicity?)
    rot_target12_param      = [11.703 + 0.71, -1.38517, 0.0384484]
    ;rot_target12_param      = [7.3643, -0.84191, 0.0236297]

    rotation_trend_target12 = time_axis_min
    rotation_trend_target12 = [[rotation_trend_target12],[                            rot_target12_param[0] + $
                                                          rotation_trend_target12   * rot_target12_param[1] + $
                                                          rotation_trend_target12^2 * rot_target12_param[2]]]

    ; -- here we apply an ad hoc shift to the target 1 rotations to return them to the
    ;    the values that were originally measured (before the adjustment to allow for a single fitted function).
    target_change_idx                               = MAX(WHERE(rotation_trend_target12[*,0] LE (15 + 10/60.)))
    rotation_trend_target12[0:target_change_idx,1] += 0.26
    rotation_trend_target12[*,1] *= -1

    rot_grid_to_sol_north           = rotation_trend_target12
    ; -- return x-axis units to Julian day (not fractional hours after midnight)
    rot_grid_to_sol_north[*,0]     = rot_grid_to_sol_north[*,0]/24. + 2457866.5d

    ; -- plate scale of WL images
    scale_ibis_wl   = [0.09448, 0.09681]
    ; -- angle of rotation between WL array axes and dot grid directions
    rot_ibis_wl     = [ -0.62 ]
    
    ; ---------- Narrowband Corrections --------
    
    filter_ids       = ['', '8542', '7773', '6563', '', '7699', '5896', '6173']
    ;filters_used     = [1,3,6,7]   ; filter wheel order
    filters_used     = [6,7,3,1]   ; wavelength order
    num_filters      = N_ELEMENTS(filters_used)

    ; -- plate scale of NB images - varies with each prefilter
    scale_ibis_5896 = [0.09531, 0.09738, 5896]
    scale_ibis_6173 = [0.09538, 0.09745, 6173]
    scale_ibis_6563 = [0.09544, 0.09752, 6563]
    scale_ibis_8542 = [0.09554, 0.09761, 8542]
    scale_ibis_nb   = [[scale_ibis_5896],[scale_ibis_6173],[scale_ibis_6563],[scale_ibis_8542]]
    rot_ibis_nb     = [ -0.87 ]
    
    ; -- internal optical offsets between NB and WL channels - as determined from target images
    shift_target_5896 = [-2.045, -0.1611, 5896]
    shift_target_6173 = [-2.824,  0.0531, 6173]
    shift_target_6563 = [-3.743,  0.0290, 6563]
    shift_target_8542 = [-7.981,  4.566,  8542]
    shift_target_nb   = [[shift_target_5896],[shift_target_6173],[shift_target_6563],[shift_target_8542]]
    
    ; atmospheric dispersion shifts
    RESTORE,calibration_location + 'atmospheric.refraction.calc.23Apr2017.sav'
    wl_wave_idx = get_closest(atm_refraction.wavelengths,630) 
    atm_dispersion_nb = FLTARR(2,num_min,num_filters)
    
    FOR filt_n = 0,num_filters-1 DO BEGIN
        filt_id     = filter_ids[filters_used[filt_n]]
        filt_atm_idx = get_closest(atm_refraction.wavelengths,FLOAT(filt_id)/10.)
        atm_dispersion_ns_full = atm_refraction.SFTS_HELIOCENT_NS[*,wl_wave_idx] - atm_refraction.SFTS_HELIOCENT_NS[*,filt_atm_idx]
        atm_dispersion_ns_min  = INTERPOL(atm_dispersion_ns_full,(atm_refraction.times_jd - time_ref)*24, time_axis_min)
        atm_dispersion_ew_full = atm_refraction.SFTS_HELIOCENT_EW[*,wl_wave_idx] - atm_refraction.SFTS_HELIOCENT_EW[*,filt_atm_idx]
        atm_dispersion_ew_min  = INTERPOL(atm_dispersion_ew_full,(atm_refraction.times_jd - time_ref)*24, time_axis_min)
        
        atm_dispersion_nb[0,*,filt_n] = atm_dispersion_ew_min
        atm_dispersion_nb[1,*,filt_n] = atm_dispersion_ns_min
    ENDFOR
    
    ; additional time-dependent optical shifts (due to beam wobble?)
    ; -- These drifts have been calculated by comparing narrowband and whitelight images and 
    ;    then subtracting off the atmospheric-dispersion-induced shifts
    ; -- This was done for both H-alpha and Na 5896 which have good granulation images in the wings
    ;    to check validity, but we use only drifts measured from H-alpha because they are 
    ;    more numerous and cover the full time range of the observations
    RESTORE,calibration_location + 'wl.optical.drifts.vs.nb.23Apr2017.sav'
    wl_optical_drifts = FLTARR(2,num_min)
    wldrift_ha_ew_smooth = SMOOTH(MEDIAN(wldrift_ha_ew_orig,9),15,/EDGE_TRUNC)
    wldrift_ha_ns_smooth = SMOOTH(MEDIAN(wldrift_ha_ns_orig,9),15,/EDGE_TRUNC)
    wl_optical_drifts[0,*] = interpol(wldrift_ha_ew_smooth,(wldrift_halpha_times-time_ref)*24,time_axis_min,/LSQUAD)
    wl_optical_drifts[1,*] = interpol(wldrift_ha_ns_smooth,(wldrift_halpha_times-time_ref)*24,time_axis_min,/LSQUAD)
    
    ; combine optical drifts and atmospheric dispersion into a single array with 
    ; both time dependent offsets. 
    ; convert from arcseconds back to nominal pixel units for convenience
    ; output array will have 3 dimensions: [4, time, filter]
    ; output array with have four values per time and per filter:
    ; [ x-direction shift, y-direction shift, filter id/wavelength, Julian date ]
    wl_to_nb_drift =  DBLARR(4,num_min,num_filters)
    FOR filt_n = 0,num_filters-1 DO BEGIN
        wl_to_nb_drift[0:1,*,filt_n] = (wl_optical_drifts - atm_dispersion_nb[*,*,filt_n]) / target_plate_scale
        wl_to_nb_drift[2,*,filt_n]   = FLOAT(filter_ids[filters_used[filt_n]])
        wl_to_nb_drift[3,*,filt_n]   = time_axis_min/24.d + time_ref
    ENDFOR
   
   ; ---------- ROSA Corrections --------

    
    ; -- plate scale and rotation of ROSA images
    scale_rosa_gband = [0.07571, 0.07748]
    rot_rosa_gband   = [-0.42]
    scale_rosa_cak   = [0.1543,  0.1542]
    rot_rosa_cak     = [+0.41]
    
    
    shift_ibis_wl_even = [0, -21]
    shift_ibis_nb_even = [0, 0]  
    ;Precalculated hmi_pixel cutouts for the 
    hmi_pixel_cutout_target_1 = [346,1345,359,1358] 
    hmi_pixel_cutout_target_2 = [375,1374,324,1323] 
    hmi_pixel_cutout          = FLTARR(2,4)
    hmi_pixel_cutout[0,*]     = hmi_pixel_cutout_target_1
    hmi_pixel_cutout[1,*]     = hmi_pixel_cutout_target_2
    
    END    

'22apr2017' : BEGIN
    rot_grid_to_sol_north = -0.48

    scale_ibis_wl   = [0.09447, 0.09679]
    rot_ibis_wl     = [ -0.57 ]

    scale_ibis_8542 = [0.09553, 0.09760, 8542]
    scale_ibis_6563 = [0.09544, 0.09750, 6563]
    scale_ibis_7699 = [0.09552, 0.09756, 7699]
    scale_ibis_5896 = [0.09531, 0.09736, 5896]
    scale_ibis_nb   = [[scale_ibis_5896],[scale_ibis_6563],[scale_ibis_7699],[scale_ibis_8542]]
    rot_ibis_nb     = [ -0.82 ]

    ; derived from 20170422_192512 target series    
    shift_target_8542 = [  -9.190,   3.170, 8542 ] 
    shift_target_6563 = [  -5.880,  -0.994, 6563 ] 
    shift_target_7699 = [  -7.353,   7.156, 7699 ] 
    shift_target_5896 = [  -3.740,  -0.872, 5896 ] 
    shift_target_nb   = [[shift_target_5896],[shift_target_6563],[shift_target_7699],[shift_target_8542]]

    scale_rosa_gband = [0.0, 0.0]
    rot_rosa_gband   = [-0.42]
    scale_rosa_cak   = [0.0, 0.0]
    rot_rosa_cak     = [+0.41]
    filter_ids       = ['', '8542', '7773', '6563', '', '7699', '5896', '6173']
    
    shift_ibis_wl_even = [0, -21]

    END    
'21apr2017' : BEGIN
    num_min       = 300
    time_start    = 13.5
    time_axis_min = (DINDGEN(num_min)/60.) + time_start
    time_ref      = (fits_date_convert('2017-04-21T15:00:00'))[0]

    rot_grid_to_sol_north = [[time_ref],[-0.48]]

    scale_ibis_wl     = [ 0.09668, 0.09680 ]
    rot_ibis_wl       = [ -0.72 ]

    scale_ibis_8542   = [0.09775, 0.09761, 8542]
    scale_ibis_7699   = [0.09775, 0.09758, 7699]
    scale_ibis_5896   = [0.09752, 0.09736, 5896]
    scale_ibis_5434   = [0.09738, 0.09721, 5434]
    scale_ibis_nb     = [[scale_ibis_5434],[scale_ibis_5896],[scale_ibis_7699],[scale_ibis_8542]]
    rot_ibis_nb       = [ -0.62 ]

    shift_target_8542 = [ -26.310,  -3.739, 8542 ] 
    shift_target_7699 = [ -27.416,   0.146, 7699 ] 
    shift_target_5896 = [ -30.538,  -6.323, 5896 ] 
    shift_target_5434 = [ -30.803,  -5.496, 5434 ] 
    shift_target_nb   = [ [shift_target_8542], [shift_target_7699], [shift_target_5896], [shift_target_5434] ]
    
    filter_ids        = ['', '8542', '7773', '6563', '', '7699', '5896', '5434']
     ;filters_used     = [1,5,6,7]   ; filter wheel order
    filters_used      = [7,6,5,1]   ; wavelength order
    num_filters       = N_ELEMENTS(filters_used)

    wl_to_nb_drift    =  DBLARR(4,num_min,num_filters)
    wl_optical_drifts = FLTARR(2,num_min)
    atm_dispersion_nb = FLTARR(2,num_min,num_filters)
    
    shift_ibis_wl_even = [0, -21]
    shift_ibis_nb_even = [0, 0]  

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

    END    

'19apr2017' : BEGIN
    num_min       = 240
    time_start    = 13.5
    time_axis_min = (DINDGEN(num_min)/60.) + time_start
    time_ref      = (fits_date_convert('2017-04-19T15:00:00'))[0]
    
    rot_grid_to_sol_north = [[time_ref],[-0.48]]

    scale_ibis_wl   = [0.09461, 0.09686]
    rot_ibis_wl     = [ -1.05 ]

    scale_ibis_8542 = [0.09552, 0.09759, 8542]
    scale_ibis_7699 = [0.09550, 0.09756, 7699]
    scale_ibis_5896 = [0.09530, 0.09737, 5896]
    scale_ibis_6563 = [0.09541, 0.09748, 6563]
    scale_ibis_nb   = [[scale_ibis_6563],[scale_ibis_5896],[scale_ibis_7699],[scale_ibis_8542]]
    rot_ibis_nb     = [ -0.96 ]

    shift_target_8542 = [ -22.059,   4.701, 8542 ] 
    shift_target_6563 = [ -18.490,   3.195, 6563 ] 
    shift_target_7699 = [ -20.293,   8.785, 7699 ] 
    shift_target_5896 = [ -15.527,   1.583, 5896 ] 
    shift_target_nb   = [ [shift_target_8542], [shift_target_6563], [shift_target_7699], [shift_target_5896] ]
    
    scale_rosa_gband = [0.0, 0.0]
    rot_rosa_gband   = [-0.42]
    scale_rosa_cak   = [0.0, 0.0]
    rot_rosa_cak     = [+0.41]
    filter_ids       = ['', '8542', '7773', '6563', '', '7699', '5896', '5434']
    ;filters_used     = [1,3,5,6]   ; filter wheel order
    filters_used     = [6,3,5,1]   ; wavelength order
    num_filters      = N_ELEMENTS(filters_used)

    wl_to_nb_drift =  DBLARR(4,num_min,num_filters)
    wl_optical_drifts = FLTARR(2,num_min)
    atm_dispersion_nb = FLTARR(2,num_min,num_filters)
    
    shift_ibis_wl_even = [0, -21]
    shift_ibis_nb_even = [0, 0]  
 
    ;Precalculated hmi_pixel cutouts for the 
    hmi_pixel_cutout_target_1 = [346,1345,359,1358] 
    hmi_pixel_cutout_target_2 = [375,1374,324,1323] 
    hmi_pixel_cutout          = FLTARR(2,4)
    hmi_pixel_cutout[0,*]     = hmi_pixel_cutout_target_1
    hmi_pixel_cutout[1,*]     = hmi_pixel_cutout_target_2
   
   ; ---------- ROSA Corrections --------

    ; -- plate scale and rotation of ROSA images
    scale_rosa_gband = [0.07571, 0.07748]
    rot_rosa_gband   = [-0.42]
    scale_rosa_cak   = [0.1543,  0.1542]
    rot_rosa_cak     = [+0.41]

    END    

ENDCASE

rot_to_solnorth_combo = rot_grid_to_sol_north

CASE STRLOWCASE(instrument_channel) OF
    'ibis_wl' : BEGIN
        rot_to_solnorth_combo[*,1] += rot_ibis_wl[0]
        output_info = CREATE_STRUCT('transpose',        ibis_wl_transpose, $
                                    'rot_to_solnorth',  rot_to_solnorth_combo, $
                                    'rot_to_grid',      rot_ibis_wl, $
                                    'plate_scale',      scale_ibis_wl, $
                                    'shift_even_scale', shift_ibis_wl_even, $
                                    'hmi_pixel_cutout', hmi_pixel_cutout)
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
                                    'filter_ids',      filter_ids)
    END
ENDCASE

RETURN,output_info

END
