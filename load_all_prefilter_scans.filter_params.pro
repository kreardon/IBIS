
; ----------------------------------------------------------------------
;  This file contains all the filter-specific setting needed in the 
;  processing carried out in "load_all_prefilter_scans.v2.pro"
;
;  This program is called by load_all_prefilter_scans.v2.pro and should 
;  be in the same directory as that program (or in the working directory?).
; ----------------------------------------------------------------------

params_quality_metric       = FLTARR(2,6)

PRINT, ' --> Loading specific settings for ' + filter_name + ' filter'

CASE filter_name OF
'5896' : BEGIN
    ; prefilter scan acceptability ranges for different metrics
    params_quality_metric(*,1)  = [2.5e4, 2.5e5]
    params_quality_metric(*,3)  = [1.04, 1.15]
    params_quality_metric(*,4)  = [0.99, 1.04]
    params_quality_metric(*,5)  = [0.0, 0.018]

    ; parameters for fringe averaging and phase evaluation

    ; for each prefilter we need to manually pick out:
    ;     [fringe_stwv, fringe_ndwv] = the optimal spectral range over which to sum the fringes
    ;     [fringe_stscn, fringe_ndscn] = the desired indexes (i.e. time) over which to sum the fringes
    ;     phase_ref = the frequency bins over which to compute the phase shift of the fringes
    fringe_stwv  = 101
    fringe_ndwv  = 199
    fringe_wvlen = fringe_ndwv - fringe_stwv + 1
    fringe_stscn = 62
    fringe_ndscn = 121
    phase_ref    = [5]

    ; paramaters for sinusoidal fringe model
    fringe_freq = 0.393
    fringe_off  = 0.047
    fringe_amp  = 0.01307

    ; the array of scan indexes (it doesn't necessarily have to be continuous) 
    ; over which to compute the final averaged/median prefilter profile
    pref_sum_range = INDGEN(350) + 149
    ; set to zero to use all valid profiles in the average/median
    ;pref_sum_range = 0
    pref_sum_range_initial = INDGEN(prefilter_scan_good_num) + 0

END
'6173' : BEGIN
    params_quality_metric(*,1)  = [2.0e4, 2.5e5]
    params_quality_metric(*,3)  = [1.01, 1.15]
    params_quality_metric(*,4)  = [0.95, 1.04]
    params_quality_metric(*,5)  = [0.0, 0.025]

    fringe_stwv  = 99
    fringe_ndwv  = 245
    fringe_wvlen = fringe_ndwv - fringe_stwv + 1
    fringe_stscn = 438
    fringe_ndscn = 488
    phase_ref    = [7]

    fringe_freq = 0.42
    fringe_off  = 0.2
    fringe_amp  = 0.01

    pref_sum_range_initial = INDGEN(prefilter_scan_good_num) + 0
    ;pref_sum_range_final   = INDGEN(450) + 85
    ; for the 2013 scans
    pref_sum_range_final   = INDGEN(120) + 368

END
'6563' : BEGIN
    params_quality_metric(*,1)  = [2.0e4, 2.e5]
    params_quality_metric(*,3)  = [1.01, 1.15]
    params_quality_metric(*,4)  = [0.95, 1.04]
    params_quality_metric(*,5)  = [0.0, 0.025]

    fringe_freq = 0.42
    fringe_off  = 0.17
    fringe_amp  = 0.0

    pref_sum_range_initial = INDGEN(prefilter_scan_good_num) + 0
    pref_sum_range_final   = INDGEN(prefilter_scan_good_num) + 0

END
'7090' : BEGIN
    params_quality_metric(*,1)  = [2.0e4, 2.e5]
    params_quality_metric(*,3)  = [1.01, 1.15]
    params_quality_metric(*,4)  = [0.95, 1.04]
    params_quality_metric(*,5)  = [0.0, 0.025]

    fringe_stwv  = 137
    fringe_ndwv  = 256
    fringe_wvlen = fringe_ndwv - fringe_stwv + 1
    fringe_stscn = 200
    fringe_ndscn = 399
    phase_ref    = [10,11]

    fringe_freq = 0.42
    fringe_off  = 0.17
    fringe_amp  = 0.0

    pref_sum_range_initial = INDGEN(460) + 176
    pref_sum_range_final   = INDGEN(460) + 176

END
'8542' : BEGIN
    params_quality_metric(*,1)  = [6.0e3, 2.e5]
    params_quality_metric(*,3)  = [1.01, 1.15]
    params_quality_metric(*,4)  = [0.95, 1.04]
    params_quality_metric(*,5)  = [0.0, 0.025]

    fringe_stwv  = 237
    fringe_ndwv  = 356
    fringe_wvlen = fringe_ndwv - fringe_stwv + 1
    fringe_stscn = 375
    fringe_ndscn = 474
    phase_ref    = [12,13]

    fringe_freq = 0.194
    fringe_off  = 0.04
    fringe_amp  = 0.000

    ;pref_sum_range = INDGEN(470) + 370
    pref_sum_range_final = 0

END
ELSE: BEGIN
    params_quality_metric(*,1)  = [1.0e4, 2.5e5]
    params_quality_metric(*,3)  = [1.01, 1.15]
    params_quality_metric(*,4)  = [0.95, 1.05]
    params_quality_metric(*,5)  = [0.0, 0.035]

    fringe_stwv  = 99
    fringe_ndwv  = 245
    fringe_wvlen = fringe_ndwv - fringe_stwv + 1
    fringe_stscn = 341
    fringe_ndscn = 488
    phase_ref    = [5,10]

    fringe_freq = 0.4
    fringe_off  = 0.0
    fringe_amp  = 0.01

    pref_sum_range_initial = INDGEN(prefilter_scan_good_num) + 0
    pref_sum_range_final   = INDGEN(prefilter_scan_good_num) + 0
    
END
ENDCASE
