;-------------------------------------------------------------
;+
; NAME:
;       prepare_prefilter_scans
; PURPOSE:
;       
; CATEGORY:
; CALLING SEQUENCE:
;       prefilter_scan_even = prepare_prefilter_scans (prefilter_scan_files)
; INPUTS:
;       prefilter_scan_files = an array of one or more prefilter scan log files
;                              for a single prefilter
; KEYWORD PARAMETERS:
;       verbose = turns on prining of each log line for the requested prefilter scan  
;
;       params_quality_metric = paramaters to tune the acceptable ranges of several 
;                               quality metrics used to select "valid" prefilter scans 
; OUTPUTS:
;       prefilter_scan_even = structure containing all unique and valid prefilter scans
;                             extracted from the input files, interpolated onto a common
;                             wavelength grid for easy comparison.
;                             see details on structure contents in NOTES section.
; OPTIONAL OUTPUTS:
;       uniq_prefilter_scans = structure containing all unique prefilter scans
;                             extracted from the input files, stored with the wavelength 
;                             grid in which they were recorded.
;       good_prefilter_scans = structure containing all unique and valid prefilter scans
;                             extracted from the input files, stored with the wavelength 
;                             grid in which they were recorded.
; COMMON BLOCKS:
;       
; NOTES:
;       The input to this function is an array of one or more prefilter scan files. For a 
;       given filter these files can be determined with commands like this:
;           IDL> filter = 6173
;           IDL> prefscan_all_files = file_search('~/IBIS','sintodata.' + STRTRIM(filter,2) + '*', $
;                                            COUNT=num_files)
;
;       IDL> HELP, prefilter_scan_even, /STRUCTURE
;       ** Structure <7078e0>, 10 tags, length=1544, data length=1544, refs=1:
;          COUNTS          - the measured intensity at each given Fabry-Perot tuning
;          WAVELENGTHS     - the wavelengths (relative to the determined peak of given filter)
;                            corresponding to each measured intensity
;          FP1_VOLTAGES    - the voltage applied to Fabry-Perot #1 for each given measurement
;          FP2_VOLTAGES    - the voltage applied to Fabry-Perot #2 for each given measurement
;          SCAN_LOCATION   - starting position (from POINT_LUN) for prefilter scan
;          SCAN_START_DATE - starting date of prefilter scan as a Julian date
;          SCAN_END_DATE   - ending date of prefilter scan as a Julian date
;          SCAN_DATE_TEXT  - the starting and ending dates given in original format
;          LOGFILE_USED    - name of the log file to which the parsing
;          FILTER_NAME     - input name of filter
;
;       COUNTS and FP[12]_VOLTAGES are interpolated from their original values, and may not be exact.
;
;
; MODIFICATION HISTORY:
;       26 Apr, 2015 - KPR - initial implementation in a standalone function
;
;-
;-------------------------------------------------------------
FUNCTION prepare_prefilter_scans, prefscan_all_files, verbose=verbose, wvstep_even=wvstep_even, $
                                  good_prefilter_scans=good_prefilter_scans, uniq_prefilter_scans=uniq_prefilter_scans, $
                                  params_quality_metric=params_quality_metric, quality_metrics=quality_metrics

num_files          = N_ELEMENTS(prefscan_all_files)

IF N_ELEMENTS(verbose) EQ 0 THEN verbose=1

IF NOT KEYWORD_SET(params_quality_metric) THEN BEGIN
    param_quality_metric       = FLTARR(2,6)
    ; the range of acceptable average PMT counts around the prefilter peak transmission
    param_quality_metric(*,1)  = [1.0e4, 2.5e5]
    ; ratio of counts in the central portion of the scan versus the average of the two wings
    param_quality_metric(*,3)  = [1.01, 1.15]
    ; ratio of counts in the left versus the right wing
    param_quality_metric(*,4)  = [0.95, 1.05]
    ; the RMS values of the small-scale, normalized residuals of the prefilter scan profile
    param_quality_metric(*,5)  = [0.0, 0.025]
ENDIF

central_wave_range = [-0.125,0.125]
left_wave_range    = [-0.75,-0.50]
right_wave_range   = [-0.75,-0.50]
; the uniform wavelength step size onto which to interpolate the prefilter profiles
; -- in Angstroms
IF NOT KEYWORD_SET(wvstep_even) THEN wvstep_even         = 0.02

; Now scan through all the located prefilter scan log files and identify all the recorded prefilter scans.
; Some log files may contain duplicate informations (i.e. the log files may be copies or overlap in time),
; but we will sort that out later.
; During the scanning, we will gather the times and the lengths of all the prefilter scans (the number 
; of points in the scan may change with time due to changes settings in the acquisition program).
 
IF num_files GE 1 THEN BEGIN
    max_numscans = 0
    FOR logn=0,N_ELEMENTS(prefscan_all_files)-1 DO BEGIN 
        spawn,'grep ''Begin Prefilter Scan'' ' + STRING(prefscan_all_files(logn)) + ' | wc -l' ,output
        max_numscans += FIX(output)
    ENDFOR

    scan_lengths = FLTARR(max_numscans + 10)
    scan_dates   = DBLARR(max_numscans + 10)
    cnt=0                                                                                                                                                 

    FOR logn=0,N_ELEMENTS(prefscan_all_files)-1 DO BEGIN                                                                                                                  
        scan_log_info = parse_prefilter_log(prefscan_all_files(logn))
        FOR scann=0,N_ELEMENTS(scan_log_info)-1 DO BEGIN
            scan0 = load_prefilter_scan(scan_log_info(scann))
            scan_lengths(cnt) = N_ELEMENTS(scan0.COUNTS)
            scan_dates(cnt)   = scan0.scan_start_date
            cnt += 1
        ENDFOR
    ENDFOR
    IF verbose THEN PRINT,' --> Found ' + STRTRIM(cnt-1,2) + ' (non-unique) profile scans for ' + STRTRIM(scan0.filter_name,2) + ' filter.'
    
    num_samp = cnt - 1
    max_length = MAX(scan_lengths)  + 3
    midpos = FIX(max_length/2.)                                                                                                                                            

    scan_dates = scan_dates(0:cnt-1)
    IF verbose THEN BEGIN
        CALDAT,MIN(scan_dates),mo,dd,yr
        PRINT,'     Earliest Scan -  ' + STRING(yr,mo,dd,FORMAT='(%"%4i-%2.2i-%2.2i")')
        CALDAT,MAX(scan_dates),mo,dd,yr
        PRINT,'     Latest Scan   -  ' + STRING(yr,mo,dd,FORMAT='(%"%4i-%2.2i-%2.2i")')
    ENDIF
    
ENDIF

; Generate a structure that can hold any the prefilter scans, regardless of their length
prefilter_scan_all = CREATE_STRUCT('counts'      , FLTARR(max_length), $ 
                               'wavelengths'     , FLTARR(max_length), $
                               'fp1_voltages'    , FLTARR(max_length), $
                               'fp2_voltages'    , FLTARR(max_length), $
                               'valid'           , FLTARR(max_length), $
                               'scan_location'   , 0L, $
                               'scan_start_date' , 0.0d, $
                               'scan_end_date'   , 0.0d, $
                               'scan_date_text'  , '', $
                               'logfile_used'    , '', $
                               'filter_name'     , '' )
; And then make an array of those structures that is long enough for all the identified prefilter scans 
prefilter_scan_all = REPLICATE(prefilter_scan_all,num_samp+1)

cnt=0   

for logn=0,N_ELEMENTS(prefscan_all_files)-1 DO BEGIN                                                                                                        
    scan_log_info = parse_prefilter_log(prefscan_all_files(logn))

    for scann=0,N_ELEMENTS(scan_log_info)-1 DO BEGIN
    
        scan0          = load_prefilter_scan(scan_log_info(scann))
        scan0_len      = N_ELEMENTS(scan0.COUNTS)
        scan0_len_half = FIX(scan0_len / 2.)

        prefilter_scan_all(cnt).FP1_voltages(midpos - scan0_len_half:midpos + scan0_len_half) = scan0.FP1_voltages
        prefilter_scan_all(cnt).FP2_voltages(midpos - scan0_len_half:midpos + scan0_len_half) = scan0.FP2_voltages
        fp_volts_valid    = WHERE((ABS(scan0.FP1_voltages) LE 2048) AND (ABS(scan0.FP2_voltages) LE 2048))
        num_valid         = N_ELEMENTS(fp_volts_valid)
        points_filled     = INDGEN(scan0_len) + midpos - scan0_len_half
        
        prefilter_scan_all(cnt).Counts = (MEAN((scan0.Counts(fp_volts_valid))(0:4)) + MEAN((scan0.Counts(fp_volts_valid))(-5:-1))) / 2.
        prefilter_scan_all(cnt).Counts(points_filled(fp_volts_valid))       = scan0.Counts(fp_volts_valid)
        ;prefilter_scan_all(cnt).Counts(0:MIN(points_filled(fp_volts_valid)-1)                           = MEAN(scan0.Counts(0:4))
        ;prefilter_scan_all(cnt).Counts(midpos + scan0_len_half:*)                             = MEAN(scan0.Counts(-5:-1))

        prefilter_scan_all(cnt).Wavelengths(0:midpos-1)                           = -10
        prefilter_scan_all(cnt).Wavelengths(midpos+1:*)                           = 10
        prefilter_scan_all(cnt).Wavelengths(midpos - scan0_len_half:midpos + scan0_len_half)  = scan0.Wavelengths
 
        prefilter_scan_all(cnt).valid(points_filled(fp_volts_valid)) = 1
        
        prefilter_scan_all(cnt).Scan_Location   = scan0.Scan_Location
        prefilter_scan_all(cnt).Scan_Start_Date = scan0.Scan_Start_Date
        prefilter_scan_all(cnt).Scan_End_Date   = scan0.Scan_End_Date
        prefilter_scan_all(cnt).scan_date_text  = scan0.scan_date_text
        prefilter_scan_all(cnt).logfile_used    = scan0.logfile_used
        prefilter_scan_all(cnt).filter_name     = scan0.filter_name

        cnt += 1
        
    ENDFOR
ENDFOR

unique_scans        = UNIQ(prefilter_scan_all.Scan_Start_Date,SORT(prefilter_scan_all.Scan_Start_Date))
prefilter_scan_uniq = prefilter_scan_all(unique_scans)
num_uniq            = N_ELEMENTS(prefilter_scan_uniq)

IF verbose THEN PRINT,' --> Identified ' + STRTRIM(num_uniq,2) + ' unique profile scans for ' + STRTRIM(scan0.filter_name,2) + ' filter.'

quality_metrics = FLTARR(7,num_uniq)
quality_metrics(0,*) = REFORM(REBIN(prefilter_scan_uniq.Counts(midpos - 20:midpos + 20),1,num_uniq))

FOR scann=0,num_uniq-1 DO BEGIN
    mid0      = get_closest(prefilter_scan_uniq(scann).Wavelengths, central_wave_range)
    left0     = get_closest(prefilter_scan_uniq(scann).Wavelengths, left_wave_range)
    right0    = get_closest(prefilter_scan_uniq(scann).Wavelengths, right_wave_range)
    int_mid   = MEAN(prefilter_scan_uniq(scann).Counts(mid0[0]:mid0[1]))
    int_left  = MEAN(prefilter_scan_uniq(scann).Counts(left0[0]:left0[1]))
    int_right = MEAN(prefilter_scan_uniq(scann).Counts(right0[0]:right0[1]))
    
    ; mean counts in the central 0.25 Angstrom of the scan (depends on central_wave_range)
    quality_metrics(1,scann) = int_mid
    ; mean counts in the central 1.5 Angstrom of the scan
    quality_metrics(2,scann) = MEAN(prefilter_scan_uniq(scann).Counts(left0[0]:right0[1]))

    ; ratio of counts in the central portion of the scan versus the average of the two wings
    quality_metrics(3,scann) = int_mid / ((int_left + int_right)/2.)
    ; ratio of counts in the left versus the right wing
    quality_metrics(4,scann) = int_left / int_right
    
    ; the RMS values of the small scale residuals of the prefilter scan profile
    int_smooth = prefilter_scan_uniq(scann).Counts / SMOOTH(prefilter_scan_uniq(scann).Counts,15,/EDGE_TRUNC)
    quality_metrics(5,scann) = STDDEV(int_smooth(left0[0]:right0[1]))
ENDFOR

qcrit1 = (quality_metrics(1,*) LE param_quality_metric(1,1)) AND (quality_metrics(1,*) GE param_quality_metric(0,1))
qcrit3 = (quality_metrics(3,*) LE param_quality_metric(1,3)) AND (quality_metrics(3,*) GE param_quality_metric(0,3))  
qcrit4 = (quality_metrics(4,*) LE param_quality_metric(1,4)) AND (quality_metrics(4,*) GE param_quality_metric(0,4)) 
qcrit5 = (quality_metrics(5,*) LE param_quality_metric(1,5)) AND (quality_metrics(4,*) GE param_quality_metric(0,5)) 
quality_metrics(6,*) = qcrit1 * qcrit3 * qcrit4 * qcrit5
quality_valid        = WHERE(REFORM(quality_metrics(6,*)), valid_count)
IF verbose THEN PRINT,' --> Found ' + STRTRIM(valid_count,2) + ' profiles scans meeting the quality criteria for the ' + STRTRIM(scan0.filter_name,2) + ' filter.'

prefilter_scan_good           = prefilter_scan_uniq(quality_valid)
prefilter_scan_good_num       = N_ELEMENTS(prefilter_scan_good)

prefilter_scan_good_allwv     = prefilter_scan_good.wavelengths
prefilter_scan_good_wvmin     = MIN(prefilter_scan_good_allwv(WHERE(ABS(prefilter_scan_good_allwv) LE 9)), max=prefilter_scan_good_wvmax)
num_wv_steps                  = FIX((prefilter_scan_good_wvmax - prefilter_scan_good_wvmin) / wvstep_even / 2.) * 2 + 1
prefilter_scan_even           = FLTARR(num_wv_steps, prefilter_scan_good_num)
wavelength_scale_new          = (FINDGEN(num_wv_steps) - FIX(num_wv_steps/2.)) * wvstep_even

; Generate a structure that will hold all the prefilter scans after interpolation onto a common wavelength grid
prefilter_scan_even = CREATE_STRUCT('counts'      , FLTARR(num_wv_steps), $ 
                               'wavelengths'      , wavelength_scale_new, $
                               'fp1_voltages'     , FLTARR(num_wv_steps), $
                               'fp2_voltages'     , FLTARR(num_wv_steps), $
                               'valid'            , FLTARR(num_wv_steps), $
                               'scan_location'   , 0L, $
                               'scan_start_date' , 0.0d, $
                               'scan_end_date'   , 0.0d, $
                               'scan_date_text'  , '', $
                               'logfile_used'    , '', $
                               'filter_name'     , '' )
; And then make an array of those structures that is long enough for all the valid prefilter scans 
prefilter_scan_even = REPLICATE(prefilter_scan_even, prefilter_scan_good_num)

FOR scann=0,prefilter_scan_good_num-1 DO BEGIN
    valid_pts     = WHERE(prefilter_scan_good(scann).valid)
    wavematch_min = get_closest(wavelength_scale_new,MIN(prefilter_scan_good(scann).Wavelengths(valid_pts)))
    wavematch_max = get_closest(wavelength_scale_new,MAX(prefilter_scan_good(scann).Wavelengths(valid_pts)))
    
    prof_interp   = INTERPOL(prefilter_scan_good(scann).Counts(valid_pts), prefilter_scan_good(scann).Wavelengths(valid_pts), wavelength_scale_new(wavematch_min:wavematch_max),/SPLINE)
    
    prefilter_scan_even(scann).Counts(wavematch_min:wavematch_max) = prof_interp
    prefilter_scan_even(scann).valid(wavematch_min:wavematch_max)  = 1
    
    prefilter_scan_even(scann).Scan_Location   = prefilter_scan_good(scann).Scan_Location
    prefilter_scan_even(scann).Scan_Start_Date = prefilter_scan_good(scann).Scan_Start_Date
    prefilter_scan_even(scann).Scan_End_Date   = prefilter_scan_good(scann).Scan_End_Date
    prefilter_scan_even(scann).scan_date_text  = prefilter_scan_good(scann).scan_date_text
    prefilter_scan_even(scann).logfile_used    = prefilter_scan_good(scann).logfile_used
    prefilter_scan_even(scann).filter_name     = prefilter_scan_good(scann).filter_name

ENDFOR

good_prefilter_scans = prefilter_scan_good
uniq_prefilter_scans = prefilter_scan_uniq

RETURN,prefilter_scan_even

END
