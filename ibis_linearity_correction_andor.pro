;-------------------------------------------------------------
;+
; NAME:
;       ibis_linearity_correction_andor
; PURPOSE:
;       to correct an IBIS image for non-linearites in the
;       system response
; CATEGORY:
; CALLING SEQUENCE:
;       ibis_linearity_correction_andor, raw_data
; INPUTS:
;       data_raw = the raw image as read in, either with or without
;                  the bias level removed
; KEYWORD PARAMETERS:
;       camera_id = serial number of Andor camera that acquired
;                   the input image
;       data_includes_bias = indicates whether input image has not
;                   had its bias level removed
;       data_bias_level = the mean bias level for the input data;
;                    equivalent to 0 for data with bias level removed
;       
; OUTPUTS:
;       file_selected = returns a string containing the name of the linearity
;                       that was selected and used for the linearity correction
; COMMON BLOCKS:
;       ibis_config
; NOTES:
;       The linearity curve calculated without bias can be converted to a 
;       correction curve for data with an arbitrary bias level according to
;       the following formula:
;       linearity_curve_withbias = (bias+COUNTS_MEASURED_RAW)/(bias+COUNTS_MEASURED_RAW/CORRECTION_LINEARITY)
;       where the relative measured counts become:
;       COUNTS_MEASURED_RAW_withbias = (bias+COUNTS_MEASURED_RAW) 
; MODIFICATION HISTORY:
;       K. Reardon, Feb 2011, adapted ibis_linearity_correction
;                             for use with new Andor cameras
;                             (and their messy behavior, including
;                             pixel-to-pixel variations in the linearity,
;                             which we try to compensate for, but not (yet)
;                             perform a full correction.
;-
;-------------------------------------------------------------
FUNCTION ibis_linearity_correction_andor, data_raw, camera_id=camera_id, $
         data_date=data_date, data_includes_bias=data_includes_bias, $
         data_bias_level=data_bias_level, $
         verbose=verbose, file_selected=file_selected,$
         exact_interpolate=exact_interpolate, caldir=caldir, $
         lin_ptr=lin_ptr

    ;COMMON ibis_config, ibis_environment

    ; some hardwired values that could eventually become input parameters if needed
    
    ; the maximum fraction of the total pixels to replace at high or low end of
    ; signal range. If more than this fraction are outside the range of acceptable
    ; values, then no correction at all will be applied.
    cutoff  = 0.05
    ; the width of the neighborhood over which to perform the median smoothing
    ; to determine the replacement value for the out-of-range pixels
    med_wid = 5
    ; replace the out-of-range pixels with the local median
    ; the alternative, for use_median=0, is to replace the pixels with the local average
    use_median = 1
    valid_count = 0
    do_linearity_cor = 1
    skip_low_cor = 1

    IF N_ELEMENTS(verbose) LT 1 THEN verbose=0
    IF NOT KEYWORD_SET(data_date) THEN data_date=SYSTIME(/JULIAN)
    ; if lin_ptr is not a single element (i.e. 1 pointer) then set it to zero to
    ; avoid problems with PTR_VALID later. 
    IF N_ELEMENTS(lin_ptr) GT 1 THEN lin_ptr=0
    ; if the lin_ptr keyword was passed in and is either a) undefined or b) a valid pointer
    ; (though not necessarily valid inearity correction data), 
    ; then we will keep the pointer and pass it back to the user. 
    ; Otherwise assume the user has no use for the pointer and free it at the end of execution.
    IF ( ARG_PRESENT(lin_ptr) EQ 1 AND N_ELEMENTS(lin_ptr) EQ 0) OR $
       ( ARG_PRESENT(lin_ptr) EQ 1 AND N_ELEMENTS(lin_ptr) EQ 1 AND PTR_VALID(lin_ptr)) $
       THEN keep_pointer=1 ELSE keep_pointer=0

    ; if input data_date is a string then try to convert it to a julian date
    ; (string is assumed to be in the format YYYY-DD-MMZHH24:MI:SS.sss
    ; if data_date is not a string, blindly assume it is already a julian date!
    data_date_type = SIZE(data_date, /TNAME)
    IF (data_date_type EQ 'STRING') THEN $    
        data_date_jd = fits_date_convert(data_date) $
    ELSE $
        data_date_jd = data_date
    data_date_jd = data_date_jd[0]
    
    ; by default we will load the linearity file, but we'll check to see 
    ; a pointer to the lineary corrected was passed in and if it is 
    ; valid for the data being corrected
    load_linearity = 1
    IF PTR_VALID(lin_ptr) THEN BEGIN
        IF N_ELEMENTS((*lin_ptr).camera_id) GE 1 THEN BEGIN
            IF (*lin_ptr).camera_id NE camera_id THEN BEGIN 
                PRINT,'Error: Data Camera ID (' + STRTRIM(camera_id,2) + $
                      ') doesn''t match input linearity data (' + $
                      STRTRIM((*lin_ptr).camera_id,2) + ').'
                 load_linearity = 1
             ENDIF ELSE BEGIN
                 load_linearity = 0
             ENDELSE
        ENDIF
    ENDIF 
                                                             
    IF load_linearity THEN BEGIN
        IF NOT(KEYWORD_SET(caldir)) THEN BEGIN
            this_func_info = ROUTINE_INFO('ibis_linearity_correction_andor',/SOURCE,/FUNCTION)
            this_func_path = FILE_DIRNAME(this_func_info.path, /MARK)
            caldir         = this_func_path
        ENDIF
    
        IF FILE_TEST(caldir, /DIRECTORY, /READ) THEN BEGIN
            linearity_files_andor = FILE_SEARCH(caldir +  'linearity.correction.*andor*sav', count=andor_linearity_count)
    
            IF (andor_linearity_count GE 1) THEN BEGIN
                date_string    = STRARR(andor_linearity_count)
                date_jd        = DBLARR(andor_linearity_count)
                camera_id_file = STRARR(andor_linearity_count)
                bias_included  = BYTARR(andor_linearity_count)
                
                FOR filenum=0,andor_linearity_count - 1 DO BEGIN
                    filename_split       = STRSPLIT(file_basename(linearity_files_andor(filenum)), '.', /EXTRACT)
                    date_string(filenum) = filename_split(2)
                    date_year            = STRMID(date_string(filenum), 0, 4)
                    date_month           = STRMID(date_string(filenum), 4, 2)
                    date_day             = STRMID(date_string(filenum), 6, 2)
                    date_jd(filenum)     = JULDAY(date_month, date_day, date_year, 0, 0, 0)
                    
                    camera_id_file(filenum) = STRMID(filename_split[4],1,4)
                    bias_included(filenum) = 0
                ENDFOR
    
                andor_linearity_file_info  = CREATE_STRUCT('linearity_count', andor_linearity_count, $
                                                   'filename'       , linearity_files_andor, $
                                                   'date_string'    , date_string, $
                                                   'date_jd'        , date_jd, $
                                                   'camera_id'      , camera_id_file, $
                                                   'bias_included'  , bias_included $
                                                   )
    
                ; identify linearity curves matching the type of data input
                ; i.e. with or without the bias level removed
                camera_linearity_curves = WHERE(andor_linearity_file_info.camera_id EQ camera_id, valid_count)
                ;HELP,andor_linearity_file_info,/ST
    
             ENDIF ELSE BEGIN
                 do_linearity_cor = 0
                 IF (verbose GE 1) THEN PRINT,'No linearity correction files found!'
             ENDELSE
                     
        ENDIF ELSE BEGIN
            do_linearity_cor = 0
            IF (verbose GE 1) THEN PRINT,'Invalid calibration directory -- ' + STRTRIM(caldir,2)
        ENDELSE
    
        IF (valid_count GE 1) THEN BEGIN
            selected_curves = camera_linearity_curves
        
            ;PRINT,andor_linearity_file_info.date_string(selected_curves)
            ;HLP,valid_linearity_curves,selected_curves
        
            ; find closest linearity scan to the input 
            linearity_dates_jd = andor_linearity_file_info.date_jd
            time_difference = MIN(ABS(linearity_dates_jd(selected_curves) - data_date_jd), closest_linearity)
            linearity_file  = andor_linearity_file_info.filename(selected_curves(closest_linearity))
        
            file_selected = linearity_file
            IF (verbose GE 1) THEN PRINT,'Using File - ' + linearity_file
            RESTORE,linearity_file

            lin_ptr = PTR_NEW(lincor)
            
        ENDIF ELSE do_linearity_cor = 0
    ENDIF ; load_linearity

    IF (do_linearity_cor GE 1) THEN BEGIN

        IF N_ELEMENTS(data_bias_level) LT 1 THEN data_bias_level = (*lin_ptr).bias_level
        
        data_raw_adj = data_raw + ((*lin_ptr).bias_level - data_bias_level)
        
        IF N_ELEMENTS((*lin_ptr).lowresp_lincor) GT 1 THEN $
            data_raw_adj = data_raw_adj / (*lin_ptr).lowresp_lincor
    
        IF (N_ELEMENTS((*lin_ptr).nonlin_low_limit) GT 1) OR (N_ELEMENTS((*lin_ptr).nonlin_high_limit) GT 1)THEN BEGIN

            ; Find points in image that lie out of the range of linearity,
            ; the value of which varies for each pixel. Need to check and correct 
            ; both high and low pixels (though having both in the same image is improbable).
            over_mask      =  (data_raw_adj - (*lin_ptr).nonlin_high_limit) GT 0
            over_mask_whr  = WHERE(over_mask, over_mask_num)
            under_mask     = (data_raw_adj>0 - (*lin_ptr).nonlin_low_limit)  LT 0
            under_mask_whr = WHERE(under_mask, under_mask_num)
            data_raw_cnt   = FLOAT(N_ELEMENTS(data_raw_adj))
            max_data_cor   = data_raw_cnt*cutoff
            IF skip_low_cor EQ 1 THEN under_mask_num = 0
            
            ; If out-of-range points are found, a median filter is applied to the images
            ; (which is computationally expensive) and the out of bound points are
            ; substituted with the median value of the surrounding area. Alternatively,
            ; the average of the surrounding area can be used, which is faster, but gives
            ; a worse (i.e. noisier) result.
            ; If there are too many bad points (e.g. more than 5% of the total), 
            ; we bail and don't apply any of this correction (and print out an error
            ; message).
            IF (over_mask_num GE 1) OR (under_mask_num GE 1) THEN BEGIN

                IF use_median THEN $
                     data_raw_adj_cor = MEDIAN(data_raw_adj,med_wid,/EVEN) $
                ELSE data_raw_adj_cor = (SMOOTH(data_raw_adj,med_wid,/EDGE) * (med_wid^2) - data_raw_adj)/(med_wid^2 -1)

                IF (over_mask_num GE 1) AND (over_mask_num LT max_data_cor) THEN BEGIN
                    data_raw_adj(over_mask_whr)  = data_raw_adj_cor(over_mask_whr)
                    IF (verbose GE 2) THEN $
                        PRINT,'ibis_linearity_correction_andor: Correcting ' + $
                               STRTRIM(over_mask_num) + ' partially saturated pixels.'
                ENDIF

                ;IF (under_mask_num GE 1) AND (under_mask_num LT max_data_cor) THEN BEGIN
                IF (under_mask_num GE 1) THEN BEGIN
                    data_raw_adj(under_mask_whr) = data_raw_adj_cor(under_mask_whr)
                    IF (verbose GE 2) THEN $
                        PRINT,'ibis_linearity_correction_andor: Correcting ' + $
                               STRTRIM(under_mask_num) + ' low intensity pixels.'
                ENDIF
            ENDIF
            
;            IF     (over_mask_num GE max_data_cor) OR (under_mask_num GT max_data_cor) $
            IF     (over_mask_num GE max_data_cor) $
               AND (verbose GE 1) THEN BEGIN
               PRINT,'ibis_linearity_correction_andor: Excessive number of out-of-range points found!'
               PRINT,'ibis_linearity_correction_andor: Full non-linearity correction not applied'
            ENDIF

        ENDIF
    
        ; The index values in correction_linearity correspond to the raw data
        ; ADU (e.g. index value 128 might be the correction for ADU values
        ; between 1024 and 1032), so we can use the data values directly to 
        ; the index value from the correction curve.
        ; We use the raw data values, scaled to the range of the correction_curve
        ; indices to construct an array of indices 
        ; corresponding to index value in correction_linearity For example, 
        ; if there were 
        data_cor = FLOAT(data_raw_adj) / $
                   (*lin_ptr).nonlin_fit_crv( FIX(data_raw_adj/(*lin_ptr).nonlin_fit_step) )
    
        data_cor = data_cor - ((*lin_ptr).BIAS_LEVEL - data_bias_level)
    
        data_cor = FIX(ROUND(data_cor))
        
    ENDIF ELSE BEGIN
        data_cor = data_raw
        IF (verbose GE 1) THEN PRINT,'No linearity correction applied!'
    ENDELSE

IF (keep_pointer NE 1) THEN PTR_FREE, lin_ptr

RETURN,data_cor

END
