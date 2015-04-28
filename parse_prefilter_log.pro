;-------------------------------------------------------------
;+
; NAME:
;       parse_prefilter_log
; PURPOSE:
;       reads a log of prefilter scans to determine the location (beginning and end)
;       of the different prefilter scans stored in the file
; CATEGORY:
; CALLING SEQUENCE:
;       prefilter_log_scans = parse_prefilter_log (filter_name)
; INPUTS:
;       logfile = the name of the IBIS prefilter scan log file that is to be read
;                 these generally are of the form "sintodata.${wavelemngth}.dat"
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       pref_log_scans = array of structures (see NOTES), one-entry per 
;                         prefilter scan in the input log file.
; OPTIONAL OUTPUTS:
;       imaging_spect_scans = array of structures containing information on
;                             the imaging spectral scan times as recorded in the same 
;                             input log file.
; COMMON BLOCKS:
;
; NOTES:
;
;       The prefilter scan log typically records multiple prefilter scans in a single
;       text file. The scan information (wavelength, PMT counts) is recorded as a 
;       series of individual lines for each wavelength tuning position of the scan.
;       This program scans through the log file and identifies the beginning and ending
;       of each individual prefilter scan recorded in the supplied log file. This way
;       it is easy to go back (using load_prefilter_scan.pro) and load individual 
;       prefilter scans as needed.
;
;       One element (structure) of the output array of structures from this 
;       function is generally passed to the load_profile_scan.pro program 
;       to extract the information on a selected prefilter scan.
;
;       e.g. IDL> prefilter_log_scans = parse_prefilter_log (filter_name)
;            IDL> prefilter_scan      = load_profile_scan (prefilter_log_scans(i))
;
;       The format of the output structure is as follows:
;       IDL> HELP,prefilter_log_scans,/STRUCTURE
;       ** Structure <70ad90>, 10 tags, length=80, data length=76, refs=1:
;          SCAN_START_DATE     - starting date of prefilter scan as a Julian date
;          SCAN_START_LINE     - starting line number for prefilter scan
;          SCAN_START_POSITION - starting position (from POINT_LUN) for prefilter scan
;          SCAN_END_DATE       - ending date of prefilter scan as a Julian date
;          SCAN_END_LINE       - ending line number for prefilter scan
;          SCAN_END_POSITION   - ending position (from POINT_LUN) for prefilter scan
;          SCAN_END_START_DATE - starting date of prefilter scan, as recorded at the end of
;                                scan, given as a Julian date; should be identical
;                                to SCAN_START_DATE
;          SCAN_DATE_TEXT      - the starting and ending dates given in original format
;          SCAN_LOG_FILE       - name of the log file to which the parsing
;          FILTER_NAME         - input name of filter
;
; MODIFICATION HISTORY:
;       19 Jul, 2003 - KPR - added documentation
;       27 Apr, 2015 - KPR - renamed/copied from parse_prefilter_scan_log.pro 
;                          - reworked to simply accept a full log file name as the input
;                            rather than constructing the logfile name and searching in
;                            a fixed directory. 
;-
;-------------------------------------------------------------

FUNCTION parse_prefilter_log, logfile, imaging_spect_scans=imaging_spect_scans

FORWARD_FUNCTION FILE_SEARCH

scan_start_text      = 'Begin Prefilter Scan'
scan_end_text        = 'End Prefilter Scan'
scan_imaging_record  = 'Imaging Spectral Scan'
input_record         = ''
current_line         = -1L

pref_log_scans       = CREATE_STRUCT('scan_start_date'               , 0.0D, $
                                     'scan_start_line'     , 0L  , $
                                     'scan_start_position' , 0L  , $
                                     'scan_end_date'       , 0.0D, $
                                     'scan_end_line'       , 0L  , $
                                     'scan_end_position'   , 0L  , $
                                     'scan_end_start_date' , 0.0D, $
                                     'scan_date_text'      , ''  , $
                                     'scan_log_file'       , ''  , $
                                     'filter_name'         , ''    $
                                     )
                                     
imaging_spect_scans     = CREATE_STRUCT('scan_date' ,        0.0d, $
                                        'scan_date_text'     , '', $
                                        'scan_directory'     , '', $
                                        'center_voltage_fp1' , 0, $
                                        'center_voltage_fp2' , 0, $
                                        'filter_name'        , ''  $
                                        )
latest_imaging_scan = imaging_spect_scans



logfile_is_valid     = FILE_TEST(logfile, /READ)

IF (NOT logfile_is_valid) THEN BEGIN

    PRINT, 'Error! Could not locate logfile: ' + STRTRIM(logfile, 2)

ENDIF ELSE BEGIN

    PRINT, 'Using logfile: ' + STRTRIM(logfile, 2)
    
    ; try to identify which filter the logfile pertains to based on 
    ; a simple pattern search in the filename.
    ; All IBIS filters should be between 5000-9000 Å
    IF STRMATCH(logfile,'*[1-9][0-9][0-9][0-9]*') THEN BEGIN
        filter_name = STRMID(logfile, STREGEX(logfile,'[5-8][0-9][0-9][0-9]'), 4)
    ENDIF ELSE IF STRMATCH(logfile,'*Dark*') THEN BEGIN
        filter_name = 'dark'
    ENDIF ELSE IF STRMATCH(logfile,'*Open*') THEN BEGIN
        filter_name = 'open'
    ENDIF ELSE BEGIN
        filter_name = 'unknown'
    ENDELSE

    GET_LUN,log_lun
    OPENR,log_lun,logfile

    WHILE NOT (EOF(log_lun)) DO BEGIN
        POINT_LUN, -log_lun, current_position
        READF,log_lun,input_record
        current_line = current_line + 1

        ; see if we are at the beginning of a prefilter scan
        IF STRMATCH(input_record, scan_start_text + '*', /FOLD_CASE) THEN BEGIN
            input_fields         = STRSPLIT(input_record, '-', /EXTRACT)
            scan_start_date      = ibis_log_date_convert(input_fields(3))
            scan_start_line      = current_line
            scan_start_position  = current_position

        ; check to see if we are at the end of prefilter scan, and if so
        ; construct a structure summarizing the scan
        ENDIF ELSE IF STRMATCH(input_record, scan_end_text + '*', /FOLD_CASE) THEN BEGIN
            POINT_LUN, -log_lun, current_position
            input_fields         = STRSPLIT(input_record, '-', /EXTRACT)
            scan_end_start_date  = ibis_log_date_convert(input_fields(3))
            scan_end_date        = ibis_log_date_convert(input_fields(5))
            scan_date_text       = 'Prefilter Scan Date: ' + input_fields(3) + ' - ' + input_fields(5)

            pref_log_scans_temp  = CREATE_STRUCT('scan_start_date' , scan_start_date, $
                                            'scan_start_line'      , scan_start_line, $
                                            'scan_start_position'  , LONG(scan_start_position), $
                                            'scan_end_date'        , scan_end_date, $
                                            'scan_end_line'        , current_line, $
                                            'scan_end_position'    , LONG(current_position), $
                                            'scan_end_start_date'  , scan_end_start_date, $
                                            'scan_date_text'       , scan_date_text, $
                                            'scan_log_file'        , logfile, $
                                            'filter_name'          , STRTRIM(filter_name,2))

            ; append new prefilter scan record to existing array of structures
            pref_log_scans       = [pref_log_scans, pref_log_scans_temp]

        ; see if the line is a record of an imaging spectral scan (i.e. with the CCD)
        ; if so, parse the line to extract information on the directory location, 
        ; time, and FPI voltages
        ENDIF ELSE IF STRMATCH(input_record, scan_imaging_record + '*', /FOLD_CASE) THEN BEGIN
            input_fields    = STRSPLIT(input_record,'--',/EXTRACT,/REGEX)
            field1          = STRSPLIT(input_fields(1), /EXTRACT)
            filter_name     = field1(1)
            fp1start        = STRPOS(input_fields(1), '=')
            fp1end          = STRPOS(input_fields(1), 'V',fp1start)
            fp1_voltage     = FIX(STRMID(input_fields(1), fp1start+1, fp1end-fp1start-1))
            fp2start        = STRPOS(input_fields(1), '=', fp1end)
            fp2end          = STRPOS(input_fields(1), ')',fp2start)
            fp2_voltage     = FIX(STRMID(input_fields(1), fp2start+1, fp2end-fp2start-1))
            date_jd         = ibis_log_date_convert(input_fields(3))
            image_directory = STRTRIM(input_fields(5),2)

            ; populate an empty structure with the parsed out information from the log record
            imaging_spect_scan_temp  = CREATE_STRUCT('scan_date' ,           date_jd, $
                                                     'scan_date_text'      , input_fields(3), $
                                                     'scan_directory'      , image_directory, $
                                                     'center_voltage_fp1'  , fp1_voltage, $
                                                     'center_voltage_fp2'  , fp2_voltage, $
                                                     'filter_name'         , filter_name  $
                                                     )
 
            ; append new imaging spectral scan record to existing array of structures
            imaging_spect_scans  = [imaging_spect_scans, imaging_spect_scan_temp]
            latest_imaging_scan  = imaging_spect_scan_temp
        ENDIF

    ENDWHILE

    FREE_LUN, log_lun
        
    ; if any prefilter scans were found in the file, then we will discard the 
    ; initial (zero-index) "dummy" element of the array of prefilter scan 
    ; structures. Otherwise, we will return the empty structure we initially created.
    nel_pref_log_scans   = N_ELEMENTS(pref_log_scans)
    IF nel_pref_log_scans GT 1 THEN $
        pref_log_scans = pref_log_scans(1:*)

    ; same thing for the imaging spectral scan structures
    nel_imaging_spect_scans = N_ELEMENTS(imaging_spect_scans)
    IF nel_imaging_spect_scans(0) GT 1 THEN $
        imaging_spect_scans = imaging_spect_scans(1:*)

ENDELSE

RETURN, pref_log_scans

END
