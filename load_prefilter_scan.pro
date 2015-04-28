;-------------------------------------------------------------
;+
; NAME:
;       load_prefilter_scan
; PURPOSE:
;       
; CATEGORY:
; CALLING SEQUENCE:
;       prefilter_scan = load_profile_scan (prefilter_scan_location)
; INPUTS:
;       prefilter_scan_location = structure as output by parse_prefilter_scan_log.pro
; KEYWORD PARAMETERS:
;       verbose = turns on prining of each log line for the requested prefilter scan   
; OUTPUTS:
;       prefilter_scan = structure containing information on the selected prefilter scan
;                        see details on structure contents in NOTES section.
; COMMON BLOCKS:
;       
; NOTES:
;       The input to this function is a single element of the structure returned by
;       the parse_prefilter_scan_log.pro program.
;
;       e.g. IDL> prefilter_log_scans = parse_prefilter_scan_log (filter_name)
;            IDL> prefilter_scan      = load_profile_scan (prefilter_log_scans(i))
;
;       IDL> HELP, prefilter_scan, /STRUCTURE
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
;          Scan lengths or wavelength coverage may vary between scans, even those located in
;          the same log file, and this program does nothing to normalize these differences.
;
; MODIFICATION HISTORY:
;       19 Jul, 2003 - KPR - added documentation
;       26 Apr, 2015 - KPR - added additional documentation
;
;-
;-------------------------------------------------------------

FUNCTION load_prefilter_scan, prefilter_scan_location, verbose=verbose

IF N_ELEMENTS(verbose) LT 1 THEN verbose=0

scan_start               = 'Begin Prefilter Scan'
scan_end                 = 'End Prefilter Scan'
measurement_line         = 'Continuum Lamp'

input_record             = ''
current_line             = -1L
prefilter_counts         = 0.
wavelengths              = 0.
fp1_voltages             = 0.
fp2_voltages             = 0.

logfile                  = prefilter_scan_location.scan_log_file
scan_location            = prefilter_scan_location.scan_start_position

GET_LUN,log_lun
OPENR,log_lun,logfile
POINT_LUN, log_lun, scan_location

WHILE NOT (EOF(log_lun)) AND NOT (STRMATCH(input_record, scan_end + '*', /FOLD_CASE)) DO BEGIN
    READF,log_lun,input_record
    current_line = current_line + 1

    IF STRMATCH(input_record, measurement_line + '*', /FOLD_CASE) THEN BEGIN
        input_fields = STRSPLIT(input_record, /EXTRACT)
        v2start = STRPOS(input_fields(9),'=')
        v2end   = STRPOS(input_fields(9),')')
        v2      = STRMID(input_fields(9), v2start+1, v2end-v2start-1)
        
        prefilter_counts = [prefilter_counts, input_fields(2)]
        wavelengths      = [wavelengths, input_fields(5)]
        fp1_voltages     = [fp1_voltages, input_fields(8)]
        fp2_voltages     = [fp2_voltages, v2]

    ENDIF ELSE BEGIN
        IF verbose GE 1 THEN PRINT,input_record
    ENDELSE

ENDWHILE

FREE_LUN, log_lun

prefilter_scan = CREATE_STRUCT('counts'          , prefilter_counts(1:*), $ 
                               'wavelengths'     , wavelengths(1:*), $
                               'fp1_voltages'    , fp1_voltages(1:*), $
                               'fp2_voltages'    , fp2_voltages(1:*), $
                               'scan_location'   , scan_location, $
                               'scan_start_date' , prefilter_scan_location.scan_start_date, $
                               'scan_end_date'   , prefilter_scan_location.scan_end_date, $
                               'scan_date_text'  , prefilter_scan_location.scan_date_text, $
                               'logfile_used'    , logfile, $
                               'filter_name'     , prefilter_scan_location.filter_name)

RETURN, prefilter_scan

END ; load_prefilter_scan
