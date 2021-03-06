;-------------------------------------------------------------
;+
; NAME:
;       fits_date_convert
; PURPOSE:
;       Convert date text from FITS DATE-OBS format into numerical date
; CATEGORY:
; CALLING SEQUENCE:
;       converted_date = fits_date_convert(input_date)
; INPUTS:
;       input_date = data formatted according to FITS standard
; KEYWORD PARAMETERS:    
;       julian_day = output Julian date (default)
; OUTPUTS:
;       julian_date = DBLARR(1), generated by JULDAY
; COMMON BLOCKS:
;       
; NOTES:
;       The input date format is: YYYY-MM-DDTHH24:MI:SS.sss
;       The input date format is: YYYY-MM-DDZHH24:MI:SS.sss
;       The input date format is: YYYY-MM-DD HH24:MI:SS.sss
;       Example: 2003-06-20Z13:14:29.921
;       
; MODIFICATION HISTORY:
;       18 Nov, 2003 - KPR - initial implementation
;-
;-------------------------------------------------------------

FUNCTION fits_date_convert, fits_date

; default times if hour:min:sec components are not found
input_hour = 0
input_minute = 0
input_sec = 0

fits_date     = STRTRIM(fits_date, 2)
fits_date_len = N_ELEMENTS(fits_date)
julian_dates_all = DBLARR(fits_date_len)

sample_date   = fits_date(0)
IF STRPOS(sample_date, 'Z') GE 1 THEN separator='Z' ELSE $
IF STRPOS(sample_date, 'T') GE 1 THEN separator='T' ELSE $
                                    separator=' '

FOR daten = 0LL,fits_date_len - 1 DO BEGIN
    fits_date_sel   = fits_date(daten)
    IF STRMATCH(fits_date_sel, '*Not Available*') EQ 0 THEN BEGIN
    split_date_time = STRSPLIT(fits_date_sel, separator, /EXTRACT, COUNT=split_date_num)
    
    date_only = split_date_time(0)
    time_only = split_date_time(1)
    
    split_date      = STRSPLIT(date_only, '-', /EXTRACT)
    input_year      = split_date(0)
    IF (STRLEN(input_year) EQ 2) THEN input_year = 1900 + FIX(input_year) ELSE $
                                      input_year = FIX(input_year) 
    input_month     = FIX(split_date(1))
    input_day       = FIX(split_date(2))
    
    IF (split_date_num GE 2) THEN BEGIN
        time_only       = split_date_time(1)
        split_time      = STRSPLIT(time_only, ':', /EXTRACT, COUNT=split_time_num)
        input_hour      = FIX(split_time(0))
        IF (split_time_num GE 2) THEN input_minute    = FIX(split_time(1))
        IF (split_time_num GE 3) THEN input_sec       = FLOAT(split_time(2))
    ENDIF
    
    julian_date          = JULDAY(input_month, input_day, input_year, input_hour, input_minute, input_sec)

    julian_dates_all(daten) = julian_date
    ENDIF
    
ENDFOR

julian_dates_all = REFORM(julian_dates_all)
IF N_ELEMENTS(julian_dates_all) EQ 1 THEN julian_dates_all = julian_dates_all[0]

RETURN, julian_dates_all

END                  ; END fits_date_convert
