PRO extract_dst_header_info, input_files, output_save_file=output_save_file

num_files     = N_ELEMENTS(input_files)

; test to see how many extensions there are in the first file
; DANGER! we will assume all the files have the same number of extensions.
fits_open,input_files[0],fitsblock
extn_total = fitsblock.nextend
fits_close,fitsblock

; generate the arrays into which to store the data
date_obs_str        = STRARR(extn_total, num_files)
date_end_str        = STRARR(extn_total, num_files)
dst_times_str       = STRARR(extn_total, num_files)
dst_seeing_str      = STRARR(extn_total, num_files)
dst_lightlevel_str  = STRARR(extn_total, num_files)
dst_elevation_str   = STRARR(extn_total, num_files)
dst_azimuth_str     = STRARR(extn_total, num_files)
dst_tableangle_str  = STRARR(extn_total, num_files)
dst_guiderangle_str = STRARR(extn_total, num_files)

for filn=0,num_files-1 do begin
    fits_open,input_files[filn],fitsblock
    ;IF (filn MOD 20) EQ 0 THEN PRINT,SYSTIME(),filn
    for extn=0,extn_total-1 do begin
	fits_read,fitsblock,exten=extn+1,data,header,/header_only
	date_obs_str[extn,filn]        = sxpar(header,'DATE-OBS',/Silent)
	date_end_str[extn,filn]        = sxpar(header,'DATE-OBS',/Silent)
	dst_times_str[extn,filn]       = sxpar(header,'DST_TIME',/Silent)
	dst_seeing_str[extn,filn]      = sxpar(header,'DST_SEE',/Silent)
	dst_lightlevel_str[extn,filn]  = sxpar(header,'DST_LLVL',/Silent)
	dst_elevation_str[extn,filn]   = sxpar(header,'DST_EL',/Silent)
	dst_azimuth_str[extn,filn]     = sxpar(header,'DST_AZ',/Silent)
	dst_tableangle_str[extn,filn]  = sxpar(header,'DST_TBL',/Silent)
	dst_guiderangle_str[extn,filn] = sxpar(header,'DST_GDRN',/Silent)
    endfor
    fits_close,fitsblock
endfor
; convert the FITS date strings to Julian dates
date_obs_jd           = fits_date_convert(date_obs_str)
date_end_jd           = fits_date_convert(date_end_str)
dst_times_jd_all      = fits_date_convert(dst_times_str)
; find unique samples of DST telescope parameters from the 
;     unique values of the timestamps of the parameters
;     DST parameters were sampled every 10.05 seconds
dst_times_uniq_idx    = UNIQ(dst_times_jd_all)
; now select only the unique values of each parameter, based 
;     on the unique values of the DST_TIME.
dst_times_uniq        = dst_times_jd_all[dst_times_uniq_idx]
dst_seeing_uniq       = FLOAT(dst_seeing_str[dst_times_uniq_idx])
dst_lightlevel_uniq   = FLOAT(dst_lightlevel_str[dst_times_uniq_idx])
dst_elevation_uniq    = FLOAT(dst_elevation_str[dst_times_uniq_idx])
dst_azimuth_uniq      = FLOAT(dst_azimuth_str[dst_times_uniq_idx])
dst_tableangle_uniq   = FLOAT(dst_tableangle_str[dst_times_uniq_idx])
dst_guiderangle_uniq  = FLOAT(dst_guiderangle_str[dst_times_uniq_idx])
; interpolate the unique seeing and light level values onto the same 
;    time grid as the data themselves
dst_seeing_interp     = INTERPOL(dst_seeing_uniq, dst_times_uniq, date_obs_jd)
dst_lightlevel_interp = INTERPOL(dst_lightlevel_uniq, dst_times_uniq, date_obs_jd)

IF keyword_set(output_save_file) THEN BEGIN
    IF TYPENAME(output_save_file) NE 'STRING' THEN output_save_file = 'FITS.DST.header.info.parsed.sav'
    SAVE,date_obs_jd, date_end_jd, dst_times_jd, dst_times_jd_all, dst_seeing_uniq, $
	dst_lightlevel_uniq, dst_azimuth_uniq, dst_elevation_uniq, dst_tableangle_uniq, $
	dst_guiderangle_uniq, dst_seeing_interp, dst_lightlevel_interp, filename=output_save_file
ENDIF

END
