FUNCTION align_sum_profiles, flat_scan_ave, blueshift_map, wavelength_scale, $
                             plot_realtime=plot_realtime, radial_mask=radial_mask, $
			     line_start=line_start, line_end=line_end, $
			     interp_wavestep=interp_wavestep, aligned_profiles_slice=aligned_profiles_slice

size_ave_flats             = SIZE(flat_scan_ave)
xsize_ave_flats            = size_ave_flats(1)
ysize_ave_flats            = size_ave_flats(2)
zsize_ave_flats            = size_ave_flats(3)
output_text_length         = 1

IF NOT KEYWORD_SET(line_start) THEN line_start = 0
IF NOT KEYWORD_SET(line_end)   THEN line_end   = zsize_ave_flats - 1
IF NOT KEYWORD_SET(interp_wavestep)   THEN interp_wavestep   = 0.002

wavelengths_orig           = wavelength_scale(line_start:line_end)
blue_cor_range             = (MAX(blueshift_map) - MIN(blueshift_map))
line_wave_range            = MAX(wavelengths_orig) - MIN(wavelengths_orig) + blue_cor_range
line_wave_range_numpoint   = ROUND(line_wave_range / interp_wavestep) + 1
wavelengths_new            = FINDGEN(line_wave_range_numpoint) * interp_wavestep + (MIN(wavelengths_orig) - MAX(blueshift_map))

alignedprofile_flats       = DBLARR(line_wave_range_numpoint)
alignedprofile_flats_var   = DBLARR(line_wave_range_numpoint)
alignedprofile_flats_ext   = DBLARR(line_wave_range_numpoint)
alignedprofile_flats_cnt   = DBLARR(line_wave_range_numpoint)
aligned_profiles_slice     = FLTARR(MAX([xsize_ave_flats,ysize_ave_flats]),line_wave_range_numpoint,2) 

counter                    = 0L
mask_points_total          = N_ELEMENTS(WHERE(radial_mask EQ 1))
start_time                 = SYSTIME(1,/SECONDS)

;blueshift_map              = analytical_blueshift_map
;blueshift_map              = fitted_blueshift_map/1000.
;blueshift_map              = (blueshift_fit_7090_prefcor - $
;                             MIN(blueshift_fit_7090_prefcor(400:600,400:600))) * 0.002

PRINT,'Creating Average Profile - '

FOR xx=0L,xsize_ave_flats-1 DO BEGIN
    IF ((xx MOD 10) EQ 0) THEN BEGIN
        percent_finished   = FLOAT(counter)/mask_points_total*100
	seconds_elapsed    = SYSTIME(1,/SECONDS) - start_time
	seconds_remaining  = seconds_elapsed/(percent_finished>0.01) * (100 - percent_finished)
        IF seconds_remaining GE 61 THEN $
	    time_remaining = STRING(seconds_remaining/60,FORMAT='(F5.1)') + ' minutes' $
	ELSE $
	    time_remaining = STRING(seconds_remaining,FORMAT='(F4.1)') + ' seconds'
        ibis_overwrite_line, output_text_length, $
	    'Column # ' + STRING(xx, FORMAT='(I4.1)') + $
	    '  -- Percent Completed : ' + STRING(percent_finished, FORMAT='(F5.1)') + '%' + $
	    '  -- Time Remaining : ' + time_remaining, $
             new_overwrite_length=output_text_length
	IF KEYWORD_SET(plot_realtime) THEN BEGIN
            PLOT, wavelengths_new, alignedprofile_flats/(counter>1),/yn
            OPLOT, wavelengths_new, alignedprofile_flats_ext/(alignedprofile_flats_cnt>1),LI=1,th=3
	ENDIF
    ENDIF

    FOR yy=0,ysize_ave_flats-1 DO BEGIN

        IF (radial_mask(xx,yy) EQ 1) THEN BEGIN

            flat_single_spectrum      = REFORM(flat_scan_ave(xx,yy,line_start:line_end))
            wavelengths_new_and_blue  = wavelengths_new + blueshift_map(xx,yy)
            wavelengths_orig_and_blue = wavelengths_orig - blueshift_map(xx,yy)

            ; interpolate single flat field spectral line to common wavelength
            ; scale (i.e. no blueshift) and increase spectral sampling
            flat_spectrum_shifted     = INTERPOL(flat_single_spectrum, wavelengths_orig_and_blue, $
                                                 wavelengths_new, SPLINE=0, QUAD=0)
						 
	    validwv = WHERE((wavelengths_new GE MIN(wavelengths_orig_and_blue)) AND $
	                    (wavelengths_new LE MAX(wavelengths_orig_and_blue)))
	    alignedprofile_flats_ext(validwv)   = alignedprofile_flats_ext(validwv) + flat_spectrum_shifted(validwv)
	    alignedprofile_flats_cnt(validwv)   = alignedprofile_flats_cnt(validwv) + 1

            alignedprofile_flats      = alignedprofile_flats + flat_spectrum_shifted
            ; the following can be used to calculate RMS of averaged spectrum - is it useful?
            alignedprofile_flats_var  = alignedprofile_flats_var + (flat_spectrum_shifted^2)

            counter = counter + 1
	    
	    IF xx EQ ROUND(xsize_ave_flats/2.) THEN aligned_profiles_slice(yy,validwv,0) = flat_spectrum_shifted(validwv)
	    IF yy EQ ROUND(ysize_ave_flats/2.) THEN aligned_profiles_slice(xx,validwv,1) = flat_spectrum_shifted(validwv)
            
        ENDIF
    ENDFOR
ENDFOR

percent_finished   = FLOAT(counter)/mask_points_total*100
seconds_elapsed    = SYSTIME(1,/SECONDS) - start_time
seconds_remaining  = seconds_elapsed/(percent_finished>0.01) * (100 - percent_finished)
time_remaining     = STRING(seconds_remaining,FORMAT='(F4.1)') + ' seconds'
ibis_overwrite_line, output_text_length, $
    'Column # ' + STRING(xx, FORMAT='(I4.1)') + $
    '  - Percent Completed : ' + STRING(percent_finished, FORMAT='(F5.1)') + '%' + $
    '  - Time Remaining : ' + time_remaining, $
     new_overwrite_length=output_text_length
PRINT

alignedprofile_flats_ext = alignedprofile_flats_ext /  (alignedprofile_flats_cnt > 1)

alignedprofile_flats       = alignedprofile_flats / counter
; RMS variations of all average spectra - presently not put to any real use
alignedprofile_flats_var   = SQRT(((counter * alignedprofile_flats) - alignedprofile_flats_var) / $
                                              FLOAT(counter * (counter - 1)))

alignedprof = CREATE_STRUCT ('alignedprof',       alignedprofile_flats,$
                             'alignedprof_var',   alignedprofile_flats_var,$
                             'alignedprof_ext',   alignedprofile_flats_ext,$
                             'alignedprof_cnt',   alignedprofile_flats_cnt,$
			     'wave_scale_orig',   wavelengths_orig, $
			     'wave_scale_interp', wavelengths_new)

RETURN,alignedprof


END
