FUNCTION construct_shifted_profile_array, aligned_profile, $
             blueshift_map, wavelength_scale_orig, wavelength_scale_interp, $
	     radial_mask=radial_mask

zsize_scan                  = N_ELEMENTS(wavelength_scale_orig)

size_image                  = SIZE(blueshift_map)
xsize_image          	    = size_image(1)
ysize_image	            = size_image(2)

shiftedprofile_flats_interp = FLTARR(zsize_scan)
; the array where the shifted mean profiles will be stored
shiftedprofile_flats        = FLTARR(xsize_image, ysize_image, zsize_scan)

;blueshift_map              = analytical_blueshift_map
;blueshift_map              = fitted_blueshift_map/1000.
;blueshift_map               = (blueshift_fit_7090_prefcor - $
;                             MIN(blueshift_fit_7090_prefcor(400:600,400:600))) * 0.002

output_text_length          = 1
start_time                  = SYSTIME(1,/SECONDS)
counter                     = 0L
mask_points_total           = N_ELEMENTS(WHERE(radial_mask EQ 1))
PRINT,'Creating Shifted Mean Spectrum Array - '

FOR xx=0,xsize_image-1 DO BEGIN
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
    ENDIF	
	
    FOR yy=0,ysize_image-1 DO BEGIN
    
        IF (radial_mask(xx,yy) EQ 1) THEN BEGIN

            ; construct wavelength scale at this (xx,yy) position including blueshift
	    wavelengths_orig_with_blue  = wavelength_scale_orig - blueshift_map(xx,yy)

            ; interpolate mean flat field spectral line to shifted wavelength
	    ; range and with same spectral sampling as original data
	    shiftedprofile_flats_interp = INTERPOL(aligned_profile, wavelength_scale_interp, $
	                                           wavelengths_orig_with_blue, /SPLINE)
	
	    shiftedprofile_flats(xx,yy,*) = shiftedprofile_flats_interp
	
	    counter = counter + 1
	ENDIF
    ENDFOR
ENDFOR

percent_finished   = FLOAT(counter)/mask_points_total*100
seconds_elapsed    = SYSTIME(1,/SECONDS) - start_time
seconds_remaining  = seconds_elapsed/(percent_finished>0.01) * (100 - percent_finished)
; since we are done, we know the result should be given in seconds
time_remaining = STRING(seconds_remaining,FORMAT='(F4.1)') + ' seconds'

ibis_overwrite_line,output_text_length,$
               'Column # ' + STRING(xx, FORMAT='(I4.1)') + $
               '  -- Percent Completed : ' + STRING(percent_finished, FORMAT='(F5.1)') + '%' + $
	       '  -- Time Remaining : ' + time_remaining, $
                new_overwrite_length=output_text_length

PRINT

RETURN,shiftedprofile_flats

END
