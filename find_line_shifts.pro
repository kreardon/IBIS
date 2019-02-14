FUNCTION find_line_shifts, flat_scan_ave, wavelength_scale_orig,   $
			   prefilter_correction=prefilter_correction, $
			   wavelength_scale_interp = wavelength_scale_interp, $
			   line_start=line_start, line_end=line_end,  $
			   interp_wavestep=interp_wavestep, $
                           fit_points=fit_points, radial_mask=radial_mask
        

size_ave_flats             = SIZE(flat_scan_ave)
xsize_ave_flats            = size_ave_flats(1)
ysize_ave_flats            = size_ave_flats(2)
zsize_ave_flats            = size_ave_flats(3)

IF NOT KEYWORD_SET(interp_wavestep) THEN interp_wavestep = 0.002       ; Angstroms
; fit_points = number of pixels to fit to line center
; spectral fit range = fit_points * interp_wavestep
; (?) does this need to change with lines of different size/depth?
IF NOT KEYWORD_SET(fit_points)      THEN fit_points      = 21

IF NOT KEYWORD_SET(line_start)      THEN line_start      = 0
IF NOT KEYWORD_SET(line_end)        THEN line_end        = zsize_ave_flats-1
IF NOT KEYWORD_SET(prefilter_correction) THEN $
      prefilter_correction = FLTARR(line_end-line_start+1) + 1.0

line_wave_range            = MAX(wavelength_scale_orig(line_start:line_end),MIN=minwave) - minwave
line_wave_range_numpoint   = line_wave_range / interp_wavestep
wavelengths_new            = FINDGEN(line_wave_range_numpoint) * interp_wavestep + wavelength_scale_orig(line_start)
line_centers               = FLTARR(xsize_ave_flats, ysize_ave_flats, 2)
wavelength_scale_interp    = wavelengths_new

output_text_length         = 1
mask_points_total          = N_ELEMENTS(WHERE(radial_mask EQ 1))
counter                    = 0L
start_time                 = SYSTIME(1,/SECONDS)

PRINT,'Finding Line Center Positions - '
; -----  Calculate Line Position  -----
FOR xx=0, xsize_ave_flats-1 DO BEGIN
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

    FOR yy=0, ysize_ave_flats-1 DO BEGIN
        IF (radial_mask(xx, yy) EQ 1) THEN BEGIN
            extract_line = REFORM(flat_scan_ave(xx, yy, line_start:line_end))
	    extract_line = extract_line / prefilter_correction
            extract_line_interp = INTERPOL(extract_line, wavelength_scale_orig(line_start:line_end), wavelengths_new, /SPLINE)
            ;; result = lc_find(intensities, starting pixel num,
            ;;                  ending pixel num, number of fit points)
            ;; result = [core decimal pixel number, core intensity]
            line_center_fit = lc_find(extract_line_interp, $
                                      line_wave_range_numpoint*0.1, $
				      line_wave_range_numpoint*0.9, $
                                      fit_points)
            line_centers(xx, yy, *) = line_center_fit
	    
            counter = counter + 1

        ENDIF
    ENDFOR
ENDFOR

percent_finished   = FLOAT(counter)/mask_points_total*100
seconds_elapsed    = SYSTIME(1,/SECONDS) - start_time
seconds_remaining  = seconds_elapsed/(percent_finished>0.01) * (100 - percent_finished)
time_remaining = STRING(seconds_remaining,FORMAT='(F4.1)') + ' seconds'

ibis_overwrite_line, output_text_length, $
    'Column # ' + STRING(xx, FORMAT='(I4.1)') + $
    '  -- Percent Completed : ' + STRING(percent_finished, FORMAT='(F5.1)') + '%' + $
    '  -- Time Remaining : ' + time_remaining, $
     new_overwrite_length=output_text_length
PRINT

RETURN,line_centers

END
