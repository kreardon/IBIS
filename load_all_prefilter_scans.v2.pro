
; ----------------------------------------------------------------------
;  Script for processing and averaging all the prefilter scans for a 
;  selected prefilter. Includes the separation of a periodic fringe 
;  component from the main prefilter transmission profile.
;  Requires some tuning for different prefilters 
;  (e.g. quality metrics, fringe model).
;
;  requires load_all_prefilter_scans.filter_params.pro
;    to provide filter specific parameters
; ----------------------------------------------------------------------

version = '2.1'

; the uniform wavelength step size onto which to interpolate the prefilter profiles
; -- in Angstroms
wvstep_new                    = 0.02

; the wavelength of the filter to search for
filter      = 6563
filter_name = STRTRIM(filter,2)

; the filter specific paramaters need from here on are now saved in a separate file for clarity.
;
; we have to run the script twice. Once before executing "prepare_prefilter_scans" in order
; define the acceptable ranges for the quality metric parameters. But in this case we need to 
; provide a dummy value for the number of valid prefilter scans we will have. 
prefilter_scan_good_num = 200
do_params               = 1
do_fringes              = 0
@load_all_prefilter_scans.filter_params.pro

; ----------------------------------------------------------------------
;  locate and load all prefilter scans
; ----------------------------------------------------------------------

; Scan IBIS directory for all prefilter scan log files for the requestd filter.
; Assumes log files still have the name 'sintodata.${filter_wavelength}.dat
prefscan_all_files = file_search('~/IBIS','sintodata.' + STRTRIM(filter,2) + '*', COUNT=num_files)
PRINT,'Found ' + STRTRIM(num_files,2) + ' prefilter scan log files for ' + STRTRIM(filter,2) + ' filter.'

; load all prefilter scans for the located prefilter log files, and interpolate them onto a common wavelength grid.
prefscans_even = prepare_prefilter_scans(prefscan_all_files, /Verbose, wvstep_even=wvstep_new, $
                     params_quality_metric=params_quality_metric, quality_metrics=quality_metrics)

; for simplicity, extract the array of measured counts from the returned structure
prefilter_scan_even     = prefscans_even.Counts
wavelength_scale_new    = prefscans_even(0).Wavelengths
prefilter_scan_good_num = N_ELEMENTS(prefscans_even)
num_wv_steps            = N_ELEMENTS(wavelength_scale_new)
; +/- 1 Angstrom
scan_range                    = INDGEN(101) + FIX(num_wv_steps/2.) - 50
; +/- 0.2 Angstrom
scan_range_center             = INDGEN(21) + FIX(num_wv_steps/2.) - 10

; the filter specific paramaters need from here on are now saved in a separate file for clarity
; And now we run the file specific paramater definition file again, this time using the 
; actual value for the number of valid prefilter scans (prefilter_scan_good_num, defined above).
do_params               = 0
do_fringes              = 1
@load_all_prefilter_scans.filter_params.pro

; ----------------------------------------------------------------------
;  normalize and average located prefilter scans
; ----------------------------------------------------------------------

; generate an array with the normalization value for each prefilter scan which is 
; simply the average intensity over then central portion (typically +/- 1 Angstrom) of the profile. 
pref_sum                  = REBIN(REBIN(prefilter_scan_even(scan_range,*),1,prefilter_scan_good_num),num_wv_steps,prefilter_scan_good_num)
; normalize all the scanned profiles to their average central wavelengths
prefilter_scan_even_norm  = prefilter_scan_even/(pref_sum>0.001)
; do the same thing, but only around the center (typically +/- 0.2 Angstrom) of each profile
pref_sum_center           = REBIN(rebin(prefilter_scan_even(scan_range_center,*),1,prefilter_scan_good_num),num_wv_steps,prefilter_scan_good_num)
prefilter_scan_even_norm2 = prefilter_scan_even/(pref_sum_center>0.001)

IF N_ELEMENTS(pref_sum_range_initial) LE 3 THEN pref_sum_range_initial = INDGEN(prefilter_scan_good_num)

; compute the average number of valid points in each wavelength bin
pref_ave_cnts             = REBIN(FLOAT(prefilter_scan_even_norm(*,pref_sum_range_initial) GE 0.001),num_wv_steps,1)
; compute the average number of valid points in each wavelength bin
pref_ave                  = REBIN(prefilter_scan_even_norm(*,pref_sum_range_initial), num_wv_steps, 1)
; computing the average profile over the valid points is easy, since invalid points are set to zero
scan_average              = pref_ave/(pref_ave_cnts>0.001)
; to compute the median we have to first extract the valid points at each wavelength and compute the median over those
scan_median   = FLTARR(num_wv_steps)
FOR nn=0,num_wv_steps-1 DO BEGIN
    pref_ratio            = REFORM( prefilter_scan_even(nn,*) / (pref_sum(nn,*)>0.001))
    scan_median(nn)       = MEDIAN(pref_ratio(WHERE(pref_sum(nn,*) GE 0.01)))
    ; we could also compute the average in this same way, but that is too clunky...
    ; scan_average(nn)      = MEAN(pref_ratio(WHERE(pref_sum(nn,*) GE 0.01)))
ENDFOR

; ----------------------------------------------------------------------
;  isolate, locate, and align small-scale fringes on prefilter profile
; ----------------------------------------------------------------------

fringes = (prefilter_scan_even_norm) / (REBIN(scan_median ,num_wv_steps, prefilter_scan_good_num)>0.001)

fringes_align = fringes
phase_off = FLTARR(prefilter_scan_good_num)

; loop over fringe determination - as we locate and align the fringes, the fringe reference will improve
; which will make the fringe offset determination work better.
FOR fringe_loop = 0,2 DO BEGIN
    fringe_ref         = rebin(fringes_align(fringe_stwv:fringe_ndwv,fringe_stscn:fringe_ndscn),fringe_wvlen,1)
    fringe_ref_linfit  = LINFIT(FINDGEN(fringe_wvlen),fringe_ref,yfit=fringe_ref_fit)
    fringe_ref        /= fringe_ref_fit
    
    freq_wv            = findgen(fringe_wvlen/2.)/fringe_wvlen / wvstep_new
    apodwin            = apod(fringe_wvlen, 1, 0.025, 0.0, 2)
    fringe_ref_fft     = FFT((fringe_ref - 1) * apodwin, -1)
    
    
    FOR scann=0,prefilter_scan_good_num-1 DO BEGIN
        fringe_test_fft  = FFT((fringes(fringe_stwv:fringe_ndwv,scann) - 1) * apodwin, -1)
        cpow             = fringe_ref_fft * CONJ(fringe_test_fft) 
        cpow_re          = SMOOTH(REAL_PART(cpow),5,/EDGE)
        cpow_im          = SMOOTH(IMAGINARY(cpow),5,/EDGE) 
        ; phase_off is converted from a degree to a pixel offset
        phase_off(scann) = MEAN((ATAN(cpow_im,cpow_re))(phase_ref) /!pi/2. / freq_wv(phase_ref)) /  wvstep_new
        fringes_align(*,scann) = SHIFT_VEC_CBC(fringes(*,scann),-phase_off(scann))
    ENDFOR

ENDFOR
; ----------------------------------------------------------------------
;  generate model for small-scale fringes 
; ----------------------------------------------------------------------

; generate a new, slightly finer wavelength scale to generate our sinusoidal model for the fringes
wavelength_scale_new_fine = (FINDGEN(num_wv_steps*2 + 1) - (num_wv_steps-1)) * wvstep_new/2.
; generate a sinusoid to model that matches the observed fringe pattern.
; this probably needs to be done individually for each wavelength.
; paramaters are set above, tuned for each prefilter
model_fringe = SIN((wavelength_scale_new_fine + fringe_off) / fringe_freq * !pi * 2) * fringe_amp + 1

fringes_model_sft     = fringes * 0.0
fringes_model_sft_scl = fringes * 0.0
fringe_linfits        = FLTARR(2,prefilter_scan_good_num)

; shift fringes to wavelength appropriate for each prefilter scan
FOR scann=0,prefilter_scan_good_num-1 DO BEGIN
    ; convert pixel offset to a wavelength shift
    wave_sft                       = phase_off(scann) * wvstep_new
    fringes_model_sft(*,scann)     = INTERPOL(model_fringe,wavelength_scale_new_fine,wavelength_scale_new - wave_sft)
    ; correct any amplitude variations in the fringe pattern from one scan to the next
    ; using a linear fit between the modeled and observed fringes
    valid_fringe_pts               = WHERE(fringes(fringe_stwv:fringe_ndwv,scann) GE 0.95)
    fringes_valid                  = (fringes(fringe_stwv:fringe_ndwv,scann))(valid_fringe_pts)
    fringes_model_valid            = (fringes_model_sft(fringe_stwv:fringe_ndwv,scann))(valid_fringe_pts)
    fringe_linfits(*,scann)        = LINFIT(fringes_model_valid,fringes_valid)
    fringes_model_sft_scl(*,scann) = fringes_model_sft(*,scann) * fringe_linfits(1,scann) + fringe_linfits(0,scann)
ENDFOR

; ----------------------------------------------------------------------
;  remove model fringes from observed profiles
; ----------------------------------------------------------------------

; divide out model fringes
IF use_scaled_fringes THEN BEGIN
    prefilter_scan_even_norm_nofringe  = prefilter_scan_even_norm / fringes_model_sft_scl
ENDIF ELSE BEGIN
    prefilter_scan_even_norm_nofringe  = prefilter_scan_even_norm / fringes_model_sft
ENDELSE

; Use the center-of-gravity method to determine any shifts in the fringe-corrected profiles
prefilter_cog                           = FLTARR(prefilter_scan_good_num)
prefilter_scan_even_norm_nofringe_algn  = prefilter_scan_even_norm_nofringe
FOR scann=0,prefilter_scan_good_num-1 DO BEGIN
    prefilter_cog(scann) = cog1d(wavelength_scale_new(scan_range),prefilter_scan_even_norm_nofringe(scan_range,scann))
ENDFOR
prefilter_cog -= MEAN(prefilter_cog)

; align the profiles based on those shifts
FOR scann=0,prefilter_scan_good_num-1 DO BEGIN
    prefilter_scan_even_norm_nofringe_algn(*,scann) = SHIFT_VEC_CBC(prefilter_scan_even_norm_nofringe(*,scann),prefilter_cog(scann)/0.02)
ENDFOR

; renormalize profiles based on fringe-corrected measurements
prefilter_sum_nofringe             = REBIN(prefilter_scan_even_norm_nofringe(scan_range,*),1,prefilter_scan_good_num)
prefilter_scan_even_norm2_nofringe = prefilter_scan_even_norm_nofringe / REBIN(prefilter_sum_nofringe,num_wv_steps,prefilter_scan_good_num)
prefilter_sum_nofringe2            = REBIN(prefilter_scan_even_norm_nofringe(scan_range_center,*),1,prefilter_scan_good_num)
prefilter_scan_even_norm_cent      = prefilter_scan_even_norm/REBIN(prefilter_sum_nofringe2,num_wv_steps,prefilter_scan_good_num)

IF N_ELEMENTS(pref_sum_range_final) LE 3 THEN pref_sum_range_final = INDGEN(prefilter_scan_good_num)

; sum up the normalized profiles
pref_ave_cnts_nofringe = REBIN(FLOAT(prefilter_scan_even_norm2_nofringe(*,pref_sum_range_final) GE 0.01), num_wv_steps, 1)
pref_ave_nofringe      = REBIN(FLOAT(prefilter_scan_even_norm2_nofringe(*,pref_sum_range_final)),         num_wv_steps, 1)
; and construct an averaged, normalized profile
pref_ave_ref           = pref_ave_nofringe / (pref_ave_cnts_nofringe>0.001)
pref_ave_ref          /= MAX(pref_ave_ref)

;alternatively, construct the median profile, by selecting valid points at each wavelength position.
pref_med_ref = FLTARR(num_wv_steps)
FOR nn=0,num_wv_steps-1 DO BEGIN
    goodp = WHERE(prefilter_scan_even_norm2_nofringe(nn,pref_sum_range_final) GE 0.05,numgood) + min(pref_sum_range_final)
    IF NUMGOOD GE 1 THEN pref_med_ref(nn) = MEDIAN(prefilter_scan_even_norm2_nofringe(nn,goodp),/EVEN)
ENDFOR
pref_med_ref /= MAX(pref_med_ref)

; make copies of the average prefilter profile and the fringe profile with filter-specific variable names
valid_pts_out = WHERE( pref_ave_cnts_nofringe GE 3./prefilter_scan_good_num)
execute_out =  EXECUTE('prefilt' + filter_name + '_ref_main        = pref_ave_ref(valid_pts_out)')
execute_out =  EXECUTE('prefilt' + filter_name + '_ref_wvscl       = wavelength_scale_new(valid_pts_out)')
execute_out =  EXECUTE('prefilt' + filter_name + '_fringe          = INTERPOL(model_fringe, wavelength_scale_new_fine, wavelength_scale_new(valid_pts_out),/SPLINE)')
execute_out =  EXECUTE('prefilt' + filter_name + '_ref_interval    = [prefscans_even[MIN(pref_sum_range_final)].scan_date_text, prefscans_even[MAX(pref_sum_range_final)].scan_date_text]')
CALDAT,SysTime(/Julian,/UTC),mo,dd,yr
provenance  = '''Calculated using load_all_prefilter_scans.pro, v' + version + ', ' + STRING(yr,mo,dd,FORMAT='(%"%4i-%2.2i-%2.2i")') + ', by K. Reardon (kevinpreardon@gmail.com)'''
execute_out =  EXECUTE('prefilt' + filter_name + '_ref_source      = ' + provenance)


END
