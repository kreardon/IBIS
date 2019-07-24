FUNCTION load_calibrate_image, filename, extension, channel, wavelength_nb=wavelength_nb, dark_file=dark_file, gain_file=gain_file, $
                               rotate_array=rotate_array, keep_size=keep_size, calibration_location=calibration_location, target_scale=target_scale

;+
; NAME:
;	  load_calibrate_image
;
; PURPOSE:
;	  
; EXPLANATION:
;	  
; CALLING SEQUENCE:
;	  calibrated_image = load_calibrate_image(filename, extension)
;
; INPUTS:
;	  filename     = 
;	  extension    = 
; 
; OPTIONAL INPUT KEYWORDS:
;     rotate_array = 
;     target_scale =
;
; OPTIONAL OUTPUT KEYWORD:
;	
; RESULTS:
;
; EXAMPLE:
;	  
; COMMON BLOCKS:
;     none	
;
; PROCEDURE:
;
; NOTES:
;
; MODIFICATION HISTORY:
;	Written, Kevin Reardon, National Solar Observatory, 2018-2019
;-

; set up desired orientation and scale 
IF NOT KEYWORD_SET(rotate_array) THEN rotate_array='solar_north'
IF NOT KEYWORD_SET(target_scale) THEN target_scale=0.096

; determine whether input image is narrowband or whitelight
IF StrMatch(channel,'*nb*') THEN channel_id = 'nb' ELSE channel_id = 'wl'

; read in image to be loaded
image_array     = readfits(filename, exten=extension, image_hdr, /Silent)
image_array_sz  = SIZE(image_array,/str)
image_array_szx = image_array_sz.dimensions[0]
image_array_szy = image_array_sz.dimensions[1]

image_dateobs    = sxpar(image_hdr, 'DATE-OBS')
image_dateobs_jd = fits_date_convert(image_dateobs)
date_str         = STRMID(image_dateobs,0,10)

IF NOT keyword_set(wavelength_nb) THEN wavelength_nb = ROUND(sxpar(image_hdr, 'WAVELNTH'))

cal_params = load_calibration_info(date_str, 'ibis_' + channel_id)

IF NOT Keyword_Set(calibration_location) THEN BEGIN
    repository_location = File_Dirname(Routine_Filepath(/Either),/Mark)
    calibration_location = repository_location + 'calibration_files/'
ENDIF

; determine if calibration files are valid and accessible
IF N_ELEMENTS(dark_file) EQ 1 THEN dark_file_use = dark_file ELSE dark_file_use = cal_params.dark_file
dark_file_valid = FILE_TEST(dark_file_use)
IF FILE_TEST(calibration_location) THEN dark_file_use = calibration_location + '/' + dark_file_use
IF N_ELEMENTS(gain_file) EQ 1 THEN gain_file_use = gain_file ELSE gain_file_use = cal_params.gain_file
gain_file_valid = FILE_TEST(gain_file_use)
IF FILE_TEST(calibration_location) THEN gain_file_use = calibration_location + '/' + gain_file_use

IF dark_file_valid NE 1 THEN BEGIN
    dark_files = File_Search(calibration_location,'*dark*' + channel_id + '*', count=dark_count)
    dark_file_use  = dark_files[0]
ENDIF

IF gain_file_valid NE 1 THEN BEGIN
    gain_files = File_Search(calibration_location,'*gain*' + channel_id + '*', count=dark_count)
    gain_file_use  = gain_files[0]
ENDIF

; restore identified calibration files
RESTORE,/VE,dark_file_use
RESTORE,/VE,gain_file_use

; set defined calibration array to a common name
res = EXECUTE('dark_cal = ' + cal_params.dark_name)
res = EXECUTE('gain_cal = ' + cal_params.gain_name)

; apply dark and flat correction to data
image_array = (image_array - dark_cal) / wl_gain

; rotate image to roughly match solar Cartesian coordinates (North up, East left)
image_array = ROTATE(image_array,cal_params.transpose)

; Determine rotation angle to be applied to image

rot_fit_params     = cal_params.rot_to_solnorth
rot_fit_params_num = N_ELEMENTS(rot_fit_params_num)
image_time_scale   = image_dateobs_jd - cal_params.rot_to_sol_reftime
image_rot_angle    = 0.0

FOR nn=0,rot_fit_params_num-1 DO BEGIN
    image_rot_angle += rot_fit_params[nn] * (image_time_scale)^nn
ENDFOR

IF STRMATCH(rotate_array,'*solar*') THEN BEGIN
    image_array = ROT(image_array,image_rot_angle,CUBIC=-0.5)
    PRINT,'Rotating: ', cal_params.rot_to_solnorth[1]
ENDIF ELSE IF STRMATCH(rotate_array,'*grid*') THEN BEGIN
    image_array = ROT(image_array,cal_params.rot_to_grid[0],CUBIC=-0.5)
    PRINT,'Rotating: ', cal_params.rot_to_grid[0]
ENDIF ELSE IF (size(rotate_array,/str)).TYPE_NAME NE 'STRING' THEN BEGIN
    image_array = ROT(image_array,rotate_array,CUBIC=-0.5)
    PRINT,'Rotating: ', rotate_array
ENDIF

IF channel_id EQ 'wl' THEN BEGIN
    plate_scale_im = cal_params.plate_scale[0:1]
ENDIF ELSE BEGIN
    ; filter_idx_select = where(fix(cal_params.plate_scale[2,*])  eq wavelength)
    filter_idx_select = get_closest(cal_params.plate_scale[2,*],wavelength)
    plate_scale_im = cal_params.plate_scale[0:1, filter_idx_select]
ENDELSE

target_scale_ratio = plate_scale_im / target_scale
target_scale_pixel = ROUND([image_array_szx,image_array_szy] * target_scale_ratio)
target_scale_diff  = target_scale_pixel - [image_array_szx,image_array_szy]
image_array        = CONGRID(image_array, target_scale_pixel[0], target_scale_pixel[1], CUBIC=-0.5)

IF KEYWORD_SET(keep_size) THEN BEGIN

    pix_rangex_arr = [0,target_scale_pixel[0]]
    pix_rangex_out = [0,image_array_szx-1]
    IF target_scale_diff[0] GT 0 THEN pix_rangex_arr = [ROUND((target_scale_diff[0])/2.),(target_scale_diff[0])/2. + image_array_szx-1]
    IF target_scale_diff[0] LT 0 THEN pix_rangex_out = [ROUND(target_scale_diff[0]/2.),       ROUND(target_scale_diff[0]/2.) + target_scale_pixel[0]]

    pix_rangey_arr = [0,target_scale_pixel[1]]
    pix_rangey_out = [0,image_array_szy-1]
    IF target_scale_diff[1] GT 0 THEN pix_rangey_arr = [ROUND((target_scale_diff[1])/2.),(target_scale_diff[1])/2. + image_array_szy-1]
    IF target_scale_diff[0] LT 0 THEN pix_rangey_out = [ROUND(target_scale_diff[1]/2.),       ROUND(target_scale_diff[1]/2.) + target_scale_pixel[1]]

    image_array_keep = FLTARR(image_array_szx, image_array_szy) + MEDIAN(image_array)
    image_array_keep[pix_rangex_out[0]:pix_rangex_out[1],pix_rangey_out[0]:pix_rangey_out[1]] = $
         image_array[pix_rangex_arr[0]:pix_rangex_arr[1],pix_rangey_arr[0]:pix_rangey_arr[1]]
    image_array = image_array_keep
ENDIF

RETURN,image_array

END
