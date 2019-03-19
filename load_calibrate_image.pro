FUNCTION load_calibrate_image, filename, extension, channel, wavelength_nb=wavelength_nb, dark_file=dark_file, gain_file=gain_file, $
                               rotate_array=rotate_array, keep_size=keep_size, calibration_location=calibration_location, target_scale=target_scale

IF NOT KEYWORD_SET(rotate_array) THEN rotate_array='solar_north'
IF NOT KEYWORD_SET(target_scale) THEN target_scale=0.096

IF StrMatch(channel,'*nb*') THEN channel_id = 'nb' ELSE channel_id = 'wl'

cal_params = load_calibration_info(date_str, 'ibis_' + channel_id)

IF NOT Keyword_Set(calibration_location) THEN BEGIN
    repository_location = File_Dirname(Routine_Filepath(/Either),/Mark)
    calibration_location = repository_location + 'calibration_files/'
ENDIF

IF N_ELEMENTS(dark_file) EQ 1 THEN dark_file_valid = FILE_TEST(dark_file) ELSE dark_file_valid=0
IF N_ELEMENTS(gain_file) EQ 1 THEN gain_file_valid = FILE_TEST(gain_file) ELSE gain_file_valid=0

IF dark_file_valid NE 1 THEN BEGIN
    dark_files = File_Search(calibration_location,'*dark*' + channel_id + '*', count=dark_count)
    dark_file  = dark_files[0]
ENDIF

IF gain_file_valid NE 1 THEN BEGIN
    gain_files = File_Search(calibration_location,'*gain*' + channel_id + '*', count=dark_count)
    gain_file  = gain_files[0]
ENDIF

RESTORE,/VE,dark_file
RESTORE,/VE,gain_file

image_array = readfits(filename, exten=extension, image_hdr, /Silent)
image_array = (image_array - wl_dark) / wl_gain
image_array = ROTATE(image_array,cal_params.transpose)
IF STRMATCH(rotate_array,'*solar*') THEN BEGIN
    image_array = ROT(image_array,cal_params.rot_to_solnorth,CUBIC=-0.5)
ENDIF ELSE IF STRMATCH(rotate_array,'*grid*') THEN BEGIN
    image_array = ROT(image_array,cal_params.rot_to_grid,CUBIC=-0.5)
ENDIF ELSE IF (size(rotate_array,/str)).TYPE_NAME NE 'STRING' THEN BEGIN
    image_array = ROT(image_array,rotate_array,CUBIC=-0.5)
ENDIF

IF channel_id EQ 'wl' THEN BEGIN
    plate_scale_im = cal_params.plate_scale[0:1]
ENDIF ELSE BEGIN
    ;filter_idx_select = where(cal_params.filter_ids eq wavelength)
    filter_idx_select = where(fix(cal_params.plate_scale[2,*])  eq wavelength)
    plate_scale_im = cal_params.plate_scale[0:1, filter_idx_select]
ENDELSE

target_scale_ratio = target_scale / plate_scale_im
target_scale_pixel = ROUND([1000,1000] * target_scale_ratio)
target_scale_diff  = target_scale_pixel - [1000,1000]
image_array        = CONGRID(image_array, target_scale_pixel[0], target_scale_pixel[1], CUBIC=-0.5)

IF KEYWORD_SET(keep_size) THEN BEGIN

    pix_rangex_arr = [0,target_scale_pixel[0]]
    pix_rangex_out = [0,999]
    IF target_scale_diff[0] GT 0 THEN pix_rangex_arr = [ROUND((target_scale_diff[0])/2.),(target_scale_diff[0])/2. + 999]
    IF target_scale_diff[0] LT 0 THEN pix_rangex_out = [ROUND(target_scale_diff[0]/2.),       ROUND(target_scale_diff[0]/2.) + target_scale_pixel[0]]

    pix_rangey_arr = [0,target_scale_pixel[1]]
    pix_rangey_out = [0,999]
    IF target_scale_diff[1] GT 0 THEN pix_rangey_arr = [ROUND((target_scale_diff[1])/2.),(target_scale_diff[1])/2. + 999]
    IF target_scale_diff[0] LT 0 THEN pix_rangey_out = [ROUND(target_scale_diff[1]/2.),       ROUND(target_scale_diff[1]/2.) + target_scale_pixel[1]]

    image_array_keep = FLTARR(1000,1000) + MEDIAN(image_array)
    image_array_keep[pix_rangex_out[0]:pix_rangex_out[1],pix_rangey_out[0]:pix_rangey_out[1]] = $
         image_array[pix_rangex_arr[0]:pix_rangex_arr[1],pix_rangey_arr[0]:pix_rangey_arr[1]]
    image_array = image_array_keep
ENDIF

RETURN,image_array

END
