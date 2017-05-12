PRO average_series_images, datapath, datatype, detector_name, skip_series=skip_series

;datapath  = T110820Apr2017

;DataTypes = 'DarkCalibration'
;DataTypes = 'FlatFieldCalibration'
;DataTypes = 'GridImages'
;DataTypes = 'TargetImages'
;DataTypes = 'OtherCalibration'
IF NOT KEYWORD_SET(datatype) THEN DataTypes = 'DarkCalibration' ELSE DataTypes = datatype

;Detector  = 'andor1'
IF NOT KEYWORD_SET(detector_name) THEN Detector = 'andor1' ELSE Detector = detector_name

IF Detector EQ 'andor1' THEN Channel = 'spectral'
IF Detector EQ 'andor2' THEN Channel = 'whitelight'

data_mask = BYTARR(1000,1000)
data_mask(50:950,50:950) = 1

CASE DataTypes OF
    'DarkCalibration'      : selection_mode        = 'combineall'
    'FlatFieldCalibration' : IF Detector EQ 'andor1' THEN selection_mode = 'byfilter+wave' ELSE selection_mode = 'combineall'
    'GridImages'           : IF Detector EQ 'andor1' THEN selection_mode = 'byfilter' ELSE selection_mode = 'combineall'
    'TargetImages'         : IF Detector EQ 'andor1' THEN selection_mode = 'byfilter' ELSE selection_mode = 'combineall'
    'OtherCalibration'     : IF Detector EQ 'andor1' THEN selection_mode = 'byfilter' ELSE selection_mode = 'combineall'
    'ScienceObservation'   : IF Detector EQ 'andor1' THEN selection_mode = 'byfilter+wave' ELSE selection_mode = 'combineall'
    ELSE                   : IF Detector EQ 'andor1' THEN selection_mode = 'byfilter' ELSE selection_mode = 'combineall'
ENDCASE
;selection_mode        = 'byfilter+wave'
;selection_mode        = 'byfilter'
;selection_mode        = 'combineall'

basedir = '/net/bonneville/export/dstdata/' + Detector + '/' + datapath + '/'
data_series_all = FILE_SEARCH(basedir + DataTypes + '/201*_*',count=num_series)

data_series_cor = data_series_all
FOR mm=0,N_ELEMENTS(skip_series)-1 DO BEGIN
    bad_series = WHERE(STRMATCH(data_series_cor, '*' + skip_series[mm] + '*'))
    IF MIN(bad_series) GE 0 THEN BEGIN
        good_index = good_frame_index(N_ELEMENTS(data_series_cor), bad_series)
        data_series_cor = data_series_cor[good_index]
    ENDIF
ENDFOR
num_series = N_ELEMENTS(data_series_cor)
data_series = data_series_cor
print,'Found the following data series -'
print,FORMAT='("    ",A0)',data_series

;num_series = 1
do_write = 1

modulation_options = ['I', 'I+Q', 'I+U', 'I+V', 'I-Q', 'I-U', 'I-V', '4S1', '4S2', '4S3', '4S4']

FOR ds=0,num_series-1 DO BEGIN
    vee_log        = read_vee_log( FILE_SEARCH(data_series[ds],'log*txt') )
    wavelen        = vee_log.wavelength
    filter_pos     = vee_log.filter_wheelpos
    modulation     = vee_log.modulation
    modulation_num = INTARR(N_ELEMENTS(modulation)) - 1
    for mo=0,N_ELEMENTS(modulation_options) - 1 DO BEGIN
	modulation_num(WHERE(modulation EQ  modulation_options(mo),/NULL)) = mo
    ENDFOR

    series_combo = LONARR(N_ELEMENTS(vee_log.wavelength))
    IF STRMATCH(selection_mode,'*wav*')  THEN series_combo  += wavelen * 1
    IF STRMATCH(selection_mode,'*filt*') THEN series_combo  += filter_pos * 10
    IF STRMATCH(selection_mode,'*mod*')  THEN series_combo  += modulation_num * 100
    filter_wave_mod_combo = modulation_num * 100 + filter_pos * 10 + wavelen
    filter_wave_combo	  =			   filter_pos * 10 + wavelen
    filter_combo	  =			   filter_pos * 10

    series_uniq          = series_combo(UNIQ(series_combo,SORT(series_combo)))
    filter_wave_mod_uniq = filter_wave_mod_combo(UNIQ(filter_wave_mod_combo,SORT(filter_wave_mod_combo)))
    filter_wave_uniq	 = filter_wave_combo(UNIQ(filter_wave_combo,SORT(filter_wave_combo)))
    filter_uniq	         = filter_combo(UNIQ(filter_combo,SORT(filter_combo)))
    num_uniq    = N_ELEMENTS(series_uniq)
    filter_list = vee_log.filter(UNIQ(filter_combo,SORT(filter_combo)))
    PRINT,'Found ' + STRTRIM(N_ELEMENTS(filter_uniq), 2) + ' unique Filter(s).'
    PRINT,REFORM(filter_list,1,N_ELEMENTS(filter_list)),FORMAT='(I8.4)'
    PRINT,'Found ' + STRTRIM(N_ELEMENTS(filter_wave_uniq), 2) + ' unique Filter/Wavelength combinations.'
    PRINT,'Found ' + STRTRIM(N_ELEMENTS(filter_wave_mod_uniq), 2) + ' unique Filter/Wavelength/Modulation combinations.'
    PRINT,'==>  Will sum ' + STRTRIM(N_ELEMENTS(series_uniq), 2) + ' unique series combination(s).'

    series_cnt = INTARR(num_uniq)
    series_cnt_log = INTARR(num_uniq)
    FOR nn=0,num_uniq - 1 DO series_cnt_log(nn) = TOTAL((series_combo EQ series_uniq(nn)))

    data_files = FILE_SEARCH(data_series[ds],'s*fits*',COUNT=num_data_files)
    ;sometimes it's necessary do skip a data file on an adhoc basis
    ;data_files = data_files[0:8]
    num_data_files = N_ELEMENTS(data_files)

    fits_open,data_files(0),fcb
    fits_read,fcb,main_data,main_header,exten=0
    fits_read,fcb,data,header,exten=1
    FITS_CLOSE,fcb
    imsiz = SIZE(data)
    hdrsiz = SIZE([main_header, header])

    series_ave = FLTARR(imsiz(1),imsiz(2),num_uniq)
    series_ave_input = STRARR(8, num_uniq, MAX(series_cnt_log))
    series_ave_stats = STRARR(5, num_uniq, MAX(series_cnt_log))
    images_info = STRARR(2,fcb.NEXTEND, num_data_files)
    headers_all = STRARR(hdrsiz(1),fcb.NEXTEND,num_data_files)
    series_ave_wv = FLTARR(num_uniq)
    series_ave_mod = STRARR(num_uniq)
    series_seq	   = FLTARR(1000,1000,num_data_files)
    incnt = 0

    data_includes_bias = 1
    data_bias	       = 307


    FOR fil=0,num_data_files-1 DO BEGIN
        data_file_info = file_info(data_files(fil))
        IF data_file_info.size GE 2880 THEN BEGIN
	fits_open,data_files(fil),fcb
        IF data_file_info.size GE fcb.NEXTEND * 2e6 THEN BEGIN
	fits_read,fcb,main_data,main_header,exten_no=0
	;lun = fxposit(data_files[fil],0,/ReadOnly)
	;nodata = mrdfits(lun, 0, main_header, status=status,/SILENT)
	camera_id	   = SXPAR(main_header,'SER_NUM', /SILENT)
	print,'    ' + data_files(fil)
	

	FOR ext=1,fcb.NEXTEND DO BEGIN
            fits_read,fcb,data,header,exten=ext,/NOSCALE
            IF N_ELEMENTS(data) GE 1000 THEN status=1
            ;data = mrdfits(lun, 0, header, status=status,/SILENT)
            ; if status from mrdfits is less than zero, no data was read, so skip this extension
            IF status GE 0 THEN BEGIN
		data_date	   = SXPAR(header, 'DATE-OBS', /SILENT)
		linearity_corrected = 0
		IF (camera_id EQ 2748) THEN BEGIN
		    data = ibis_linearity_correction_andor(data, camera_id=camera_id, $
                                                     data_date=data_date, verbose=verbose, $
                                                     data_includes_bias=data_includes_bias,data_bias=data_bias, $
                                                     lin_ptr=lin_ptr)
                    linearity_corrected = 1
                ENDIF

		rel_wave_in   = sxpar(header,'REL_WAVE')
		filter_pos_in = sxpar(header,'FILTER')
		stokes_in     = sxpar(header,'STOKES')
		image_time    = sxpar(header,'DATE-OBS')
		image_exptime = sxpar(header,'EXPTIME')
		stokes_num = WHERE(STRMATCH(modulation_options,stokes_in))

		;image_filt_wave_combo = FLOAT(stokes_num)*100 + FLOAT(filter_pos_in)*10. + FLOAT(rel_wave_in)
		image_combo = 0L
		IF STRMATCH(selection_mode,'*wav*')  THEN image_combo  += FLOAT(rel_wave_in) * 1
		IF STRMATCH(selection_mode,'*filt*') THEN image_combo  += FLOAT(filter_pos_in) * 10
		IF STRMATCH(selection_mode,'*mod*')  THEN image_combo  += FLOAT(stokes_num) * 100
		im_match = get_closest(series_uniq, image_combo)

		series_ave(*,*,im_match) += data
		series_ave_input(*,im_match, series_cnt(im_match)) =  [data_files(fil), STRTRIM(ext,2), $
                                                        STRTRIM(filter_pos_in,2), $
                                                        STRTRIM(rel_wave_in,2), $
                                                        STRTRIM(stokes_in,2), $
                                                        STRTRIM(image_time, 2), STRTRIM(image_exptime,2), STRTRIM(linearity_corrected,2)]
	        
	        ibis_area_statistics,data,data_mask,all_stats=image_stats
		series_ave_stats(*,im_match, series_cnt(im_match)) =  image_stats[0:4]

		;IF im_match EQ 12 THEN series_seq(*,*,series_cnt(im_match)) = data
		series_cnt(im_match)	 += 1
		images_info(*,ext-1,fil) = [data_files(fil), STRTRIM(ext,2)]
		headers_all(*,ext-1,fil) = [main_header, header]
            ENDIF
	ENDFOR
	FITS_CLOSE,fcb
	;Free_LUN,lun
	TV,GmaScl(series_ave(*,*,0),gamma=4)

        ENDIF
        ENDIF

    ENDFOR

    FOR unq=0,num_uniq-1 DO BEGIN
	series_ave(*,*,unq) /= series_cnt(unq)
    ENDFOR

    IF do_write THEN BEGIN
	series_split = (STRSPLIT(data_series[ds],'/',/EXTRACT))
	series_name = series_split(WHERE(STRMATCH(series_split,'*20[0-9]*[0-9]_[0-9]*')))
	series_type = series_split(WHERE(STRMATCH(series_split,'[FDTGSIO]*[nsg]')))
	output_file = series_type + '.' + Channel + '.' + selection_mode + '.' + series_name + '.series.ave.sav'
	PRINT,'Writing output file: ' + output_file
	SAVE,FILENAME=output_file,$
             series_ave, series_cnt, series_ave_input, series_ave_stats, images_info, headers_all,vee_log, /COMP
   ENDIF

ENDFOR

END

