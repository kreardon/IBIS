FUNCTION frame_transfer_cor, imagein, exptime, vert_sft_time, $
                             readout_transfer_only=readout_transfer_only, $
                             binning=binning, ratio_correction=ratio_correction
;+
; NAME:
;	FRAME_TRANSFER_COR
; PURPOSE:
;	Return the current operating system as in !VERSION.OS_FAMILY 
;
; CALLING SEQUENCE
;	corrected_image = FRAME_TRANSFER_COR(input_image,exptime,vert_sft_time)
; INPUTS: 
;	input_image = raw image from camera, typically with dark/bias already 
;                     subtracted off
;   exptime     = nominal exposure time for the image
;   vert_sft_time = transfer time for a single row of the image array
; KEYWORD PARAMETERS:
;   readout_transfer_only = set this keyword to calculate smearing only
;                           due to illumination after exposure is finished.
;                           Default is to correct smearing across the full chip
;   binning = set the binning factor (in the readout direction) used 
;             in acquiring the image. Should be an integer value (e.g. 1,2,4,...)
;   ratio_correction = an ad hoc correction factor by which to adjust the 
;                      exposure-time-to-transfer-time ratio to achieve the best correction. 
; OUTPUTS:
;	result - image with smearing due to illumination during shift
;            of frame into masked area removed
; NOTES:
;    The units used for the exposure time and vertical shift time are 
;    not critical (seconds, milliseconds, etc.), but they must be 
;    in the _same_ units for the exposure-time-to-transfer-time ratio 
;    to be determined
;
;    The smearing due to frame transfer may occur only after the exposure,
;    which produces smears on just one "side" of a bright object, or 
;    can occur both prior to and after exposure, which produces smears
;    on both sides of a bright object. Set readout_transfer_only=1 for the
;    former case.
;	
; REVISION HISTORY:
;	Written,  K. Reardon    
;-

IF (N_ELEMENTS(binning) EQ 0)             THEN binning=1
IF (binning LT 1)                         THEN binning=1
IF (N_ELEMENTS(ratio_correction) EQ 0)    THEN ratio_correction =1
IF N_ELEMENTS(readout_transfer_only) EQ 0 THEN readout_transfer_only=0

sft_time_ratio  = vert_sft_time/FLOAT(exptime)
sft_time_ratio *= ratio_correction
exp_time_ratio  = 1./sft_time_ratio 
IF exp_time_ratio LE 200 THEN $
    PRINT,'Warning! Image exposure-to-transfer-time ratio is strangely small - ' + STRING(exp_time_ratio, FORMAT='(F7.2)')

num_row       = N_ELEMENTS(imagein(0,*))
num_col       = N_ELEMENTS(imagein(*,0))
num_row_array = num_row * binning

imagein_cor = imagein
row_cor_sum = FLTARR(num_col)

IF KEYWORD_SET(readout_transfer_only) THEN BEGIN

    FOR row = 0,num_row-1 DO BEGIN
        row_cor            = imagein(*,row)
        row_cor            = row_cor - (row_cor_sum * sft_time_ratio)
        row_cor_sum       += row_cor
        imagein_cor(*,row) = row_cor
    ENDFOR

ENDIF ELSE BEGIN
    row_ave          = REBIN(imagein, num_col, 1)
    imagein_row_ave  = REBIN(row_ave, num_col, num_row)
    imagein_row_ave *= sft_time_ratio * num_row_array
    imagein_cor      = imagein - imagein_row_ave
ENDELSE

RETURN,imagein_cor

END
