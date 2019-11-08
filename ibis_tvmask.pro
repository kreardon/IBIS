;-------------------------------------------------------------
;+
; NAME:
;       ibis_tvmask
; PURPOSE:
;       To display an image using a scaling determined only from
;       points located within the supplied mask
; CATEGORY:
; CALLING SEQUENCE:
;       ibis_tvmask, input_array, input_mask
; INPUTS:
;       input_array = image to display
;       input_mask  = mask defining points (NE 0) which should be used in 
;                     calculating minimum and maximum values for display
;       X           = position in x-axis for image
;       Y           = position in y-axis for image
; KEYWORD PARAMETERS:
;       scale_range = the nmaximum and minimum values determined for the display
; OUTPUTS:
; COMMON BLOCKS:
;       
; NOTES:
; MODIFICATION HISTORY:
;       K. Reardon, Feb 2004, Initial Implementation
;-
;-------------------------------------------------------------
PRO ibis_tvmask, input_array, input_mask, X, Y, $
                 scale_range=scale_range, _EXTRA=_EXTRA

IF NOT KEYWORD_SET(scale_range) THEN scale_range=1.
IF N_ELEMENTS(input_mask) LT N_ELEMENTS(input_array) THEN input_mask = input_array GT 0

mask_points = WHERE(input_mask, mask_point_count)

IF mask_point_count GE 1 THEN BEGIN
    minval = MIN(input_array(WHERE(input_mask)), max=maxval)
ENDIF ELSE BEGIN
    minval = MIN(input_array, max=maxval)
ENDELSE

full_range = maxval-minval

IF (N_ELEMENTS(scale_range) EQ 2) THEN BEGIN
    min_scale = minval + full_range * scale_range(0)
    max_scale = maxval - full_range * scale_range(1)
ENDIF ELSE BEGIN
    use_range = full_range * scale_range
    min_scale = minval + full_range * ((1-scale_range)/2.)
    max_scale = maxval - full_range * ((1-scale_range)/2.)
ENDELSE

IF KEYWORD_SET(X) AND KEYWORD_SET(Y) THEN $
    TV,BYTSCL(input_array, MIN=min_scale, max=max_scale), X, Y, _EXTRA=_EXTRA $
ELSE IF KEYWORD_SET(X) AND (NOT KEYWORD_SET(Y)) THEN $
    TV,BYTSCL(input_array, MIN=min_scale, max=max_scale), X, _EXTRA=_EXTRA $
ELSE $
    TV,BYTSCL(input_array, MIN=min_scale, max=max_scale), _EXTRA=_EXTRA 

END
