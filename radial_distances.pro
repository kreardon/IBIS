;-------------------------------------------------------------
;+
; NAME:
;       radial_distances
; PURPOSE:
;       To compute a 2-D map of distances from a given central point
; CATEGORY:
; CALLING SEQUENCE:
;       distance_map = radial_distances(array_size, central_point)
; INPUTS:
;       array_size    = a two-element array defining the X and Y dimensions
;                       of the desired output map
;       central_point = a two-element array defining the coordinated of the
;                       central point from which to measure distances. If the 
;                       central_point is not given, the center of the output
;                       map is taken as the central point.
;                       Can be fractional or even outside of the output map.
;                       Distances are measured from the lower left corner of
;                       the first pixel (0,0). 
; KEYWORD PARAMETERS: 
;       fft_shift     = If this keyword is set (as a scalar), the output map 
;                       is shifted to place the central point for the output
;                       distances at the four corners of the array, as is 
;                       appropriate when looking the the calculated FFT power
;                       for a 2-D array.                     
; OUTPUTS:
; COMMON BLOCKS:
;       
; NOTES:
;       If 'array_size' has only one element, a square array, with both axes
;       equal to the array_size parameter, is generated.
;
;       If 'array_size' has more than two elements, the program will use the 
;       second and third elements for calculating the array size. This facilitates
;       simply passing the output from SIZE() to this program.
;
;       The shift applied by fft_shift is only an integer shift, not fractional.
;
; MODIFICATION HISTORY:
;       K. Reardon, 2002, Initial Implementation
;-
;-------------------------------------------------------------
FUNCTION radial_distances, array_size, central_point, fft_shift=fft_shift

IF N_ELEMENTS(array_size) GT 2 THEN BEGIN
    sx = array_size(1)
    sy = array_size(2)
ENDIF ELSE IF N_ELEMENTS(array_size) EQ 2 THEN BEGIN
    sx = array_size(0)
    sy = array_size(1)
ENDIF ELSE IF N_ELEMENTS(array_size) EQ 1 THEN BEGIN
    sx = array_size(0)
    sy = array_size(0)
ENDIF ELSE BEGIN
    PRINT, 'ERROR - array size must be defined!'
    RETURN,-1
ENDELSE

IF (N_ELEMENTS(central_point) EQ 0) THEN central_point = [(sx/2.), (sy/2.)]

; Calculate the X and Y distances from the central point
xdist = REBIN(FINDGEN(sx)-central_point(0), sx, sy)
ydist = REBIN(REFORM(FINDGEN(sy)-central_point(1), 1, sy), sx, sy)

radist = SQRT(xdist^2 + ydist^2)

; Apply the shift into FFT coordinates if desired
IF N_ELEMENTS(fft_shift) EQ 1 THEN BEGIN
    IF KEYWORD_SET(fft_shift) THEN BEGIN
        radist = SHIFT(radist, ROUND(-central_point(0)), ROUND(-central_point(1)))
    ENDIF
ENDIF ELSE IF N_ELEMENTS(fft_shift) EQ 2 THEN BEGIN
        radist = SHIFT(radist,ROUND(fft_shift(0)), ROUND(fft_shift(1)))
ENDIF

RETURN,radist

END  ; radial_distances
