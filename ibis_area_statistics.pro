;-------------------------------------------------------------
;+
; NAME:
;       ibis_area_statistics
; PURPOSE:
;       Calculate statistical properties for some given subset  
;       of points in the input array
; CATEGORY:
; CALLING SEQUENCE:
;       ibis_area_statistics, input_array, input_mask
; INPUTS:
;       input_array = array of input data values
;       input_mask  = definition of points over which to calculate
;                     the statistical properties. Can either be given
;                     as:
;                     a) a mask in an array with same dimensions as 
;                        input_array. Any points with non-zero values
;                        in the mask will be used in calculating the
;                        statisical values.
;                     b) a list of indices (generally in a 1-D array)
;                        of the points to be included in calculating
;                        the statistical values. In this case, the
;                        MASK_POINTS_INDICES keyword must be set.
; KEYWORD PARAMETERS:
;       mask_points_indices  = indicates that input_mask is a list of indices
;                              desccribing the points over which to calculate
;                              the statiscal values
;       print_stats          = print out the calculated results either in a short
;                              form (print_stats = 1) or a longer form (print_stats=2).
; OUTPUT                       
;       mean_calculated      = the calculated mean of the selected points - see MOMENT man page
;       median_calculated    = the calculated median of the selected points - see MEDIAN man page
;       min_calculated       = the minimum value in the set of selected points
;       max_calculated       = the maximum value in the set of selected points
;       stdev_calculated     = the calculated standard deviation of the selected points - see MOMENT man page
;       num_selected_points  = the number of selected points over which the statistics are calculated
;       all_stats_calculated = a single 6-element array holding all the above values
; COMMON BLOCKS:
;       
; NOTES:
; MODIFICATION HISTORY:
;       K. Reardon, May 2004, Initial Implementation
;       K. Reardon, March 2005, Changed output order, add structure output
;-
;-------------------------------------------------------------

PRO ibis_area_statistics, input_array, input_mask, $
                          mask_points_indices  = mask_points_indices,  $
                          mean_calculated      = mean_calculated,   $
                          median_calculated    = median_calculated, $
                          min_calculated       = min_calculated,    $
                          max_calculated       = max_calculated,    $
                          stdev_calculated     = stdev_calculated,  $
                          num_selected_points  = num_selected_points, $
                          all_stats_calculated = all_stats_calculated, $
                          print_stats          = print_stats, $
                          structure_output     = structure_output

IF NOT KEYWORD_SET(print_stats) THEN print_stats = 0

num_points_mask    = 0
num_elements_array = N_ELEMENTS(input_array)
selected_points    = 0
CASE SIZE(input_array,/TNAME) OF
    'DOUBLE': use_double=1 
    ELSE:     use_double=0
ENDCASE

; if there is any error in the input mask or array we will return all zeroes
mean_calculated      = 0
median_calculated    = 0
min_calculated       = 0
max_calculated       = 0
stdev_calculated     = 0
num_selected_points  = 0

IF KEYWORD_SET(mask_points_indices) THEN BEGIN
    mask_points      = input_mask
    num_points_mask  = N_ELEMENTS(mask_points)
ENDIF ELSE BEGIN
    mask_points      = WHERE(input_mask, num_points_mask)
ENDELSE
mask_points_valid = WHERE((mask_points GE 0) AND (mask_points LT num_elements_array), valid_points)

IF (valid_points GE 1) THEN BEGIN
    mask_points = mask_points(mask_points_valid)
    selected_points  = input_array(mask_points)
    finite_points    = WHERE(FINITE(selected_points),finite_points_count)
    IF finite_points_count GE 1 THEN BEGIN
        selected_points  = selected_points(finite_points)

       ; need to use a call to MOMENT to calculate mean in order to avoid conflicts
       ; with JHU/APL IDL library's 'mean.pro', which doesn't support DOUBLE or NAN
       moment_calculated    = MOMENT(selected_points, SDEV=stdev_calculated,$
           DOUBLE=use_double, MAXMOMENT = 2)
       mean_calculated      = moment_calculated(0)
       median_calculated    = MEDIAN(selected_points, /EVEN)
       min_calculated       = MIN(selected_points,     MAX=max_calculated)
       num_selected_points  = finite_points_count
    ENDIF
ENDIF

IF KEYWORD_SET(structure_output) THEN BEGIN
    all_stats_calculated = CREATE_STRUCT (            $
                     'mean',        mean_calculated,    $
                     'median',      median_calculated,  $
                     'min',         min_calculated,     $
                     'max',         max_calculated,     $
                     'stdev',       stdev_calculated,   $
                     'numpoints',   num_selected_points $
                     )
ENDIF ELSE BEGIN
    all_stats_calculated = [mean_calculated, median_calculated, $
                            min_calculated, max_calculated, $
                            stdev_calculated, num_selected_points]
ENDELSE

IF ((print_stats GT 0) AND (print_stats LT 2)) THEN BEGIN
    
    PRINT, FORMAT='(%"Mean: %f, Median: %f, Min: %f, Max: %f, Std Dev: %f")', $
           mean_calculated, median_calculated, min_calculated, max_calculated, $
           stdev_calculated

ENDIF ELSE IF (print_stats GE 2) THEN BEGIN

    PRINT, 'Number Mask Points : ', STRING(N_ELEMENTS(selected_points))
    PRINT, 'Mean, Median       : ', mean_calculated, median_calculated
    PRINT, 'Min, Max           : ', min_calculated, max_calculated
    PRINT, 'Standard Deviation : ', stdev_calculated

ENDIF


END
