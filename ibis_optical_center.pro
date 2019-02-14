FUNCTION ibis_optical_center, blueshift, mask=mask, coeff=coeff, $
                              blueshift_fit_matrix=blueshift_fit_matrix, $
                              use_squares=use_squares, max_degree=max_degree
;+
; Name 
;         IBIS_OPTICAL_CENTER
; Purpose 
;         Find optical center, e.g. the minimum of the smoothed blueshift array
; Calling Sequence 
;         vector = ibis_optical_center(blueshiftarray)
; Inputs
;         blueshiftarray: (float) 2d-array of blueshifts
; Keywords
;         mask   (byte)   2d-array of points to use for calculation. If
;                         not set, all the blushiftarray is used
;         coeff  (double) 3x3 array of fit-coefficients 
;         blueshift_fit_matrix (float) 2-D array containing fitted surface
;         use_squares (byte) flag indicating that squared distances should
;                           be used (approximating a parabolic fit)
; Outputs
;         vector (integer) vector of x and y position of fit-minimum 
; Common Blocks 
; Side Effect  
; Restriction  
; Procedure 
; Modification History 
;         28-Nov-03          K. Janssen, Arcetri
;         22-Dec-03 KJ renamed to find_optical_center, to avoid using
;                      names twice.
;
; Comment
;         Used by ibis_make_gain.pro
;         Seems not to need any variables out of varfile
;-




xsize = (size(blueshift))(1)
ysize = (size(blueshift))(2)

IF NOT keyword_set(mask) THEN mask = replicate(1, xsize, ysize) 
IF NOT keyword_set(use_squares) THEN use_squares = 0 
IF NOT keyword_set(max_degree) OR (use_squares EQ 1) THEN max_degree = 0 


; Parabolic fit, using irregular points

xmatrix = rebin(dindgen(xsize), xsize, ysize)                      
ymatrix = rebin( reform(dindgen(ysize), 1, ysize), xsize, ysize)
zmatrix = blueshift 

index = where(mask EQ 1) 
vectors_n = [reform(xmatrix(index), 1, n_elements(index)) , $
             reform(ymatrix(index), 1, n_elements(index)) , $
             reform(zmatrix(index), 1, n_elements(index))]

vectors = vectors_n
nterms = 2

fit = sfit(vectors, nterms, /irr, kx=coeff, MAX_DEGREE=max_degree)   

IF (max_degree GE 1) THEN BEGIN
    coeff_new = FLTARR(3,3)
    coeff_new = [[coeff(0),coeff(1),coeff(2)],[coeff(3),coeff(4),0.0],[coeff(5),0.0,0.0]]
    coeff     = coeff_new
ENDIF

;blueshift_fit_matrix = dblarr(xsize, ysize)
;blueshift_fit_matrix(vectors_n(0, *), vectors_n(1, *)) = fit
blueshift_fit_matrix = coeff(0,0) + xmatrix*coeff(0,1) + xmatrix^2*coeff(0,2) + $
        		       ymatrix*coeff(1,0) + ymatrix^2*coeff(2,0) + $
        		       xmatrix * ymatrix * coeff(1,1) + xmatrix^2 * ymatrix^2 * coeff(2,2) + $
        		       xmatrix^2 * ymatrix * coeff(1,2) + xmatrix * ymatrix^2 * coeff(2,1)



;; Find center of Blueshift, depending of the orientation of the parabola

IF coeff(0, nterms) LT 0 THEN $
   center_index=where(fit EQ max(fit)) $
ELSE $
   center_index=where(fit EQ min(fit)) 
   
center_x = ROUND(vectors_n(0, center_index))
center_y = ROUND(vectors_n(1, center_index))
center_z = (center_index)

IF use_squares EQ 1 THEN BEGIN
    xmatrix_cent = xmatrix - center_x[0]
    ymatrix_cent = ymatrix - center_y[0]
    vectors = [reform((xmatrix_cent(index))^2, 1, n_elements(index)) , $
               reform((ymatrix_cent(index))^2, 1, n_elements(index)) , $
               reform(zmatrix(index), 1, n_elements(index))]
    nterms = 1
    fit = sfit(vectors, nterms, /irr, kx=coeff, MAX_DEGREE=max_degree)   
    blueshift_fit_matrix = coeff(0,0) + xmatrix_cent^2 * coeff(0,1) + $
                           ymatrix_cent^2 * coeff(1,0) + xmatrix_cent^2 * ymatrix_cent^2 * coeff(1,1)
    IF coeff(0, nterms) LT 0 THEN BEGIN
        center_index=where(fit EQ max(fit))
    ENDIF ELSE BEGIN
        center_index=where(fit EQ min(fit)) 
    ENDELSE
    center_x = ROUND(vectors_n(0, center_index))
    center_y = ROUND(vectors_n(1, center_index))
    center_z = (center_index)
ENDIF

;; Perhaps delete this lines later, they were inserted to check and
;; compare the blueshift...

;save, blueshift_fit_matrix, file='~/ibis/blueshift_data_fit.sav'
;save, blueshift, file='~/ibis/blueshift_data.sav'


return, [center_x, center_y]
END
