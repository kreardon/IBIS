
;+
; NAME:
;        get_closest
; PURPOSE:
;        for each element of the reference vector, search in vector vec for its
;        closest value.  Returns a vector of indices.
; CATEGORY:
; CALLING SEQUENCE:
;        index = get_closest(vec,ref)
; INPUTS:
;        vec = vector to be searched.
;        ref = reference vector
; OUTPUTS:
;        index = index pointing to vec for each element of ref.
; KEYWORDS:
;        out_of_range = flags points where the reference vector is outside of range
;                       of search vector. 
;                       1 if reference vector value is greater than maximum value of search vector
;                      -1 if reference vector value is less than minimum value of search vector
;                       0 if reference vector value is within limits of search vector
;                       
; COMMON BLOCKS:
;        None.
; SIDE EFFECTS:
; RESTRICTIONS:
;        If more than one index in the search vector matches an element of 
;            reference vector, then only the first/lowest matching index 
;            will be returned.
; PROCEDURE:
; MODIFICATION HISTORY:
;        Jean-Pierre Wuelser, many years ago in 1992-ish
;        KPR - 2019, cleaned up, took out date processing options 
;-

function get_closest, vec0, ref0, out_of_range=out_of_range

;vec0 = reform(array(vec0))
;ref0 = reform(array(ref0))
vec0 = reform(vec0)
ref0 = reform(ref0)
svec = size(vec0)
sref = size(ref0)

index        = long(ref0)              ; array with same size and shape as ref0, but type long
out_of_range = FIX(ref0) * 0          ; array with same size and shape as ref0, but type integer
vec0_min     = MIN(vec0, max=vec0_max)

for i=0LL,n_elements(ref0)-1 do begin
    IF ref0(i) LT vec0_min THEN out_of_range[i] = -1
    IF ref0(i) GT vec0_max THEN out_of_range[i] = 1
    
    delt = abs(vec0 - ref0[i])
    dmin = where(delt eq min(delt))
    index[i] = dmin[0]
endfor

return,index

end

