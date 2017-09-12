FUNCTION calc_dest_vect_noflows, image_seq, kernels=kernels, destr_vect_raw=destr_vect_raw

;+
; NAME:
;       calc_dest_vect_noflows
; PURPOSE:
;       calculate the destretching vectors for a time series of images, removing
;       the long term trends (i.e. solar flows) from the measured subfield motions  
; CATEGORY:
; CALLING SEQUENCE:
;       destr_vect = calc_dest_vect_noflows(image_seq)
; INPUTS:
;       image_seq = a time series of (speckle-reconstructed) images.     in
;                   either an image array or an array of filenames.
; KEYWORD PARAMETERS:
;       -- kernel = kernel to be used in destretch calculations = BYTARR(kx,ky)
;       kernels = list of dimensions of kernels for nested destretching
; OUTPUTS:
;       destr_vect =     structure containing the time series of destretch 
;                         vectors with flows removed.        out
;       destr_vect_raw = structure containing the time series of destretch
;                         vectors as measured.               out
; COMMON BLOCKS:
; NOTES:
;       Uses reg_fft.pro algorithm to calculate destretch vectors. 
;       Only able to use a single iteration of the destretching
;
; MODIFICATION HISTORY:
;       Kevin Reardon,  October, 2010
;
;-
;-------------------------------------------------------------

; the size of the smoothing and median range may need to be adjusted for each given 
; data set - the smoothing should probably be done over a number of images 
; corresponding to several minutes of time.
smooth_num        = 21
median_num        = 5

image_seq_size    = SIZE(image_seq, /ST)
IF image_seq_size.TYPE_NAME EQ 'STRING' THEN BEGIN
    file_list     = 1
    numims        = image_seq_size.DIMENSIONS[0]
    im_ref        = load_speckle_im(image_seq[0])
    im_cor        = load_speckle_im(image_seq[1])
    nx            = N_ELEMENTS(im_ref[*,0])
    ny            = N_ELEMENTS(im_ref[0,*])
ENDIF ELSE BEGIN
    file_list     = 0
    numims        = image_seq_size.DIMENSIONS[2]
    nx            = N_ELEMENTS(image_seq[*,0,0])   ; or image_seq_size.DIMENSIONS[0]
    ny            = N_ELEMENTS(image_seq[0,*,0])   ; or image_seq_size.DIMENSIONS[1]
    im_ref        = image_seq[*,*,0]
    im_cor        = image_seq[*,*,1]
ENDELSE
num_spk_files     = numims

; large array to hold all destretched images
; spkim_cor          = FLTARR(nx,ny,numims,/NOZ)

; if provided as keyword, "kernel" should be full 2-D kernel, not the dimension of the kernel
IF NOT KEYWORD_SET(kernel) THEN kernel= BYTARR(ROUND(nx/20.),ROUND(ny/20.))
;test_reg           = reg(im_cor, im_ref, kernel,disp=wl_disp,rdisp=wl_rdisp)
test_reg           = reg_loop(im_cor, im_ref, kernels,disp=wl_disp,rdisp=wl_rdisp)

disp_sz_x          = N_ELEMENTS(wl_disp(0,*,0))
disp_sz_y          = N_ELEMENTS(wl_disp(0,0,*))
spkim_wl_disp      = FLTARR(2, disp_sz_x, disp_sz_y, numims)
spkim_wl_rdisp     = FLTARR(2, disp_sz_x, disp_sz_y, numims)
; make a mask to only use central portion of calculated displacements
; especially suitable if IBIS images with circular mask are used, otherwise
; it may be necessary to allow disp_mask to be optionaly provided through a keyword.
disp_radist        = radial_distances([1,disp_sz_x,disp_sz_y],[disp_sz_x/2.,disp_sz_y/2.])
disp_mask          = disp_radist LE FIX(MAX([disp_sz_x,disp_sz_y])/2. * 0.95)

;************************************************************************
; step through images calculating destretch vectors between successive images
;************************************************************************
FOR im = 1,numims-1 DO BEGIN

    IF (im MOD 25) EQ 0 THEN PRINT,'Image - ' + STRTRIM(im,2)

    IF im EQ 1 THEN BEGIN
        ;spkim_ref                = REBIN(image_seq(*,*,1:7),nx,ny,1)
        spkim_ref                = im_ref
    ENDIF ELSE BEGIN
        spkim_ref                = im_cor
    ENDELSE

    IF file_list THEN im_cor = load_speckle_im(image_seq[im]) ELSE im_cor = image_seq(*,*,im)

;    spkim_dest               = reg(im_cor, spkim_ref, kernel, disp=wl_disp, rdisp=wl_rdisp)
    spkim_dest               = reg_loop(im_cor, spkim_ref, kernels, disp=wl_disp, rdisp=wl_rdisp)
    spkim_wl_disp(*,*,*,im)  = wl_disp
    spkim_wl_rdisp(*,*,*,im) = wl_rdisp
ENDFOR

spkim_wl_disp(*,*,*,0)   = wl_rdisp
spkim_wl_rdisp(*,*,*,0)  = wl_rdisp

; create structure of measured destretch vectors
destr_vect_raw = CREATE_STRUCT('disp', spkim_wl_disp, $
                               'rdisp', spkim_wl_rdisp, $
                               'kernels', kernels)

;************************************************************************
; calculate the bulk shifts, the shifts averaged over the whole field
; due to rigid motion of the image. Most of this should have already been
; corrected before destretching
;************************************************************************

; this is the vector magnitude, after removing the local values of the reference positions
shifts_all              = spkim_wl_disp - spkim_wl_rdisp
shifts_bulk             = FLTARR(2,numims)
shifts_bulk_sum         = FLTARR(2,numims)
shifts_bulk_cor         = FLTARR(2, disp_sz_x, disp_sz_y, numims)
shifts_cor_sum          = FLTARR(2, disp_sz_x, disp_sz_y, numims)
shifts_cor_sum(*,*,*,0) = shifts_bulk_cor(*,*,*,0)
step_size               = 1

FOR i=1,numims-1 DO BEGIN
    ; these are the average of the destretch displacements in the two directions
    ibis_area_statistics,shifts_all(0,*,*,i),disp_mask,median=median_timex
    ibis_area_statistics,shifts_all(1,*,*,i),disp_mask,median=median_timey
    shifts_bulk(*,i)  = [median_timex, median_timey]
    
    ; remove these bulk shifts from the destretch vectors
    shifts_bulk_cor(0,*,*,i) = shifts_all(0,*,*,i) - shifts_bulk(0,i)
    shifts_bulk_cor(1,*,*,i) = shifts_all(1,*,*,i) - shifts_bulk(1,i)

    ; and calculate the cumulative displacement at each vector position
    IF i GE 1 THEN $
        shifts_cor_sum(*,*,*,i)  = shifts_cor_sum(*,*,*,i-1) + shifts_bulk_cor(*,*,*,i)
        shifts_bulk_sum(*,i)     = shifts_bulk_sum(*,i-1) + shifts_bulk(*,i)
ENDFOR

spkim_wl_disp_cor = FLTARR(2, disp_sz_x, disp_sz_y, numims)
spkim_wl_disp_der = FLTARR(2, disp_sz_x, disp_sz_y, numims)
spkim_wl_disp_use = FLTARR(2, disp_sz_x, disp_sz_y, numims)
spkim_wl_disp_fit = FLTARR(2, disp_sz_x, disp_sz_y, numims)
spkim_wl_disp_med = FLTARR(2, disp_sz_x, disp_sz_y, 1)

;************************************************************************
; For each reference position in the grid of destretch vectors, calculate
; a running mean of the cumulative displacements (a polynomial fit can also
; be used, but may be less robust). 
; calculate the changes in this running mean between time steps and subtract
; this from the measured shifts in order to remove these smooth variations 
; in the destretch vectors.
;************************************************************************
FOR xx=0,disp_sz_x-1 DO BEGIN
    FOR yy=0,disp_sz_y-1 DO BEGIN
        FOR dir=0,1 DO BEGIN
             disp_seq = REFORM(shifts_cor_sum(dir,xx,yy,*))

             spkim_wl_disp_fit(dir,xx,yy,*)      = SMOOTH(MEDIAN(disp_seq,median_num),smooth_num,/EDGE)
             spkim_wl_disp_use(dir,xx,yy,*)      = shifts_cor_sum(dir,xx,yy,*) - spkim_wl_disp_fit(dir,xx,yy,*)

             ;spkim_wl_disp_der(dir,xx,yy,1:*)    = first_der(REFORM(spkim_wl_disp_fit(dir,xx,yy,*)))
             ;spkim_wl_disp_use(dir,xx,yy,*)      = shifts_bulk_cor(dir,xx,yy,*) - spkim_wl_disp_der(dir,xx,yy,*)
    
         ENDFOR
    ENDFOR
ENDFOR 

; add the bulk displacements back into the destretch vector magnitudes
; this should be fine as long as these bulk displacements are small
FOR i=0,numims-1 DO BEGIN
    ; add the bulk shifts back into the destretch vectors
    spkim_wl_disp_use(0,*,*,i) = spkim_wl_disp_use(0,*,*,i) + shifts_bulk_sum(0,i)
    spkim_wl_disp_use(1,*,*,i) = spkim_wl_disp_use(1,*,*,i) + shifts_bulk_sum(1,i)
ENDFOR

; add the local values of the reference positions to the destretch vector magnitudes, as 
; needed by the "reg_fft" program
spkim_wl_disp_use = spkim_wl_disp_use + spkim_wl_rdisp


destr_vect = CREATE_STRUCT('disp',            spkim_wl_disp_use, $
                           'rdisp',           spkim_wl_rdisp, $
                           'disp_sum',        shifts_cor_sum, $
                           'kernels',         kernels, $
                           'shifts_bulk',     shifts_bulk,$
                           'shifts_bulk_sum', shifts_bulk_sum )

RETURN,destr_vect

END
