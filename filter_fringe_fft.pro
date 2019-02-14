FUNCTION filter_fringe_fft, imin, halpha=halpha, ca8542=ca8542, na5896=na5896, fft_power=fft_power, fft_mask=fft_mask

nx = 1000
ny = 1000
apod_cutoff = 0.15

mask1_pos = [70,150,0]
mask2_pos = [175,890,0]
mask5_pos = [nx/2.,ny/2.,0]

IF KEYWORD_SET(halpha) THEN BEGIN

;    mask1_pos = [77,161,29]
;    mask2_pos = [187,847,25]
;    mask5_pos = [nx/2.,ny/2.,3]

    mask1_pos = [80,170,29]
    mask2_pos = [200,875,29]
    mask5_pos = [nx/2.,ny/2.,3]

ENDIF

IF KEYWORD_SET(ca8542) THEN BEGIN
    ;mask1_pos = [145,90,25]
    mask1_pos = [63,134,31]
    ;mask2_pos = [890,134,25]
    mask2_pos = [153,914,31]
    ;mask5_pos = [88,215,3]
    mask5_pos = [91,778,3]
ENDIF

IF KEYWORD_SET(na5896) THEN BEGIN
    mask1_pos = [90,188,35]
    mask2_pos = [113,230,15]
    mask5_pos = [nx/2.,ny/2.,3]
ENDIF

    radist1 = radial_distances([1,nx,ny],mask1_pos[0:1])
    mask1   = radist1 LE mask1_pos[2]
    mask1   = mask1 + ROTATE(mask1,2)

    radist2 = radial_distances([1,nx,ny],mask2_pos[0:1])
    mask2   = radist2 LE mask2_pos[2]
    mask2   = mask2 + ROTATE(mask2,2)

    radist5 = radial_distances([1,nx,ny],mask5_pos[0:1])
    mask5   = radist5 LE mask5_pos[2]
    mask5   = mask5 + ROTATE(mask5,2)


mask1a  = (radist1 GE (mask1_pos[2]*1.05)) AND (radist1 LE (mask1_pos[2]*1.2))
mask1a  = mask1a + ROTATE(mask1a,2)
mask2a  = (radist2 GE (mask2_pos[2]*1.05)) AND (radist2 LE (mask2_pos[2]*1.2))
mask2a  = mask2a + ROTATE(mask2a,2)
mask5a  = (radist5 GE (mask5_pos[2]*1.1)) AND (radist5 LE (mask5_pos[2]*1.5))
mask5a  = mask5a + ROTATE(mask5a,2)

imin_mean = mean(imin)
abs_offset = ABS(MIN(imin - imin_mean)) * 2
apodwin = apod(nx,ny,0.01,0.01,2)
;apodwin = INTARR(nx,ny) + 1

imin_fft  = FFT((imin - imin_mean) * apodwin,-1)
imin_fft_cor = imin_fft
imin_fft_abs = ABS(imin_fft)
imin_fft_cor = (imin_fft_cor * ABS(1 - mask1)) + (mask1 * MEDIAN(imin_fft_abs(WHERE(mask1a))))
imin_fft_cor = (imin_fft_cor * ABS(1 - mask2)) + (mask2 * MEDIAN(imin_fft_abs(WHERE(mask2a))))
imin_fft_cor = (imin_fft_cor * ABS(1 - mask5)) + (mask5 * MEDIAN(imin_fft_abs(WHERE(mask5a))))

imin_out = ABS(FFT(imin_fft_cor,1) / (apodwin > apod_cutoff) + imin_mean + abs_offset) - abs_offset

apod_win_low = WHERE(apodwin LE apod_cutoff)
imin_out(apod_win_low) = imin(apod_win_low)

fft_mask  = mask1 + mask2 + mask5
fft_power = imin_fft

RETURN,imin_out

END
