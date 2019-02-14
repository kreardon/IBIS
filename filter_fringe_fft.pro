FUNCTION filter_fringe_fft, imin, halpha=halpha, ca8542=ca8542, na5896=na5896, k7699=k7699, fe5434=fe5434, fft_power=fft_power, fft_mask=fft_mask

nx = 1000
ny = 1000

mask1_pos = [nx/2.,ny/2.,2]
mask2_pos = [nx/2.,ny/2.,2]
mask5_pos = [nx/2.,ny/2.,2]

IF KEYWORD_SET(halpha) THEN BEGIN

    mask1_pos = [77,161,29]
    mask2_pos = [187,847,25]
    mask5_pos = [nx/2.,ny/2.,3]

ENDIF

IF KEYWORD_SET(ca8542) THEN BEGIN
    ;mask1_pos = [145,90,25]
    mask1_pos = [60,137,29]
    ;mask2_pos = [890,134,25]
    mask2_pos = [151,914,25]
    ;mask5_pos = [88,215,3]
    mask5_pos = [91,778,3]
ENDIF

IF KEYWORD_SET(na5896) THEN BEGIN

    mask1_pos = [112,227,19]
    mask2_pos = [246,900,15]
    mask5_pos = [nx/2.,ny/2.,3]

ENDIF

IF KEYWORD_SET(k7699) THEN BEGIN

    mask1_pos = [92,180,29]
    mask2_pos = [112,220,19]
    mask5_pos = [nx/2.,ny/2.,3]

ENDIF

IF KEYWORD_SET(fe5434) THEN BEGIN

    mask1_pos = [100,212,49]
    mask2_pos = [nx/2.,ny/2.,3]
    mask5_pos = [nx/2.,ny/2.,3]
    ;mask2_pos = [762,182,9]
    ;mask5_pos = [856,351,5]

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


mask1a  = (radist1 GE (mask1_pos[2]*1.05)) AND (radist1 LE (mask1_pos[2]*1.3))
mask1a  = mask1a + ROTATE(mask1a,2)
mask2a  = (radist2 GE (mask2_pos[2]*1.05)) AND (radist2 LE (mask2_pos[2]*1.3))
mask2a  = mask2a + ROTATE(mask2a,2)
mask5a  = (radist5 GE (mask5_pos[2]*1.1)) AND (radist5 LE (mask5_pos[2]*1.5))
mask5a  = mask5a + ROTATE(mask5a,2)

imin_mean = mean(imin)
abs_offset = ABS(MIN(imin - imin_mean)) * 2

imin_fft  = FFT(imin - imin_mean,-1)
imin_fft_cor = imin_fft
imin_fft_abs = ABS(imin_fft)
imin_fft_cor = (imin_fft_cor * ABS(1 - mask1)) + (mask1 * MEDIAN(imin_fft_abs(WHERE(mask1a))))
imin_fft_cor = (imin_fft_cor * ABS(1 - mask2)) + (mask2 * MEDIAN(imin_fft_abs(WHERE(mask2a))))
imin_fft_cor = (imin_fft_cor * ABS(1 - mask5)) + (mask5 * MEDIAN(imin_fft_abs(WHERE(mask5a))))

imin_out = ABS(FFT(imin_fft_cor,1) + imin_mean + abs_offset) - abs_offset

fft_mask  = mask1 + mask2 + mask5
fft_power = imin_fft_abs

RETURN,imin_out

END
