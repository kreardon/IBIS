;-------------------------------------------------------------
;+
; NAME: ibis_obs_modes_list
;       
; PURPOSE: returns a list of the defined IBIS observation modes
;       
; CATEGORY:
; CALLING SEQUENCE:
;       obs_modes = ibis_obs_modes_list()
; INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
; COMMON BLOCKS:
;       
; NOTES:
;       the output array will be updated as the list of defined
;       IBIS observing modes is modified (as manifested in changes
;       in the observing mode identifiers written in the IBIS VEE logs).
; MODIFICATION HISTORY:
;       K. Reardon, May 2007, Initial Implementation
;-
;-------------------------------------------------------------
FUNCTION ibis_obs_modes_list, obs_modes_desc=obs_modes_desc

; the different supported observing mode
obs_modes               = STRARR(12)

obs_modes(0)            = 'Unknown'
obs_modes(1)            = 'Science Observation' 
obs_modes(2)            = 'Imaging Spectral Scan'
; Modes 3, 4, 5 and 6 added in log format version 2,0 (Jan 2004)
obs_modes(3)            = 'Dark Calibration'
obs_modes(4)            = 'Flat Field Calibration'
obs_modes(5)            = 'Other Calibration'
obs_modes(6)            = 'Testing'
; Mode 7 added in Dec 2004
obs_modes(7)            = 'LCVR Polarimeter Calibration'
; Mode 8 added in June 2006
obs_modes(8)            = 'Manual Instrumental Polarization'
; Mode 9 added in May 2007
obs_modes(9)            = 'Target Images'
obs_modes(10)           = 'Grid Images'
; Mode 11 added in May 2009
obs_modes(11)           = 'Automatic Instrumental Polarization'

; When being used for descriptive, human-readable purpose, some of 
; the terms for the observing modes can be slightly modified. We also
; output an array with these modified terms, where necessary.
obs_modes_desc          = obs_modes
obs_modes_desc(2)       = 'Prefilter Scan'
obs_modes_desc(7)       = 'LCVR Calibration'
obs_modes_desc(8)       = 'Manual X-Matrix Calibration'
obs_modes_desc(11)      = 'Automatic X-Matrix Calibration'

RETURN, obs_modes

END
