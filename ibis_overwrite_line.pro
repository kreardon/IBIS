;-------------------------------------------------------------
;+
; NAME:
;       ibis_overwrite_line
; PURPOSE:
;       a program to "backspace" a certain number of characters on the
;       current line and print the input string on that same line
; CATEGORY:
; CALLING SEQUENCE:
;       ibis_overwrite_line, overwrite_length, overwrite_string
; INPUTS:
;       overwrite_length - the number of characters to cancel
;       overwrite_string - the new string to print on the line
; KEYWORD PARAMETERS:    
;       new_overwrite_length - the length of the input string
; OUTPUTS:
; COMMON BLOCKS:
;       
; NOTES:
;       this program can be called from within a loop using the following
;       calling sequence, where the value of the input parameter 
;       "output_text_length" is the length of the string previously written
;       to the terminal. Upon completion, the value of "output_text_length"
;       becomes the length of the string written by the program 
;       (i.e output_text). The subsequent call to ibis_overwrite_line as the
;       loop is executed again then knows the length of the previously written
;       string to cancel before writing the new text.
;
;       ibis_overwrite_line, output_text_length, output_text, $
;                            new_overwrite_length=output_text_length
;
;       If output_text_length is less than 1, then no characters are cancelled
;       and the desired text is simply written to the terminal.
;
;       It is unknown whether this program works on non-unix systems.
;
; MODIFICATION HISTORY:
;       18 Nov, 2003 - KPR - initial implementation - copyied backspace 
;                            character ('8B') from one of Katja's programs
;-
;-------------------------------------------------------------
PRO ibis_overwrite_line, overwrite_length, overwrite_string, $
                         new_overwrite_length = new_overwrite_length

backspace               = 8B

overwrite_length_size  = SIZE(overwrite_length, /STRUCTURE)
overwrite_length_type  = overwrite_length_size.TYPE
overwrite_length_tname = overwrite_length_size.TYPE_NAME
overwrite_length_nelem = overwrite_length_size.N_ELEMENTS

; Make sure 'overwrite_length' is some type of number with at least one element
IF (overwrite_length_nelem GE 1) AND $
   (((overwrite_length_type GE 1)  AND (overwrite_length_type LE  6)) OR $
    ((overwrite_length_type GE 12) AND (overwrite_length_type LE 15))) THEN BEGIN

    IF (overwrite_length GE 1) THEN BEGIN
        IF overwrite_length_nelem EQ 1 THEN BEGIN
            PRINT,FORMAT = '(a,$)', STRING(REPLICATE(backspace,overwrite_length))
        ENDIF ELSE BEGIN
            PRINT,FORMAT = '(a,$)', STRING(REPLICATE(backspace,overwrite_length(0)))
	ENDELSE
    ENDIF
; If it isn't print an error message
ENDIF ELSE BEGIN
    PRINT, 'Error! Ibis_Overwrite_Line expects first parameter to be a number.'
    PRINT, 'Instead a ' + overwrite_length_tname + ' was found.'
ENDELSE

IF KEYWORD_SET(overwrite_string) GT 0 THEN BEGIN
    new_overwrite_length = STRLEN(overwrite_string)
    PRINT,FORMAT='(a,$)',overwrite_string
ENDIF ELSE $
    new_overwrite_length = 0

END
