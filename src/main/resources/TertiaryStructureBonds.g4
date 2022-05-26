
/*
 * Grammar for reading files containing simple Tertiary stuctures without sequence
 * 
 * Arc-Annotated Sequence (AAS) Format - - with no sequence
 * Example: (10,20);(0,80);(2,80);(2,45)
 * 
 * @author Luca Tesei
 * 
 */
grammar TertiaryStructureBonds;

@header {
	package it.unicam.cs.bdslab.tersaling;
}

bonds
:
	bond ';' bonds # bondsContinue
	| bond # bondsEnd
;

bond
:
	'(' INDEX ',' INDEX ')'
;


// Tokens


INDEX
:
	[1-9] [0-9]* 
	| '0'
;

LINE_COMMENT
:
	'#' .*? '\r'? '\n' -> skip
; // Match "#" stuff '\n' and skip it

WS
:
	[ \t\r\n]+ -> skip
; // skip spaces, tabs, newlines

