#!/bin/sh
# Linux users have to change $8 to $9
gawk '
BEGIN 	{ print "File\tOwner" } 
		{ print $9, "\t", $3}	
END   	{ print " - DONE -" } 
'
