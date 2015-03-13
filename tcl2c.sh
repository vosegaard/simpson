#!/bin/sh
# convert a wish script to an initialized C string
# Modified 10/1/93 by sau@bellcore.com
# Modified 10/2/94 by jn@berlin.snafu.de

case $# in 
    0)	echo "usage: $0 c_variable_name file1 ..." 1>&2 ; exit 1 ;;
esac

var="$1"
shift

echo '/* This file has been generated from the following Tcl source file(s):'
case $# in
    0)	echo ' *    <stdin>' ;;
    *)	for i in $* ; do
	    echo ' *    '$i
	done ;;
esac
echo ' * on '"`date` by $USER"
echo ' */'
echo ""
echo 'char '$var'[] = '

# Transformations applied to each line of the file:
#       eliminate lines that are only comments
#       Convert preexisting backslash into double backslash.
#       Precede preexisting double quote (") with backslash.
#	Add double quote to the beginning of the line.
#       Add backslash, n, double quote to the end of the line.
    
cat $* |
sed -e '/^[ 	]*#/d'	\
    -e 's/\\/\\\\/g'	\
    -e 's/"/\\\"/g'	\
    -e 's/$/\\n/'	\
    -e 's/^\(.*\)$/    "\1"/'
# Final line of transformed file:

echo '    ;
/* EOF */'

# We don't remove blank lines anymore, because this changes embedded
# text strings.

# EOF
