#!/bin/bash
# This script will toggle the model code from double precision to single 
# precision. It only works setting from single to double if the single
# precision variables were previously changed from double to single.
# In other words, it sets real*8 to real*4, and real*4 to real*8, it will not
# turn "real" without precision specification into real*8. This is done to
# facilitate the process of keeping real*4 variables that should keep single 
# precision if the model is running at double precision in general.

# usage: toggleprecision.sh 8to4|4to8

function syntax_error() {
 echo "Usage: ${0} 8to4|4to8"
 exit 1
}

if [ -z "${1}" ]; then
 syntax_error
fi

if [ "${1}" == "8to4" ]; then
 todouble=false
elif [ "${1}" == "4to8" ]; then
 todouble=true
else
 syntax_error
fi

if [ ! -d include -o ! -d src ]; then
 echo "ERROR: You should run this script from the model source's root directory."
 exit 1
fi

if $todouble; then
 vardeclpat="s/\(^\|[^a-z]\)real\*4\(\$\|[^0-9]\)/\1real*8\2/gi"
 nfpat="s/\(\(^\|[^A-Z]\)NF_[A-Z_]*\)REAL/\1DOUBLE/gi"
 castpat="s/\(^[^\!]\+[^a-z\!]\)real(/\1dble(/gi"
 modpat="s/\(^[^\!]\+[^a-z\!]\)amod(/\1dmod(/gi"
else
 vardeclpat="s/\(^\|[^a-z]\)real\*8\(\$\|[^0-9]\)/\1real*4\2/gi"
 nfpat="s/\(\(^\|[^A-Z]\)NF_[A-Z_]*\)DOUBLE/\1REAL/gi"
 castpat="s/\(^[^\!]\+[^a-z\!]\)dble(/\1real(/gi"
 modpat="s/\(^[^\!]\+[^a-z\!]\)dmod(/\1amod(/gi"
fi

cd include
sed "${vardeclpat};${nfpat};${castpat};${modpat}" -i *.h 
cd -
cd src
sed "${vardeclpat};${nfpat};${castpat};${modpat}" -i *.F90 
echo "All done."
