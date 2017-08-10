#!/bin/bash

getcompilelist() {
 file="${1}"
 grep -i "^# *include *\"" "${file}" | \
      sed "s/^# *include *\"\([^\"]\+\)\".*/\1/i" | \
      while read line; do
  doecho=true

  # include files of the project itself are usually on ../include/
  if [ ! -e "${line}" ]; then
   if [ -e "../include/${line}" ]; then
    line="../include/${line}"
   fi
  fi
  if [ ! -z "${thingsdone}" ]; then
   for thing in ${thingsdone}; do
    if [ "${line}" == "${thing}" ]; then
     doecho=false
    fi
   done
  fi
  if $dolist; then
   if [ -e "${line}" ]; then
    echo -n "${line} " 
   fi
   thingsdone="${thingsdone} ${line}"
  fi
 done
}

getorderlist() {
 file="${1}"
 grep -i "^ *use " "${file}" | \
      sed "s/^ *use \([^,]\+\)\(\$\|,.*\)/\1/i" | \
      while read line; do 
  dolist=true
  if [ ! -z "${thingsdone}" ]; then
   for thing in ${thingsdone}; do
    if [ "${line}" == "${thing}" ]; then
     dolist=false
    fi
   done
  fi
  if $dolist; then
   echo -n "$(echo ${line} | sed "s/\([^ !]\+\).*/\1/").o " 
   thingsdone="${thingsdone} ${line}"
  fi
 done
}

for file in *.F90; do
 fileobj="$(echo "${file%.F90}.o")"
 makeline="${fileobj}: $file "

 # files here just set the compilation order but does not require the current
 # file to be recompiled.
 orderlist="$(getorderlist "${file}")"

 # files here if changed require the current file to be recompiled
 compilelist="$(getcompilelist "${file}")" 

 # update result if any occurrence is found
 if [ ! -z "${compilelist}" ]; then
  makeline="${makeline}${compilelist}"
 fi
 if [ ! -z "${orderlist}" ]; then
  makeline="${makeline}| ${orderlist}"
 fi
 echo "${makeline}"
done
