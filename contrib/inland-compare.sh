#!/bin/bash

#set -x

###########################
# usage: inland-compare-xxx.sh [all | run_ref | run_comp | diff]
#
# run from a directory which contains the source directories to compare
# any variables here will override the global defaults in inland-compare-functions.sh

# "all" does run_ref, run_comp and diff
# "run_ref compiles and executes the reference run(s)
# "run_comp compiles and executes the compare run(s)
# "diff" computes ans plots differences between the ref and each of the comp runs

arg_command="all"
if [[ $# -gt 0 ]]; then
    arg_command=$1
fi
arg_scriptname=$0

source `dirname $arg_scriptname`/inland-compare-functions.sh

## Now do the actual stuff

func_init `dirname $arg_scriptname`
func_print

#temp vars (for testing script)

do_run_ref=1
do_run_comp=1
do_clean=0
do_cdo=1
basedir_ref=inland-ref1
basedirs_comp="inland-ref2"
compile_compilers=( "gfortran" )
compile_numflags=1
compile_flags_names=( "O2" )
compile_flags_gfortran=( "-g -O2" )
#configure_flags_all="--disable-openmp"
configure_flags_all=""

nruns=("2")
isimvegs=("0" "1")
# use this for restart run
#isimvegs=("1" "1982")

# plotting options
do_plot=1
do_plot1d=1
do_plot2d=0
do_plot2dclim=1
domains=("small1")

# change this for default # of threads to use (default is 2)
numthreads_ref="3"

# add this for openmp runs
#do_openmp=1
#numthreads_comp=( "4" "2")

func_exec $arg_command
