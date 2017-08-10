#!/bin/bash

#set -x
set -e

###########################
# functions

# func_init - create default variables and environment, run before overriding any vars.
# pass it the directory of this script, so we can get the path of other scripts
function func_init()
{

# check if we are running in CRAY environment
if [[ "$CRAY_CPU_TARGET" == "" ]]; then env_cray=0 ; else env_cray=1; fi

if [[ "$env_cray" == "1" ]]; then echo "running in CRAY environment"; fi

#dir/file vars
datadir=`pwd`/inland-data
basedir_prefix=`pwd`
basedir_ref=inland-ref
basedirs_comp="inland-ref2 inland"
paylogicsubdir="clim" #subdir to exec e.g. daily, monthly, clim, etc.
#FIXME if not absolute path
scriptdir=$1
plot1d=`pwd`/$scriptdir/inland-compare-plot1D.py
plot2d=`pwd`/$scriptdir/inland-compare-plot2D.py
ifile_compar=`pwd`/$scriptdir/inland-compare_inland_compar.h
ifile_infile=`pwd`/$scriptdir/inland-compare_inland-grid.infile
vars="aet anpptot awc co2mic co2root co2soi drainage fl fu neetot nmintot npptot rootbio srunoff totalit totanlit totbiol totbiou totcmic totcsoi totlail totlaiu totnsoi totrlit totrnlit trunoff tsoi vegtype0 vwc wisoi wsoi"

#openmp
numthreads_ref="2"
numthreads_comp=( "4" "2" "1" )
#export OMP_NUM_THREADS=1
unset INLAND_RANDOMVAL


#compilers
if [[ "$env_cray" == "1" ]]; then
source /opt/modules/default/etc/modules.sh 
compile_env=( "PrgEnv-gnu" )
compile_envnum=( "8" )
compile_compilers=( "gfortran" )
else
compile_compilers=( "gfortran" "ifort" )
fi
compile_numflags=2
compile_flags_names=( "O2" "O0")
compile_flags_gfortran=( "-g -O2" "-g -O0" )
compile_flags_intel=( "-g -O2" "-g -O0" )
#configure_flags_all="--disable-openmp"
configure_flags_all=""


#conf
do_clean=1
do_compile=1
do_configure=1
do_run_ref=1
do_run_comp=1
do_run=1
do_cdo=1
do_cdo_diff=1
do_plot=1
do_plot2d=0
do_plot1d=0
do_plot2dclim=1
do_plot_pft=0
do_plot_ref=1 #plot ref+comp (not dif)
#do_plot_comp=1
do_plot_openmp=1 #plot openmp diffs
do_openmp=0
do_run_parallel=0 #should we do all the runs in parallel ? watch out if many runs!!!
do_qsub=0 #should we submit the batch job or just run interactively (only in cray env.)
nice_str="/usr/bin/time -f \"%Uuser %Ssystem %Eelapsed %PCPU (%Kavg %Xtext+%Ddata %Mmax)k\" nice -n 10 ionice -c 3"
#nice_str="valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes" 
#nice_str="valgrind -v --tool=memcheck --leak-check=full --track-origins=yes"

domains=("small1" "small2" "amz" "toc" "sam" "global")
# define domain parms: snorth ssouth swest seast LPT idailyout
# note: idailyout recommended for short runs and/or small domains only
domain_params_small1=( "-7.75" "-9.25" "-50.25" "-48.75" "16" "1" )
domain_params_small2=( "-14.25" "-14.74" "-57.75" "-57.25" "4" "1")
domain_params_amz=( "5" "-15" "-75" "-55" "1680" "0" )
domain_params_toc=( "-3" "-12" "-56" "-42" "550" "0" )
domain_params_sam=( "10" "-30" "-90" "-30" "5080" "0" )
domain_params_global=( "90" "-90" "-180" "180" "58920" "0" )

nruns=("2" "5")
isimvegs=("0" "1")
iyear0="1981"
iyrdaily="9999"
iyrmon="9999"
}

# func_print - print run configuration (TBD)
function func_print()
{
    echo "run configuration=================="
    echo "script parameters: $*"
    echo command: $command
}


function func_restore_conf()
{
    make_dir=$1
    # restore domain conf files in source
    if [[ "$is_sage" == "0" ]]; then
        pushd .
        cd $make_dir
        echo "restoring conf files in $make_dir"
	    ofile=include/inland_compar.h
	    cp -f ${ofile}.comp ${ofile}
	    ofile=./data/offline/grid/conf/inland-grid.infile
	    cp -f ${ofile}.comp ${ofile}
        popd
    fi        
}

# func_run - compile and run all configurations
function func_run()
{
    basedir=$1
    if [[ "${basedir:0:4}" == "sage"  ]] ; then is_sage=1 ; else is_sage=0; fi
    echo "do_run, basedir = "$basedir" subdir= "$subdir" is_sage="$is_sage
    cd ${basedir_prefix}/$basedir
    
    #loop domain
    for (( i_domains = 0 ; i_domains < ${#domains[*]} ; i_domains++ )) ; do      
        domain=${domains[$i_domains]}
        #domain_param_domain=domain_params_${domain}[0]
	    domain_param_snorth=domain_params_${domain}[0];
	    domain_param_ssouth=domain_params_${domain}[1]
	    domain_param_swest=domain_params_${domain}[2]
	    domain_param_seast=domain_params_${domain}[3]
	    domain_param_lpt=domain_params_${domain}[4]
	    domain_param_idailyout=domain_params_${domain}[5]
        
	    echo ""
        echo ++++++++++ domain: $domain
	    echo params: ${!domain_param_snorth} ${!domain_param_ssouth} ${!domain_param_swest} ${!domain_param_seast} ${!domain_param_lpt} ${!domain_param_idailyout}
	    #echo params: ${!domain_params_domain} ${!domain_param_lpt} ${!domain_param_idailyout}
        
        # copy domain conf files
        make_dir=${basedir_prefix}/$basedir
        cd $make_dir

        if [[ "$is_sage" == "0" ]]; then

        #ifile=${basedir_prefix}/compare/conf/inland_compar.h.$domain
        #ifile=${basedir_prefix}/${ifile_compar}
        ifile=${ifile_compar}
        ofile=include/inland_compar.h
        echo cp -f ${ofile} ${ofile}.comp
        cp -f ${ofile} ${ofile}.comp
        echo cp -f $ifile $ofile
        cp -f $ifile $ofile
        # this sed expression also modified the LPT value for 1D config... 
        sed -i -e"s/[ ]*\#define LPT [0-9]*/\#define LPT ${!domain_param_lpt}/" $ofile
        echo ""
	    echo +++: npoi:
	    grep "\#define LPT" $ofile | head -1

        else
            ifile=compar.h
            ifile_compar=compar.h
            ofile=compar.h
            sed -i -e"s/npoi = [0-9]*, ! change here/npoi = ${!domain_param_lpt}, ! change here/" $ofile
	        echo +++: npoi:
	        grep "npoi = " $ofile | grep "change here" | head -1
        fi

        if [[ "$is_sage" == "0" ]]; then
            infile_tmp=${basedir_prefix}/$basedir/contrib/inland-compare_inland-grid.infile
            if [ -f $infile_tmp ]; then
	            echo "using infile $infile_tmp"
	            ifile_infile=$infile_tmp
            fi
            echo "using infile "$ifile_infile
            #ifile=${basedir_prefix}/compare/conf/inland-grid.infile.${domain}
            #ifile=${basedir_prefix}/${ifile_infile}
            ifile=${ifile_infile}
            infile=./data/offline/grid/conf/inland-grid.infile
            echo cp -f ${ofile} ${ofile}.comp
            cp -f ${ofile} ${ofile}.comp
            #echo cp -f $ifile $infile
            #cp -f $ifile $infile
        else
            ifile=ibis.infile
            ifile_infile=ibis.infile
	        infile=ibis.infile
        fi

        #sed -i -e"s/[+-]*[0-9]*[.,]*[0-9]*[ ]*! snorth/"${!domain_param_snorth}"         ! snorth/" $infile
        #sed -i -e"s/[+-]*[0-9]*[.,]*[0-9]*[ ]*! ssouth/"${!domain_param_ssouth}"         ! ssouth/" $infile
        #sed -i -e"s/[+-]*[0-9]*[.,]*[0-9]*[ ]*! swest/"${!domain_param_swest}"         ! swest/" $infile
        #sed -i -e"s/[+-]*[0-9]*[.,]*[0-9]*[ ]*! seast/"${!domain_param_seast}"         ! seast/" $infile
        #sed -i -e"s/[+-]*[0-9]*[.,]*[0-9]*[ ]*! seast/"${!domain_param_seast}"         ! seast/" $infile
        #sed -i -e"s/[0-9][ ]*! idailyout/"${!domain_param_idailyout}"         ! idailyout/" $infile
        #sed -i -e"s/[0-9]*[ ]*! iyear0/"${iyear0}"         ! iyear0/" $infile
        #sed -i -e"s/[0-9]*[ ]*! iyrdaily/"${iyrdaily}"         ! iyrdaily/" $infile
        #sed -i -e"s/[0-9]*[ ]*! iyrmon/"${iyrmon}"         ! iyrmon/" $infile
        
        echo "! this is inland grid configuration file (namelist)" > $infile
        echo "! this file generated by inland-compare-functions.sh" >> $infile
        echo "&INLAND_GRID" >> $infile
        #echo "domain      = '${!domain_param_domain} ', " >> $infile
        echo "snorth      = ${!domain_param_snorth} , " >> $infile
        echo "ssouth      = ${!domain_param_ssouth} , " >> $infile
        echo "swest       = ${!domain_param_swest} , " >> $infile
        echo "seast       = ${!domain_param_seast} , " >> $infile

        #echo ""
        #echo +++: $infile:
        #echo "--"
        #cat $infile
        #echo "--"

	    echo using base conf files ${ifile_compar} ${ifile_infile}

    #loop compilers
    for (( i_compilers = 0 ; i_compilers < ${#compile_compilers[*]} ; i_compilers++ )) ; do      
        compiler=${compile_compilers[$i_compilers]}
        compiler_dir=$compiler

        # prepare cray build env.
        # this doesn't work reliably (July 6)...
        #if [[ "$env_cray" == "1" ]]; then
        #    echo ${compile_envnum[$i_compiler]} | source /usr/bin/development_config > /dev/null 2> /dev/null 
	    #    module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-pathscale PrgEnv-cray netcdf
        #    module load netcdf
        #    module list
        #fi
        
    # loop compiler flags
    for (( i_flags = 0 ; i_flags < $compile_numflags ; i_flags++ )) ; do      
        compile_flag=compile_flags_${compiler}[i_flags]; 
        echo == compiler: $compiler  flags: ${!compile_flag}       
        
        # compile
        if [[ "$do_compile" == "1" ]]; then
            if [[ "$env_cray" == "1" ]]; then compiler1="ftn"; else compiler1=$compiler; fi
            cd $make_dir
            if [[ "$is_sage" == "0" ]]; then
                if [[ "$do_configure" == "1" ]]; then
                    ./configure FC=$compiler1 FCFLAGS="${!compile_flag}" $configure_flags_all > /dev/null
                    if [ "$?" -ne 0 ]; then echo "ERROR during configure"; exit $?; fi 
                fi
            fi
            if [[ "$do_clean" == "1" ]]; then make clean > /dev/null ; fi
            nice make -j6 > /dev/null
            if [ "$?" -ne 0 ]; then echo "ERROR during make"; exit $?; fi 
            #make dev-symlinks
            echo "+++done compiling"
        fi

        #loop nrun
        for (( i_nruns = 0 ; i_nruns < ${#nruns[*]} ; i_nruns++ )) ; do      
            nrun=${nruns[$i_nruns]}
            echo ++++++++++ nrun: $nrun

        # loop isimveg flag
        for (( i_isimvegs = 0 ; i_isimvegs < ${#isimvegs[*]} ; i_isimvegs++ )) ; do      
            isimveg=${isimvegs[$i_isimvegs]}
            echo +++: isimveg = $isimveg

        # loop numthreads
        for (( i_numthreads = 0 ; i_numthreads < ${#numthreads[*]} ; i_numthreads++ )) ; do      
            if [ "$do_openmp" -eq 1 ]; then
                numthread=${numthreads[$i_numthreads]}
            else
                numthread=$numthreads_ref
                #export OMP_NUM_THREADS=1
            fi          
	        if [ "$numthread" -eq 0 ]; then numthread=1; fi 
            export OMP_NUM_THREADS=$numthread
            echo +++: numthreads = $OMP_NUM_THREADS
            numthread_str=""
            if [ "$do_openmp" -eq 1 -a "$numthread" -ne 0 ] ; then numthread_str=_omp_$numthread ; fi

            #work_dir=${basedir_prefix}/compare/$basedir/$compiler/$domain/${nrun}yr/isimveg$isimveg/flags_${compile_flags_names[$i_flags]}
            work_dir=${basedir_prefix}/compare/$basedir/$subdir/$domain/${nrun}yr/isimveg$isimveg/${compiler}_${compile_flags_names[$i_flags]}${numthread_str}

            # make install
            echo "installing files to $work_dir"
            cd $make_dir

            if [[ "$is_sage" == "0" ]]; then
                ./configure FC=$compiler1 FCFLAGS="${!compile_flag}" $configure_flags_all --prefix=${work_dir} > /dev/null
                make install > /dev/null
                if [ "$?" -ne 0 ]; then echo "ERROR during make install"; exit $?; fi 
            else
                mkdir -p ${work_dir}
                cp -r . ${work_dir}
            fi

            echo "+++done compiling"

            # prepare input files
            cd $work_dir
            if [[ "$is_sage" == "0" ]]; then
                if [ ! -d share/inland/ ] ; then ln -s share/doc/inland/ data ; else ln -s share/inland/ data ; fi
                ln -s data/offline/grid/conf
                ln -s data/offline/grid/params
                rm input
                ln -s $datadir/input 
                #export INLAND_DATADIR=$datadir/input
                mkdir output
                #export INLAND_OUTDIR=$work_dir/output2
                #mkdir $INLAND_OUTDIR
                infile1=data/offline/grid/conf/inland-grid.infile1
                infile=data/offline/grid/conf/inland-grid.infile
            else
                infile=ibis.infile
                rm input
                ln -s $datadir/input 
            fi
            cp -f $infile $infile1
            #sed -i -e"s/[0-9][ ]*! nrun/${nrun}         ! nrun/" $infile
            echo "nrun        = $nrun , " >> $infile
            if [ "$isimveg" -gt 2 ]; then isimveg0=1 ; else isimveg0=$isimveg; fi
            #sed -i -e"s/[0-9][ ]*! isimveg/${isimveg0}         ! isimveg/" $infile
            echo "isimveg     = $isimveg0 , " >> $infile
            echo "!" >> $infile
            echo " / " >> $infile
            echo '' >> $infile

	        echo ""
            echo +++: $infile:
            echo "--"
            cat $infile
            echo "--"

            # remove old files from workdir
            logfile="log.txt"
            rm -f *.nc output/*.nc output/*.dat $logfile

            if [[ "$do_run" == "1" ]]; then
                # execute
                #$nice_str ./bin/inland-grid > $logfile 2>&1
		        if [ "$nrun" -gt 5 -a "$env_cray" -eq 1 -a "$do_qsub" -eq 1 ]; then 
                    echo "qsub!!!";
		        else                 
                    if [[ "$is_sage" == "0" ]]; then
                        echo +++ running $nice_str ./inland-grid
                        echo running dir is `pwd` output to $logfile
		                if [ "$do_run_parallel" = "1" -a "$isimveg" -lt 2  ]; then
			                eval $nice_str ./bin/inland-grid > $logfile & 2>&1
		                else
			                eval $nice_str ./bin/inland-grid > $logfile 2>&1
		                fi
                    else
                        echo +++ running $nice_str ./ibis
			            eval $nice_str ./ibis > $logfile 2>&1

                    fi
                fi
                if [ "$?" -ne 0 ]; then 
                    echo ""; echo "============"; echo "ERROR during execution, command returned $?"; 
                    cat $logfile ; 
                    # restore domain conf files in source
                    func_restore_conf $make_dir
                    exit $?; 
                fi 
                echo "+++done executing $work_dir"
                tail -10 $logfile

                # restart run
                if [ "$isimveg" -gt 2 ]; then 
                    cp -f $infile1 $infile
                    #sed -i -e"s/[0-9][ ]*! irestart/1         ! irestart/" $infile
                    echo "irestart    = 1 , " >> $infile
                    #sed -i -e"s/[0-9][0-9][0-9][0-9][ ]*! iyrrestart/${isimveg}      ! iyrrestart/" $infile
                    echo "iyrrestart  = ${isimveg} , " >> $infile
                    # TODO calculate nrun when not using last year
                    echo "nrun        = 1 , " >> $infile
                    echo "!" >> $infile
                    echo " / " >> $infile
                    echo '' >> $infile
	                echo ""
	                echo "restart run..."
	                echo ""
                    echo +++: $infile:
                    echo ""
                    cat $infile
                    echo ""

		            if [ "$nrun" -gt 5 -a "$env_cray" -eq 1 -a "$do_qsub" -eq 1 ]; then 
                        echo "qsub!!!";
		            else                 
                        if [[ "$is_sage" == "0" ]]; then
                            echo +++ running $nice_str ./inland-grid
			                eval $nice_str ./bin/inland-grid > $logfile 2>&1
		                else
                            echo +++ running $nice_str ./ibis
			                eval $nice_str ./ibis > $logfile #2>&1
                        fi
                    fi
                    if [ "$?" -ne 0 ]; then 
                        echo ""; echo "============"; echo "ERROR during execution, command returned $?"; 
                        cat $logfile ; 
                        # restore domain conf files in source
                        func_restore_conf $make_dir
                        exit $?; 
                    fi 

                    echo "+++done executing $work_dir"
                    tail -10 $logfile
                fi

                #cat $logfile
            fi

        done #numthreads
        export OMP_NUM_THREADS=1
        done #isimveg
        done #nruns
    done #compiler flags
    done #compilers
    #  restore domain conf files in source
    func_restore_conf $make_dir
    done #domains
}

# func_diff - print difference summary and make figures
function func_diff()
{
    echo func_diff $*
    filename=$4
    filename_ts=`basename $filename .nc`_ts.nc
    filename_tsnorm=`basename $filename .nc`_tsnorm.nc
    filename_clim=`basename $filename .nc`_clim.nc
    filename_climnorm=`basename $filename .nc`_climnorm.nc
#    filename_norm=`basename $filename .nc`_norm.nc
    _do_plot1d=$5
    _do_plot2d=$6
    _do_plot_ref=$7
    _do_plot_comp=$8
#    if [ "$7" -gt 0 ] ; then _do_plot_ref=0 ; else  _do_plot_ref=1 ; fi

    file1=$1/$filename
    file1_ts=$1/${filename_ts}
    file1_clim=$1/${filename_clim}
    file2=$2/${filename}
    file2_ts=$2/${filename_ts}
    file2_clim=$2/${filename_clim}
    ofile=$3/$4
    ofile_ts=$3/${filename_ts}
    ofile_tsnorm=$3/${filename_tsnorm}
    ofile_clim=$3/${filename_clim}
    ofile_climnorm=$3/${filename_climnorm}
#    ofile_norm=$3/${filename_norm}
    if [ ! -f $file1 ] ; then echo "file $file1 missing" ; return 1 ; fi
    if [ ! -f $file2 ] ; then echo "file $file2 missing" ; return 1 ; fi

    # basic cdo diff
    if [ "$do_cdo_diff" = 1  -o "$command" = "all" ] ; then
        tmpoutput=`cdo diffn $file1 $file2 2>/dev/null`  
        echo "--------------- cdo diff:"
        echo "$tmpoutput" | grep -v "0.0000" | grep -v "Date  Time    Name"
        echo "$tmpoutput" | grep "records differ"
        echo "------------------------"
    fi
    
    if [ "$do_cdo" = 1  -o "$command" = "all" ] ; then

    # generate diff file
    mkdir -p $3
    cdo -s -r -O sub $file2 $file1 $ofile

    # generate ts+clim files
    if [ "$_do_plot_ref" = 1 ] ; then
        cdo -s -r -O fldmean $file1 $file1_ts
        cdo -s -r -O ymonmean $file1 $file1_clim
    fi
    if [ "$_do_plot_comp" = 1 ] ; then
	cdo -s -r -O fldmean $file2 $file2_ts
	cdo -s -r -O ymonmean $file2 $file2_clim
    fi
    cdo -s -r -O fldmean $ofile $ofile_ts
    cdo -s -r -O ymonmean $ofile $ofile_clim
    cdo -s -r -O mulc,100 -div -sub $file2_ts $file1_ts $file1_ts $ofile_tsnorm
    cdo -s -r -O mulc,100 -div -sub $file2_clim $file1_clim $file1_clim $ofile_climnorm
#    cdo -s -r -O mulc,100 -div -sub $file2 $file1 $file1 $ofile_norm

    fi #do_cdo

    # plot files
    if [ "$do_plot" = 1 ] ; then
        echo "plotting, cwd="`pwd`
        mkdir -p $1/plot2d $2/plot2d $3/plot2d
        mkdir -p $1/plotts $2/plotts $3/plotts
        mkdir -p $1/plot2dclim $2/plot2dclim $3/plot2dclim
        mkdir -p $3/plottsnorm $3/plot2dclimnorm
        mkdir -p $3/montage $3/montage2d
	
        # 1d plots
        if [ "$do_plot1d" = 1 -a "$_do_plot1d" = 1 ] ; then
            $plot1d $ofile_ts plotts $do_plot_pft &
            $plot1d $ofile_tsnorm plottsnorm $do_plot_pft 1 &
            if [ "$_do_plot_ref" = 1 ] ; then
                $plot1d $file1_ts plotts $do_plot_pft
            fi
	        if [ "$_do_plot_comp" = 1 ] ; then
		        $plot1d $file2_ts plotts $do_plot_pft
	        fi
            wait
        fi

        #2d plots
        if [ "$do_plot2dclim" = 1 ] ; then
            echo $plot2d $ofile_clim plot2dclim $do_plot_pft
            $plot2d $ofile_clim plot2dclim $do_plot_pft &
            echo $plot2d $ofile_climnorm plot2dclimnorm $do_plot_pft 1
            $plot2d $ofile_climnorm plot2dclimnorm $do_plot_pft 1 &
            if [ "$_do_plot_ref" = 1 ] ; then
                $plot2d $file1_clim plot2dclim $do_plot_pft
            fi
	        if [ "$_do_plot_comp" = 1 ] ; then
		        echo $plot2d $file2_clim plot2dclim $do_plot_pft
		        $plot2d $file2_clim plot2dclim $do_plot_pft
	        fi
            wait
        fi
        echo plot2d: $do_plot2d / $do_plot2d
        if [ "$do_plot2d" = 1 -a  "$_do_plot2d" = 1 ] ; then
            $plot2d $ofile plot2d &
            if [ "$_do_plot_ref" = 1 ] ; then
                $plot2d $file1 plot2d
            fi
 	        if [ "$_do_plot_comp" = 1 ] ; then
		        $plot2d $file2 plot2d
	        fi
#                $plot2d $ofile_norm plot2dnorm $do_plot_pft 1
            wait
        fi
	    wait
	    
        #montage of plots
        
        echo "montage2d"
        if [ "$do_plot2d" = 1 -a  "$_do_plot2d" = 1 -a  "$_do_plot_comp" = 1 ] ; then
	        for var in $vars; do
	            echo -n $var" "
                for f in `ls $3/plot2d/inland-yearly_${var}*.png` ; do
                    f_base=`basename $f`
                    montage -geometry +1+20 -tile 2x2 -pointsize 20  -label $basedir_ref $1/plot2d/$f_base -label diff $3/plot2d/$f_base  -label $basedir_comp $2/plot2d/$f_base  $3/montage2d/$f_base
                done
            done
        fi
        
        echo "montage1d"
        if [ "$do_plot1d" = 1 -a "$_do_plot1d" = 1 -a "$do_plot2dclim" ] ; then
	        for var in $vars; do
	            echo -n $var" "
	    pngfile_clim=`basename $filename .nc`_clim_${var}.png
	    pngfile_climnorm=`basename $filename .nc`_climnorm_${var}.png
	    pngfile_ts=`basename $filename .nc`_ts_${var}.png
	    pngfile_tsnorm=`basename $filename .nc`_tsnorm_${var}.png
	    # TODO fix error message libgomp: Invalid value for environment variable OMP_NUM_THREADS
	    montage -geometry +1+20 -pointsize 20  -label $basedir_ref $1/plot2dclim/$pngfile_clim  -label diff $3/plot2dclim/$pngfile_clim  -label $basedir_ref $1/plotts/$pngfile_ts -label tsdiff $3/plotts/$pngfile_ts  -label $basedir_comp $2/plot2dclim/$pngfile_clim   -label diffnorm $3/plot2dclimnorm/$pngfile_climnorm  -label $basedir_comp $2/plotts/$pngfile_ts -label tsdiff_norm $3/plottsnorm/$pngfile_tsnorm $3/montage/`basename $filename .nc`_montage_${var}.png
	done
	fi
    fi
}

# func_diff - loop for all configurations and call func_diff
# func_diff $dir_ref $dir_comp $dir_diff $year_file $_do_plot_1d 1 $_do_plot_ref 1
function func_compare()
{
    basedir=$1
    index=$2
    echo "func_compare, basedir = "$basedir
    #cd ${basedir_prefix}/$basedir

    #loop domain
    for (( i_domains = 0 ; i_domains < ${#domains[*]} ; i_domains++ )) ; do      
        domain=${domains[$i_domains]}
        echo ++++++++++ domain: $domain

    #loop compilers
    for (( i_compilers = 0 ; i_compilers < ${#compile_compilers[*]} ; i_compilers++ )) ; do      
        compiler=${compile_compilers[$i_compilers]}
        
    # loop compiler flags
    for (( i_flags = 0 ; i_flags < $compile_numflags ; i_flags++ )) ; do      
        compile_flag=compile_flags_${compiler}[i_flags]; 
        echo == compiler: $compiler  flags: ${!compile_flag}       
        
        #loop nrun
        for (( i_nruns = 0 ; i_nruns < ${#nruns[*]} ; i_nruns++ )) ; do      
            nrun=${nruns[$i_nruns]}
            echo ++++++++++ nrun: $nrun

        # loop isimveg flag
        for (( i_isimvegs = 0 ; i_isimvegs < ${#isimvegs[*]} ; i_isimvegs++ )) ; do      
            isimveg=${isimvegs[$i_isimvegs]}
            if [ $isimveg -gt 2 ]; then isimveg0=1 ; else isimveg0=$isimveg; fi
            echo +++: isimveg = $isimveg

        # loop numthreads
        for (( i_numthreads = 0 ; i_numthreads < ${#numthreads[*]} ; i_numthreads++ )) ; do      

            numthread=${numthreads[$i_numthreads]}
            #export OMP_NUM_THREADS=$numthread          
            #if [ $numthreads_ref -ne 0 ] ; then numthread_ref=_omp_$numthreads_ref ; else numthread_ref="" ; fi
            numthread_ref=""
            if [ "$do_openmp" -eq 1 -a "$numthread" -ne 0 ] ; then numthread_comp=_omp_$numthread ; else numthread_comp="" ; fi
            echo numthread_ref: $numthread_ref
            echo numthread_comp: $numthread_comp

            #compare netcdf files
            dir_ref=${basedir_prefix}/compare/${basedir_ref}/$subdir/$domain/${nrun}yr/isimveg$isimveg0/${compiler}_${compile_flags_names[$i_flags]}${numthread_ref}
            dir_comp=${basedir_prefix}/compare/${basedir}/$subdir/$domain/${nrun}yr/isimveg$isimveg/${compiler}_${compile_flags_names[$i_flags]}${numthread_comp}
            dir_diff=${basedir_prefix}/diff/${basedir}_${basedir_ref}/$subdir/$domain/${nrun}yr/isimveg$isimveg/${compiler}_${compile_flags_names[$i_flags]}${numthread_comp}
            echo "dir_diff: "$dir_diff

	        if [ "$index" -gt 0 -o "$i_numthreads" -gt 0 ] ; then _do_plot_ref=0 ; else _do_plot_ref=1 ; fi

            year_file=inland-yearly.nc
            month_file=inland-monthly.nc
            
            # make merge files
            if [ "$do_cdo" = 1 -o "$command" = "all" ] ; then

		#if [ "$_do_plot_ref" -eq 0 ] ; then dirs="$dir_comp"; else dirs="$dir_ref $dir_comp"; fi
		if [ "$_do_plot_ref" -eq 0 ] ; then _dirs=( "$dir_comp" ); basedirs=( "basedir" ); else _dirs=( "$dir_ref" "$dir_comp" ); basedirs=( "$basedir_ref" "$basedir" ); fi
		echo "index: $index dirs: "${_dirs} #[*]}
		#for dir in $dirs; do
        for (( i_dirs = 0 ; i_dirs < ${#_dirs[*]} ; i_dirs++ )) ; do      
		    dir=${_dirs[$i_dirs]}
		    bdir=${basedirs[$i_dirs]}
		    echo "merging netcdf files in dir $dir"
		    rm -f ${dir}/$year_file 
		    rm -f ${dir}/$month_file 
                
            echo "bdir: $bdir"
            if [[ "${bdir:0:4}" != "sage"  ]] ; then 
#		        cdo -s -O -r -f nc mergetime ${dir}/ibis-yearly-????.nc ${dir}/$year_file >/dev/null 2>&1
#		        if [ ! -f ${dir}/$year_file ] ; then               
#			        cdo -s -O -r -f nc mergetime ${dir}/inland-yearly-????.nc ${dir}/$year_file >/dev/null 2>&1
#		        fi
                rm -f ${dir}/tmp1.nc
		        cdo -s -O -r -f nc mergetime ${dir}/ibis-yearly-????.nc ${dir}/tmp1.nc >/dev/null 2>&1
		        if [ ! -f ${dir}/tmp1.nc ] ; then               
			        cdo -s -O -r -f nc mergetime ${dir}/inland-yearly-????.nc ${dir}/tmp1.nc >/dev/null 2>&1
		        fi
		        if [ ! -f ${dir}/tmp1.nc ] ; then               
			        cdo -s -O -r -f nc mergetime ${dir}/output/inland-yearly-????.nc ${dir}/tmp1.nc >/dev/null 2>&1
		        fi
		        if [ ! -f ${dir}/tmp1.nc ] ; then               
			        cdo -s -O -r -f nc mergetime ${dir}/output/ibis-yearly-????.nc ${dir}/tmp1.nc >/dev/null 2>&1
		        fi
                
                mkdir -p ${dir}/tmp
                cdo -s -O -r -f nc splitname $dir/tmp1.nc ${dir}/tmp/""
                rm -f ${dir}/tmp/time_weights.nc ${dir}/tmp/caccount.nc ${dir}/tmp/totfall.nc
                cdo -s -O -r -f nc -b F64 merge ${dir}/tmp/*.nc ${dir}/$year_file # 2> /dev/null
                #rm -f ${dir}/tmp?.nc
                #rm -rf ${dir}/tmp
                
		        #cdo -s -O -r -f nc delname,time_weights $dir/tmp1.nc ${dir}/$year_file >/dev/null 2>&1
                #rm -f ${dir}/tmp1.nc
		    #cdo -s -O -r -f nc mergetime ${dir}/ibis-monthly-????.nc ${dir}/$month_file >/dev/null 2>&1
		    #if [ ! -f ${dir}/$month_file ] ; then               
			#cdo -s -O -r -f nc mergetime ${dir}/inland-monthly-????.nc ${dir}/$month_file >/dev/null 2>&1
		    #fi
            else
                #pwd
                #ls
                #rm -f ${dir}/tmp?.nc
                # TODO redirect stderr to null
                mkdir -p ${dir}/tmp
                #cdo -s -O -r -f nc merge ${dir}/output/yearly/*.nc ${dir}/tmp1.nc # 2> /dev/null
                #cdo -s -O -r -f nc -b F64 delname,zbot,ztop,rratio,tratio,time_weights_2,time_weights_3,time_weights_4,time_weights_5,time_weights_6,time_weights_7,time_weights_8,time_weights_9,time_weights_10,time_weights_11,time_weights_12,time_weights_13,time_weights_14,time_weights_15,time_weights_16 ${dir}/tmp1.nc ${dir}/$year_file               
                for f in ${dir}/output/yearly/*.nc ; do
                    cdo -s -O -r -f nc splitname $f ${dir}/tmp/""
                done
                rm -f ${dir}/tmp/zbot.nc ${dir}/tmp/ztop.nc ${dir}/tmp/rratio.nc ${dir}/tmp/tratio.nc ${dir}/tmp/time_weights.nc ${dir}/tmp/disturbf.nc ${dir}/tmp/sens.nc
                cdo -s -O -r -f nc -b F64 merge ${dir}/tmp/*.nc ${dir}/$year_file # 2> /dev/null
                #rm -f ${dir}/tmp?.nc
                rm -rf ${dir}/tmp
            fi
		done
            fi

            # compare / plot ref and comp
	    if [ "$nruns" -gt 1 ] ; then _do_plot_1d=1 ; else  _do_plot_1d=0 ; fi
            func_diff $dir_ref $dir_comp $dir_diff $year_file $_do_plot_1d 1 $_do_plot_ref 1
	    
	    # compare openp runs
           if [ "$do_plot_openmp" -eq 1 -a "$i_numthreads" -gt 0 ] ; then
               dir_ref=${basedir_prefix}/compare/${basedir}/$subdir/$domain/${nrun}yr/isimveg$isimveg/${compiler}_${compile_flags_names[$i_flags]}_omp_${numthreads[0]}
               dir_comp=${basedir_prefix}/compare/${basedir}/$subdir/$domain/${nrun}yr/isimveg$isimveg/${compiler}_${compile_flags_names[$i_flags]}_omp_${numthreads[$i_numthreads]}
               dir_diff=${basedir_prefix}/diff/${basedir}_openmp/$subdir/$domain/${nrun}yr/isimveg$isimveg/${compiler}_${compile_flags_names[$i_flags]}_omp_${numthreads[$i_numthreads]}_${numthreads[0]}
	           echo "compare openmp $i_numthreads"
               echo "basedir_ref: "$basedir_ref
               echo "basedir_comp: "$basedir_comp
               echo "dir_diff: "$dir_diff
               func_diff $dir_ref $dir_comp $dir_diff $year_file $_do_plot_1d 1 0 0
	   fi

        done #numthreads
        done #isimveg
        done #nruns
    done #compiler flags
    done #compilers
    done #domains

}

# func exec - run and diff (if requested)
function func_exec()
{

    if [ ! -f $plot1d ] ; then  echo "error, file $plot1d absent!" ; exit 1 ; fi
    if [ ! -f $plot2d ] ; then echo "error, file $plot2d absent!" ; exit 1 ; fi
    if [ ! -f $ifile_compar ] ; then echo "error, file $ifile_compar absent!" ; exit 1; fi
    #if [ ! -f $ifile_infile ] ; then echo "error, file $ifile_infile absent!" ; exit 1; fi
    
command=$1

###########################

# do the reference run

if [ "$command" = "run_ref" -o "$command" = "all" ]; then
if [ "$do_run_ref" = "1" ]; then
    echo "================================================================="
    echo "reference run"
    echo "================================================================="
    numthreads=( "${numthreads_ref}" )
    echo numthreads: $numthreads
#    if [ "$do_run_parallel" = "1" ]; then
#        func_run $basedir_ref &
#    else
        func_run $basedir_ref || exit $?
#    fi   
fi
fi

# do the comp runs

if [ "$command" = "run_comp" -o "$command" = "all" ]; then
if [ "$do_run_comp" = "1" ]; then
    echo "================================================================="
    echo "comp runs"
    echo "================================================================="
    
    if [ "$do_openmp" = 1 ]; then numthreads=( ${numthreads_comp[*]} ); else numthreads=("${numthreads_ref}"); fi

    cd $basedir_prefix
    
    for basedir_comp in $basedirs_comp; do
        echo "================================================================="
        echo "comp run: "$basedir_comp
        echo "================================================================="
        func_run $basedir_comp || exit $?
    done
fi
fi

if [ "$do_run_parallel" = "1" ]; then
    wait;
fi

cd $basedir_prefix

if [ "$command" = "diff" -o "$command" = "showdiff" -o "$command" = "all" ]; then
#if [ "$do_run_parallel" = "0" -o "$command" = "diff" ]; then
    echo "================================================================="
    echo "file differences"
    echo "================================================================="

    if [ "$command" = "showdiff" ]; then do_plot=0; fi

    if [ "$do_openmp" = 1 ]; then numthreads=( ${numthreads_comp[*]} ); else numthreads=("${numthreads_ref}"); fi

    index=0
    for basedir_comp in $basedirs_comp; do
        echo "================================================================="
        echo "comp run: "$basedir_comp
        echo "================================================================="
        func_compare $basedir_comp $index || exit $?
        let "index = $index + 1"
    done
#fi
fi

}

