#!/bin/sh
#$ -V
#$ -cwd
#$ -l mem_free=2G
#$ -t 1-70

# Short jobs, your cluster may require additional resource
# specifications (current header for SGE). The number of
# in array == number of bigWig files being used (70 for H1).

# Set bigWig directory
array=(`ls ./bigWigs/*`)

# e.g. for H1 cell type, run for each type of split boundary
for inputfile in $( ls h*b.bins ); do

    # each job picks one bigWig from directory
    filename="${array[$SGE_TASK_ID-1]}"

    exprname=$( basename $filename | sed 's/wgEncode//g' )

    export JOBOUTPUT_DIR=htads
    export JOBOUTPUT=$JOBOUTPUT_DIR/${inputfile}_$exprname

    [ -d $JOBOUTPUT_DIR ] || mkdir -p $JOBOUTPUT_DIR

    # run bigWigAverageOverBed compiled binary for your platform
    ./bigWigAv_linux_86 $filename $inputfile $JOBOUTPUT
done
