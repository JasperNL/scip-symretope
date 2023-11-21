#!/bin/bash
export TSTNAME=$1;
export SETTINGS=$2;
export TIMELIMIT=$3;
export SEED=$4;
export OUTPUTNAME=$5;

echo "TSTNAME: ${TSTNAME}"
echo "SETTINGS: ${SETTINGS}"
echo "TIMELIMIT: ${TIMELIMIT}"
echo "SEED: ${SEED}"
echo "OUTPUTNAME: ${OUTPUTNAME}"

while read filename; do
    INSTANCEPATH=`dirname $filename`
    SHORTFILENAME=`basename $filename .gz`
    SHORTFILENAME=`basename $SHORTFILENAME .mps`
    SHORTFILENAME=`basename $SHORTFILENAME .lp`
    SHORTFILENAME=`basename $SHORTFILENAME .opb`
    SHORTFILENAME=`basename $SHORTFILENAME .cip`
    SOLFILENAME="$INSTANCEPATH/$SHORTFILENAME.sol"

    # Check if sol-file exists, and if so, do not need to run.
    if [ -f "$SOLFILENAME" ]
    then
        echo "$SOLFILENAME exists already, skipping."
    else
        echo "sbatch run_tue.sh $filename $SETTINGS $TIMELIMIT $SEED $OUTPUTNAME"
        sbatch run_tue.sh $filename $SETTINGS $TIMELIMIT $SEED $OUTPUTNAME
        # ./run_tue.sh $filename $SETTINGS $TIMELIMIT $SEED $OUTPUTNAME
    fi

done < testset/$TSTNAME.test
