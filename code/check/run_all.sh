#!/bin/bash
export TSTNAME=$1;
export SETTINGS=$2;
export TIMELIMIT="${3:-240}";
export OUTPUTNAME=$4;
export BINNAME="bin/sbcs.linux.x86_64.gnu.opt.spx2"
export NUMTHREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}" # take SLURM_CPUS_PER_TASK, otherwise run "nproc".

echo "TSTNAME: ${TSTNAME}"
echo "SETTINGS: ${SETTINGS}"
echo "TIMELIMIT: ${TIMELIMIT}"
echo "OUTPUTNAME: ${OUTPUTNAME}"
echo "BINNAME: ${BINNAME}"
echo "NUMTHREADS: ${NUMTHREADS}"

echo ""
echo "Starting..."

function runfile()
{
    FILENAME=$1
    BASENAME=$(basename $FILENAME)
    echo "RUNNING $FILENAME; $BASENAME; $SETTINGS; $OUTPUTNAME;"
    echo @01 $FILENAME ===========      >> results/$OUTPUTNAME$BASENAME.out
    echo @02 $BASENAME                  >> results/$OUTPUTNAME$BASENAME.out
    date +"@03 %s"                      >> results/$OUTPUTNAME$BASENAME.out
    ../$BINNAME $FILENAME -s "../settings/"$SETTINGS".set" -t $TIMELIMIT -w results/$OUTPUTNAME$BASENAME.sol >> results/$OUTPUTNAME$BASENAME.out 2>&1
    date +"@04 %s"                      >> results/$OUTPUTNAME$BASENAME.out
}
export -f runfile

xargs -P $NUMTHREADS -d $'\n' -n 1 bash -c '
  runfile $1
  ' _ < testset/$TSTNAME.test

echo "Done..."
