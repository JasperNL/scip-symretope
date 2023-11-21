#!/bin/bash
export TSTNAME=$1;
export SETTINGS=$2;
export OUTPUTNAME=$3;
export BINNAME="bin/sbcs.linux.x86_64.gnu.opt.spx2"

function runfile()
{
    FILENAME=$1
    BASENAME=$(basename $FILENAME)
    echo "RUNNING $FILENAME; $BASENAME; $SETTINGS; $OUTPUTNAME;"
    echo @02 $BASENAME > results/$OUTPUTNAME$BASENAME.out;
    ../$BINNAME $FILENAME -s "../settings/"$SETTINGS".set" -O >> results/$OUTPUTNAME$BASENAME.out 2>&1
}
export -f runfile

xargs -P 12 -d $'\n' -n 1 bash -c '
  runfile $1
  ' _ < testset/$TSTNAME.test
