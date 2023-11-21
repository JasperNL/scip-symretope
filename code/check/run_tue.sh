#!/bin/bash
#SBATCH --partition=mcs.default.q
#SBATCH --cpus-per-task=1

source $HOME/prepare_scip_dev.sh

FILENAME=$1
SETTINGS=$2
TIMELIMIT=$3
SEED=$4
OUTPUTNAME=$5

BINNAME="bin/sbcs.linux.x86_64.gnu.opt.cpx"
BASENAME=$(basename $FILENAME)


echo "FILENAME:   ${FILENAME}"
echo "SETTINGS:   ${SETTINGS}"
echo "TIMELIMIT:  ${TIMELIMIT}"
echo "OUTPUTNAME: ${OUTPUTNAME}"
echo "BINNAME:    ${BINNAME}"
echo "BASENAME:   ${BASENAME}"
echo "SEED:       ${SEED}"

echo "FILENAME:   ${FILENAME}"      >> results/$OUTPUTNAME$BASENAME.out
echo "SETTINGS:   ${SETTINGS}"      >> results/$OUTPUTNAME$BASENAME.out
echo "TIMELIMIT:  ${TIMELIMIT}"     >> results/$OUTPUTNAME$BASENAME.out
echo "OUTPUTNAME: ${OUTPUTNAME}"    >> results/$OUTPUTNAME$BASENAME.out
echo "BINNAME:    ${BINNAME}"       >> results/$OUTPUTNAME$BASENAME.out
echo "BASENAME:   ${BASENAME}"      >> results/$OUTPUTNAME$BASENAME.out
echo "SEED:       ${SEED}"          >> results/$OUTPUTNAME$BASENAME.out

echo @01 $FILENAME ===========      >> results/$OUTPUTNAME$BASENAME.out
echo @02 $BASENAME                  >> results/$OUTPUTNAME$BASENAME.out
date +"@03 %s"                      >> results/$OUTPUTNAME$BASENAME.out
../$BINNAME $FILENAME -s "../settings/"$SETTINGS".set" -t $TIMELIMIT -w results/$OUTPUTNAME$BASENAME.sol -seed $SEED >> results/$OUTPUTNAME$BASENAME.out 2>&1
date +"@04 %s"                      >> results/$OUTPUTNAME$BASENAME.out
