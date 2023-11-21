#!/usr/bin/env bash

# check if tmp-path exists
if test ! -d $CLIENTTMPDIR
then
    echo Skipping test since the path for the tmp-dir does not exist.
    exit
fi

OUTFILE=$CLIENTTMPDIR/$BASENAME.out
ERRFILE=$CLIENTTMPDIR/$BASENAME.err
TMPFILE=$SOLVERPATH/results/$BASENAME.tmp

uname -a                            > $OUTFILE
uname -a                            > $ERRFILE
echo @01 $FILENAME ===========      >> $OUTFILE
echo @01 $FILENAME ===========      >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
date                                >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date +"@03 %s"                      >> $OUTFILE

# set memory limits
# compute memory limit in KB (HARDMEMLIMIT is in MB):
HARDMEMLIMIT=`expr $HARDMEMLIMIT / 1000`
ulimit -m $HARDMEMLIMIT
ulimit -v $HARDMEMLIMIT

PERMSTR = ""
if test $PERM -gt 0
then
    PERMSTR="-p "$PERM
fi

SETTINGSSTR=""
if test "$SETTINGS" != ""
then
    SETTINGSSTR="-s "$SETTINGS
fi

$CUTOFFSTR = ""
if test "$CUTOFF" != ""
then
    CUTOFFSTR="-setcutoff "$CUTOFF
fi

echo $SOLVERPATH/../$BINNAME ${FILENAME}X.txt ${FILENAME}Y.txt -t $TIMELIMIT -m $MEMLIMIT -n $NODELIMIT -d $DISPFREQ $PERMSTR $SETTINGSSTR $CUTOFFSTR >>$OUTFILE 2>>$ERRFILE
$SOLVERPATH/../$BINNAME ${FILENAME}X.txt ${FILENAME}Y.txt -t $TIMELIMIT -m $MEMLIMIT -n $NODELIMIT -d $DISPFREQ $PERMSTR $SETTINGSSTR $CUTOFFSTR >>$OUTFILE 2>>$ERRFILE

date +"@04 %s"                      >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $ERRFILE
echo                                >> $OUTFILE
echo =ready=                        >> $OUTFILE

mv $OUTFILE $SOLVERPATH/results/$BASENAME.out
mv $ERRFILE $SOLVERPATH/results/$BASENAME.err

rm -f $TMPFILE
#chmod g+r $ERRFILE
#chmod g+r $SCIPPATH/results/$BASENAME.out
#chmod g+r $SCIPPATH/results/$BASENAME.set
