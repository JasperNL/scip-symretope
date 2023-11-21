#!/usr/bin/env bash

# check if tmp-path exists
if test ! -d $CLIENTTMPDIR/${USER}-tmpdir
then
    mkdir $CLIENTTMPDIR/${USER}-tmpdir
    echo Creating directory $CLIENTTMPDIR/${USER}-tmpdir for temporary outfile
fi

OUTFILE=$CLIENTTMPDIR/${USER}-tmpdir/$BASENAME.out
ERRFILE=$CLIENTTMPDIR/${USER}-tmpdir/$BASENAME.err
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

PERMSTR=""
if test $PERM -gt 0
then
    PERMSTR="-p "$PERM
fi

SEEDSTR=""
if test $SEED -gt 0
then
    SEEDSTR="-seed "$SEED
fi

SETTINGSSTR=""
if test "$SETTINGS" != ""
then
    SETTINGSSTR="-s "$SETTINGS
fi

CUTOFFSTR=""
if test "$CUTOFF" != ""
then
    CUTOFFSTR="-setcutoff "$CUTOFF
fi

USESOLSTR=""
if test "$USESOL" != ""
then
    USESULSTR="-l "$USESOL
fi

echo $SOLVERPATH/../$BINNAME $FILENAME -t $TIMELIMIT -m $MEMLIMIT -n $NODELIMIT -d $DISPFREQ $PERMSTR $SEEDSTR $SETTINGSSTR $CUTOFFSTR $USESULSTR >>$OUTFILE 2>>$ERRFILE
$SOLVERPATH/../$BINNAME $FILENAME -t $TIMELIMIT -m $MEMLIMIT -n $NODELIMIT -d $DISPFREQ $PERMSTR $SEEDSTR $SETTINGSSTR $CUTOFFSTR $USESULSTR >>$OUTFILE 2>>$ERRFILE

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
