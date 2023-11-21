#!/bin/bash
TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
DISPFREQ=$8
CONTINUE=$9
LOCK=${10}
VERSION=${11}
LPS=${12}
SETCUTOFF=${13}
ONLYPRE=${14}
SEEDS=${15}

# check if all variables defined (by checking the last one)
if test -z $ONLYPRE
then
    echo Skipping test since not all variables are defined
    echo "TSTNAME       = $TSTNAME"
    echo "BINNAME       = $BINNAME"
    echo "SETNAME       = $SETNAME"
    echo "BINID         = $BINID"
    echo "TIMELIMIT     = $TIMELIMIT"
    echo "NODELIMIT     = $NODELIMIT"
    echo "MEMLIMIT      = $MEMLIMIT"
    echo "DISPFREQ      = $DISPFREQ"
    echo "CONTINUE      = $CONTINUE"
    echo "LOCK          = $LOCK"
    echo "VERSION       = $VERSION"
    echo "LPS           = $LPS"
    echo "SETCUTOFF     = $SETCUTOFF"
    echo "ONLYPRE       = $ONLYPRE"
    echo "SEEDS         = $SEEDS"
    exit
fi


SETDIR=../settings

if test ! -e results
then
    mkdir results
fi
if test ! -e locks
then
    mkdir locks
fi

LOCKFILE=locks/$TSTNAME.$SETNAME.$VERSION.$LPS.lock
RUNFILE=locks/$TSTNAME.$SETNAME.$VERSION.$LPS.run.$BINID
DONEFILE=locks/$TSTNAME.$SETNAME.$VERSION.$LPS.done

if test "$LOCK" = "true"
then
    if test -e $DONEFILE
    then
	echo skipping test due to existing done file $DONEFILE
	exit
    fi
    if test -e $LOCKFILE
    then
	if test -e $RUNFILE
        then
	    echo continuing aborted run with run file $RUNFILE
	else
	    echo skipping test due to existing lock file $LOCKFILE
	    exit 1
	fi
    fi
    date > $LOCKFILE
    date > $RUNFILE
fi


for ((s = 0; ${s} <= ${SEEDS}; s++))
do
    FILEBASENAME=check.$TSTNAME.$BINID.$SETNAME
    if [ ${SEEDS} -ne 0 ]
    then
        FILEBASENAME=$FILEBASENAME.s${s}
    fi
    OUTFILE=results/$FILEBASENAME.out
    ERRFILE=results/$FILEBASENAME.err
    RESFILE=results/$FILEBASENAME.res
    TEXFILE=results/$FILEBASENAME.tex
    SETFILE=results/$FILEBASENAME.set

    if test ! -e $OUTFILE
    then
        CONTINUE=false
    fi

    if test "$CONTINUE" = "true"
    then
        MVORCP=cp
    else
        MVORCP=mv
    fi

    DATEINT=`date +"%s"`
    if test -e $OUTFILE
    then
        $MVORCP $OUTFILE $OUTFILE.old-$DATEINT
    fi
    if test -e $ERRFILE
    then
        $MVORCP $ERRFILE $ERRFILE.old-$DATEINT
    fi

    # if cutoff should be passed, check for solu file
    SOLUFILE=""
    if test $SETCUTOFF == "true"
    then
        SOLUFILE="testset/"$TSTNAME.solu
        if test ! -e $SOLUFILE
        then
            echo "Solu file ("$SOLUFILE") not available"
            exit
        fi
    fi

    # check for only preprocessing option
    PRESTR=""
    if test $ONLYPRE == "true"
    then
        PRESTR="-O"
    fi

    # check for seeds option
    SEEDSSTR="-seed "${s}

    if test "$CONTINUE" = "true"
    then
        LASTPROB=`./getlastprob.awk $OUTFILE`
        echo Continuing benchmark. Last solved instance: $LASTPROB
        echo "" >> $OUTFILE
        echo "----- Continuing from here. Last solved: $LASTPROB -----" >> $OUTFILE
        echo "" >> $OUTFILE
    else
        LASTPROB=""
    fi

    uname -a >>$OUTFILE
    uname -a >>$ERRFILE
    date >>$OUTFILE
    date >>$ERRFILE

    # set settings name if necessary
    if test "$SETNAME" = "default"
    then
        SETTINGS=""
    else
        if test -f $SETDIR/$SETNAME".set"
        then
        SETTINGS="-s "$SETDIR"/"$SETNAME".set"
        else
        echo skipping test due to non-existing settings file $SETDIR"/"$SETNAME".set"
        exit 1
        fi
    fi

    # we add 100% to the hard time limit and additional 600 seconds in case of small time limits
    HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 600\` + $TIMELIMIT`

    # values of MEMLIMIT are in MB
    # values of HARDMEMLIMIT are in kilobyte
    HARDMEMLIMIT=`expr $MEMLIMIT \* 1024`


    # process test instances
    STOREFILENAME=""
    for i in `cat testset/$TSTNAME.test` DONE
    do
        if test "$i" = "DONE"
        then
        date > $DONEFILE
        break
        fi

        FILENAME=$i

        if test "$LASTPROB" = ""
        then
        LASTPROB=""
        if test -f $FILENAME
        then
            SHORTFILENAME=`basename $i .gz`
            SHORTFILENAME=`basename $SHORTFILENAME .mps`
            SHORTFILENAME=`basename $SHORTFILENAME .lp`
            SHORTFILENAME=`basename $SHORTFILENAME .opb`

            echo @01 $FILENAME ===========
            echo @01 $FILENAME ===========                >> $ERRFILE

            # set objective limit: optimal solution value from solu file, if existent
            CUTOFFSTR=""
            if test $SETCUTOFF == "true"
            then
            CUTOFF=`grep "$SHORTFILENAME " $SOLUFILE | grep -v =feas= | grep -v =inf= | tail -n 1 | awk '{print $3}'`
            CUTOFFSTR="-setcutoff "$CUTOFF
            fi

            echo -----------------------------
            date
            date >>$ERRFILE
            echo -----------------------------
            date +"@03 %s"
            echo " ulimit -f 200000; ulimit -t $HARDTIMELIMIT; ulimit -v $HARDMEMLIMIT; ../$BINNAME $FILENAME $SETTINGS -t $TIMELIMIT -m $MEMLIMIT -n $NODELIMIT -d $DISPFREQ $CUTOFFSTR $SEEDSSTR $PRESTR" 2>>$ERRFILE
            bash -c " ulimit -f 200000; ulimit -t $HARDTIMELIMIT; ulimit -v $HARDMEMLIMIT; ../$BINNAME $FILENAME $SETTINGS -t $TIMELIMIT -m $MEMLIMIT -n $NODELIMIT -d $DISPFREQ $CUTOFFSTR $SEEDSSTR $PRESTR" 2>>$ERRFILE
            date +"@04 %s"
            echo -----------------------------
            date
            date >>$ERRFILE
            echo -----------------------------
            echo
            echo =ready=
        else
            echo @02 FILE NOT FOUND: $FILENAME ===========
            echo @02 FILE NOT FOUND: $FILENAME =========== >>$ERRFILE
        fi
        else
        echo skipping $FILENAME
        if test "$LASTPROB" = "$FILENAME"
        then
            LASTPROB=""
            fi
        fi
    done | tee -a $OUTFILE

    date >>$OUTFILE
    date >>$ERRFILE
done #end for SEEDS



if test -e $DONEFILE
then
    ./evalcheck.sh $OUTFILE

    if test "$LOCK" = "true"
    then
	rm -f $RUNFILE
    fi
fi
