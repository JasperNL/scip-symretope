#!/usr/bin/env bash
#
# Call with "make testcluster"
#
# The queue is passed via $QUEUE (possibly defined in a local makefile in scip/make/local).
#
# For each run, we can specify the number of nodes reserved for a run via $PPN. If tests runs
# with valid time measurements should be executed, this number should be chosen in such a way
# that a job is run on a single computer, i.e., in general, $PPN should equal the number of cores
# of each computer. Of course, the value depends on the specific computer/queue.
#
# To get the result files call "./evalcheck_cluster.sh
# results/check.$TSTNAME.$BINNAME.$SETNAME.eval in directory check/
# This leads to result files
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.out
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.res
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.err

TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
DISPFREQ=$8
CONTINUE=$9
VERSION=${10}
LPS=${11}
QUEUE=${12}
QUEUETYPE=${13}
PPN=${14}
CLIENTTMPDIR=${15}
NOWAITCLUSTER=${16}
PERMUTE=${17}
SETCUTOFF=${18}
SEEDS=${19}
SETUSESOL=${20}

# check if all variables defined (by checking the last one)
if test -z $SETUSESOL
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
    echo "VERSION       = $VERSION"
    echo "LPS           = $LPS"
    echo "QUEUE         = $QUEUE"
    echo "PPN           = $PPN"
    echo "CLIENTTMPDIR  = $CLIENTTMPDIR"
    echo "NOWAITCLUSTER = $NOWAITCLUSTER"
    echo "PERMUTE       = $PERMUTE"
    echo "SETCUTOFF     = $SETCUTOFF"
    echo "SEEDS         = $SEEDS"
    echo "SETUSESOL     = $SETUSESOL"
    exit 1;
fi


# get current SCIP path
SCIPPATH=`pwd`

if test ! -e $SCIPPATH/results
then
    mkdir $SCIPPATH/results
fi

# check if the settings file exists
if test $SETNAME != "default"
then
    if test ! -e $SCIPPATH/../settings/$SETNAME.set
    then
        echo Skipping test since the settings file $SCIPPATH/../settings/$SETNAME.set does not exist.
        exit
    fi
    SETTINGS=$SCIPPATH/../settings/$SETNAME.set
else
    SETTINGS=""
fi

# check if binary exists
if test ! -e $SCIPPATH/../$BINNAME
then
    echo Skipping test since the binary $BINNAME does not exist.
    exit
fi

# check if queue has been defined
if test "$QUEUE" = ""
then
    echo Skipping test since the queue name has not been defined.
    exit
fi

# check if number of nodes has been defined
if test "$PPN" = ""
then
    echo Skipping test since the number of nodes has not been defined.
    exit
fi

# if cutoff should be passed, check for solu file
if test $SETCUTOFF == "true"
then
    SOLUFILE="testset/"$TSTNAME.solu
    if test ! -e $SOLUFILE
    then
        echo "Solu file ("$SOLUFILE") not available"
        exit
    fi
fi

# we add 100% to the hard time limit and additional 600 seconds in case of small time limits
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 600\` + $TIMELIMIT`

# we add 10% to the hard memory limit and additional 100mb to the hard memory limit
# the srun queue requires a format duration HH:MM:SS (and optionally days) and megabytes,
# whereas the qsub requires the memory limit in kB
if test "$QUEUETYPE" != "qsub"
then
    # values of MEMLIMIT are in MB
    HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``

    #format is (d-)HH:MM:SS
    TMP=`expr $HARDTIMELIMIT`
    HARDTIMELIMIT=""
    DIVISORS=(60 60 24)
    for((i=0; i<=2; i++))
    do
        printf -v HARDTIMELIMIT "%02d${HARDTIMELIMIT}" `expr ${TMP} % ${DIVISORS[i]}`
        # separate the numbers by colons except for the last (HH hours)
        if test $i -lt 2
        then
            HARDTIMELIMIT=":${HARDTIMELIMIT}"
        fi
        TMP=`expr ${TMP} / ${DIVISORS[i]}`
    done
    if test ${TMP} -gt 0
    then
        HARDTIMELIMIT=${TMP}-${HARDTIMELIMIT}
    fi
else
    # values of MEMLIMIT are in MB
    # values of HARDMEMLIMIT are in byte for qsub
    HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``
    HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`
fi

# counter to define file names for a test set uniquely
COUNT=1

# loop over permutations
for ((s=0; $s <= $SEEDS; s++))
do
    for ((p = 0; $p <= $PERMUTE; p++))
    do
        FILEPOSTFIX=""
        # if number of permutations is positive, add postfix
        if test $PERMUTE -gt 0
        then
        FILEPOSTFIX="$FILEPOSTFIX.p$p"
        fi

        # if number of seeds is positive, add postfix
        if test $SEEDS -gt 0
        then
        FILEPOSTFIX="$FILEPOSTFIX.s$s"
        fi

        # if USESOL, add postfix.
        if test $SETUSESOL == "true"
        then
        FILEPOSTFIX="$FILEPOSTFIX.usesol"
        fi

        EVALFILE=$SCIPPATH/results/check.$QUEUE.$TSTNAME.$BINID.$SETNAME$FILEPOSTFIX.eval
        echo > $EVALFILE

        for i in `cat testset/$TSTNAME.test` DONE
        do
            if test "$i" = "DONE"
            then
                break
            fi

                # check if problem instance exists
            if test -f $SCIPPATH/$i
            then

                    # the cluster queue has an upper bound of 2000 jobs; if this limit is
                    # reached the submitted jobs are dumped; to avoid that we check the total
                    # load of the cluster and wait until it is save (total load not more than
                    # 1900 jobs) to submit the next job.
                    if test "$NOWAITCLUSTER" != "1"
                    then
                ./waitcluster.sh 1600 $QUEUE 200
                fi

                INSTANCEPATH=`dirname $i`
                SHORTFILENAME=`basename $i .gz`
                SHORTFILENAME=`basename $SHORTFILENAME .mps`
                SHORTFILENAME=`basename $SHORTFILENAME .lp`
                SHORTFILENAME=`basename $SHORTFILENAME .opb`
                SHORTFILENAME=`basename $SHORTFILENAME .cip`
                JOBNAME="sym-"$SHORTFILENAME

                # if number of permutations is positive, add postfix
                FILENAME=$USER.$QUEUE.$TSTNAME.$COUNT"_"$SHORTFILENAME.$BINID.$SETNAME$FILEPOSTFIX
                BASENAME=$SCIPPATH/results/$FILENAME

                TMPFILE=$BASENAME.tmp
                SETFILE=$BASENAME.set

                echo $BASENAME >> $EVALFILE

                COUNT=`expr $COUNT + 1`

                # in case we want to continue we check if the job was already performed
                if test "$CONTINUE" != "false"
                then
                    if test -e results/$FILENAME.out
                    then
                        echo skipping file $i due to existing output file $FILENAME.out
                        continue
                    fi
                fi

                # set objective limit: optimal solution value from solu file, if existent
                CUTOFF=""
                if test $SETCUTOFF == "true"
                then
                    if test $SOLUFILE == ""
                    then
                        echo Exiting test because no solu file can be found for this test
                        exit
                    fi
                    CUTOFF=`grep "$SHORTFILENAME " $SOLUFILE | grep -v =feas= | grep -v =inf= | tail -n 1 | awk '{print $3}'`
                fi

                USESOL=""
                if test $SETUSESOL == "true"
                then
                    USESOL="$INSTANCEPATH/$SHORTFILENAME.sol"
                    # echo $USESOL;
                    if [ ! -f $USESOL ]
                    then
                        echo "ERROR: $SHORTFILENAME: $USESOL does not exist.";
                        exit
                    fi
                fi

                # check queue type
                if test  "$QUEUETYPE" != "qsub"
                then
                    # additional environment variables needed by runcluster.sh
                    export SOLVERPATH=$SCIPPATH
                    export BINNAME
                    export BASENAME=$FILENAME
                    export FILENAME=$SCIPPATH/$i
                    export CLIENTTMPDIR
                    export HARDTIMELIMIT
                    export HARDMEMLIMIT
                    export CHECKERPATH=$SCIPPATH/solchecker
                    export SETTINGS
                    export TIMELIMIT
                    export NODELIMIT
                    export MEMLIMIT
                    export CUTOFF
                    export DISPFREQ
                    export PERM=$p
                    export SEED=$s
                    export USESOL
                    # bash run.sh
                    # echo sbatch --job-name=$JOBNAME --mem=$HARDMEMLIMIT -p $QUEUE -A dopt --exclusive --time=$HARDTIMELIMIT --output=/dev/null run.sh
                    sbatch --job-name=$JOBNAME --mem=$HARDMEMLIMIT -p $QUEUE -A dopt --exclusive --time=$HARDTIMELIMIT --output=/dev/null run.sh
                else
                    # -v to set local variables explicitly -V to copy all environment variables
                    qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -l nodes=1:ppn=$PPN -N sym_$SHORTFILENAME -v SOLVERPATH=$SCIPPATH,BINNAME=$BINNAME,FILENAME=$SCIPPATH/$i,,BASENAME=$FILENAME,CLIENTTMPDIR=$CLIENTTMPDIR,SETTINGS=$SETTINGS,TIMELIMIT=$TIMELIMIT,NODELIMIT=$NODELIMIT,MEMLIMIT=$MEMLIMIT,CUTOFF=$CUTOFF,HARDMEMLIMIT=$HARDMEMLIMIT,DISPFREQ=$DISPFREQ,PERM=$p -V -q $QUEUE -o /dev/null -e /dev/null runcluster.sh
                fi
            else
                echo "input file "$SCIPPATH/$i" not found!"
            fi
        done
    done
done
