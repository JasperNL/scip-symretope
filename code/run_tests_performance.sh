#!/bin/bash
if [ -z "$3" ]
then
    echo "./run_tests_performance.sh {start, >= 0} {end, <= number of instances} {number of seeds, optional, default: 1}"
    exit 1
fi

IFS=$'\r\n' GLOBIGNORE='*' command eval  'INSTANCES=($(cat testsfinal.inst))'
echo "${#INSTANCES[@]} settings."
if [ $1 -lt "0" ]
then
    echo "First argument $1 is smaller than 0, out of range."
    exit 1
fi
if [ $2 -ge ${#INSTANCES[@]} ]
then
    echo "Second argument $2 is larger or equal to ${#INSTANCES[@]}}, out of range."
    exit 1
fi

NSEEDS=${3:-1}

for i in $(seq $1 $2)
do
    elem=${INSTANCES[$i]}
    read -a arr <<< "$elem"  # uses default whitespace IFS
    TESTSET=$(basename ${arr[0]} .test)
    SETTINGS=$(basename ${arr[1]} .set)
    USESOL=${arr[2]}
    printf "# Running instance %3d \n TESTSET  : %s\n SETTINGS : %-50s\n USESOL   : %-5s\n" \
        "$i" "${TESTSET}" "${SETTINGS}" "${USESOL}"

    # SCIPs default testcluster make command variant
    # make \
    #     OPT=opt \
    #     LPS=spx \
    #     SYM=bliss \
    #     TEST=$TESTSET \
    #     SETTINGS=$SETTINGS \
    #     TIME=7200 \
    #     SETUSESOL=$USESOL \
    #     testcluster

    # TU/e cluster
    cd check
    for ((SEED=0; SEED < $NSEEDS; SEED++))
    do
        ./run_all_tue.sh $TESTSET $SETTINGS 7200 $SEED ${TESTSET}__${SETTINGS}__${SEED}__
    done
    cd ..
done
