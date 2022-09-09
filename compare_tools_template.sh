#! /bin/bash

# Copy this file to your directory and fill variables with your current file names and locations

# If multiple - insert it in a space-separated manner
CHECKM_RESULTS=''
TOOLS_LABELS=''
MIN_COMPLETENESS='0.95'
MIN_PURITY='0.95'
BINNINGS=''
DEPTHS_FILE=''
MINDEPTH='100'

while getopts ":o:" opt; do
    case $opt in
        o)
            OUTPUT=$OPTARG;;
        \?)
            echo "Invalid parameter!"
            exit;;
    esac
done

python3 compare_checkm_results.py \
-i $CHECKM_RESULTS \
-l $TOOLS_LABELS \
-c $MIN_COMPLETENESS \
-p $MIN_PURITY \
-b $BINNINGS \
-d $DEPTHS_FILE \
--mindepth $MINDEPTH \
-o $OUTPUT
