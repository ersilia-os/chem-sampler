INPUT_FILE=$1
OUTPUT_FILE=$2
TRACE_FILE=$3
CHECKPOINT_FILE=$4

if [ $# -eq 3 ] ; then
    echo "Filename not supplied."
else

python molecule_generation/cli/train.py MoLeR $TRACE_FILE \--load-saved-model $CHECKPOINT_FILE \ --load-weights-only

fi
