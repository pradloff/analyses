DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

unset TMP
unset TMPDIR

export ANALYSISHOME=$DIR
export ROOTCOREDIR=$DIR/external
export PYTHONPATH=$ANALYSISHOME:$PYTHONPATH
