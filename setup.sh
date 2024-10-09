#/bin/bash
if [ -d '/cvmfs' ]; then
    if [[ `uname -r` == *"el9"* ]]; then
	source /cvmfs/sft.cern.ch/lcg/views/LCG_104b/x86_64-el9-gcc13-opt/setup.sh
	export mg5dir=/nfs_scratch/dasu/2024-10/MG5_aMC_v3_5_6/
    else
	echo "Currently Dasu is only supporting only el9 version on UW cluster"
	exit
    fi
fi
if [ -d "$ROOTSYS" ]; then
    echo "Using ROOT from $ROOTSYS";
else
    echo "ROOT is needed for this to work";
    exit
fi
if [ -d "$mg5dir" ]; then 
    echo "Using Madgraph5 from $mg5dir";
else
    echo "Madgraph is needed for this to work - export mg5dir=... and rerun";
    echo "Otherwise, you can just do analysis using root"
fi
export workdir=`dirname -- "${BASH_SOURCE[0]}"`
export basedir=`dirname $workdir`;
export confdir=$workdir/conf;
export datadir=$workdir/data;
export rootdir=$workdir/root;
mkdir -p $datadir
echo "Using workdir=$workdir";
echo "Using datadir=$datadir";
echo "Using confdir=$confdir";
echo "Using rootdir=$rootdir";
cd $workdir
echo "Current directory set to $PWD";
