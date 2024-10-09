#/bin/bash

# If the node is busy give up running

# Setup Madgraph-Pythia-Delphes
source /nfs_scratch/dasu/2023-02-14/Madgraph-Pythia-Delphes-Studies/setup.sh

# Loop over all input files received
export inpdir=$HOME/public/Madgraph-Pythia-Delphes-Data
datestr=`date +%Y%m%d-%H%M%S`
export datadir=$inpdir/$datestr
cd $inpdir
for file in *.inp; do
    if [ -f "$file" ]
    then
	echo "Processing $file..."
	# Create unique output data directory for this run
	mkdir -p $datadir
	mv $file $datadir
	cd $datadir
	stdout=`basename $file .inp`.out
	stderr=`basename $file .inp`.err
	nohup time python $mg5dir/bin/mg5_aMC $file < /dev/null > $stdout 2> $stderr&
	cd $inpdir
    fi
done
