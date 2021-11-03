# C3-Delphes-Studies
First time only

```
export workdir=/nfs_scratch/$USER/`date +%Y-%m-%d`/C3
mkdir -p $workdir
cd $workdir
cp /afs/hep.wisc.edu/home/dasu/public/public/c3-zhh-pythia8-delphes.txt .
wget https://launchpad.net/mg5amcnlo/3.0/3.2.x/+download/MG5_aMC_v3.2.0.tar.gz
tar zxf MG5_aMC_v3.2.0.tar.gz 
cd MG5_aMC_v3_2_0/
source /cvmfs/sft.cern.ch/lcg/views/LCG_97a/x86_64-centos7-gcc8-opt/setup.sh 
python bin/mg5_aMC
```

From within the MG5_aMC prompt execute the following. They take a long time 10 mins to finish.
Optionally, you can track the log files in the secondary login window, if you wish.

```
      install pythia8
      install Delphes
      exit
```

Content of c3-zhh-pythia8-delphes.txt

```
import model sm
generate e+ e- > z h h
output c3-zhh-pythia8-delphes
launch
shower=Pythia8
detector=Delphes
set mh 125.0
set wh 0.004
set ebeam1 275.
set ebeam2 275.
set nevents 1000
set iseed 1823211
```

To produce data later - create more .txt files with different configurations

```
export workdir=/nfs_scratch/$USER/`date +%Y-%m-%d`/C3/MG5_aMC_v3_2_0
cd $workdir
source /cvmfs/sft.cern.ch/lcg/views/LCG_97a/x86_64-centos7-gcc8-opt/setup.sh
export PYTHIA8DATA=$workdir/HEPTools/pythia8/share/Pythia8/xmldoc/
python bin/mg5_aMC c3-zhh-pythia8-delphes.txt
```

With above .txt file the event files will be in the directory $workdir/c3-zhh-pythia8-delphes/run_01/
However, by default Delphes does not use ILC detector configuation

For ILC detector configuration do the following to rerun with correct delphes configuration

```
cp $workdir/Delphes/cards/delphes_card_ILCgen.tcl $workdir/c3-zhh-pythia8-delphes/Cards/delphes_card.dat
cp -r $workdir/Delphes/cards/ILCgen $workdir/c3-zhh-pythia8-delphes/Cards/
cd $workdir/c3-zhh-pythia8-delphes/
python bin/madevent
```

From within madevent prompt do:

```
   delphes run_01
```

There should now be two root files in your directory
tag_1_delphes_events.root is generated using default configuration of Delphes
tag_2_delphes_events.root is generated using ILC configuration of Delphes
Use root to analyze the data, with the following as an example to follow:

https://twiki.cern.ch/twiki/bin/view/CMSPublic/MadgraphTutorial


