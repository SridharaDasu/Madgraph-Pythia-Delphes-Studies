# C3-Delphes-Studies
First time only

```
if [ -d '/nfs_scratch' ]; then export basedir=/nfs_scratch/$USER/`date +%Y-%m-%d`; else basedir=$PWD/`date +%Y-%m-%d`; fi
mkdir -p $basedir
cd $basedir
git clone https://github.com/SridharaDasu/C3-Delphes-Studies.git
export workdir=$basedir/C3-Delphes-Studies
cd $workdir
mkdir -p $workdir/data
export datadir=$workdir/data
wget https://launchpad.net/mg5amcnlo/3.0/3.2.x/+download/MG5_aMC_v3.2.0.tar.gz
tar zxf MG5_aMC_v3.2.0.tar.gz 
export mg5dir=$workdir/MG5_aMC_v3_2_0/
cd $mg5dir
if [ -d '/cvmfs' ]; then source /cvmfs/sft.cern.ch/lcg/views/LCG_97a/x86_64-centos7-gcc8-opt/setup.sh; fi
if [ `uname` == 'Darwin' ]; then echo export MACOSX_DEPLOYMENT_TARGET=10.15; fi
python $mg5dir/bin/mg5_aMC
```

From within the MG5_aMC prompt execute the following. They take a long time 10 mins to finish.
Optionally, you can track the log files in the secondary login window, if you wish.

```
      install pythia8
      install Delphes
      exit
```

On relogin cd to the base directory, i.e., the directory with the date of creation above, e.g., /Users/dasu/2021-11-03/

```
export basedir=$PWD
export workdir=$basedir/C3-Delphes-Studies
export datadir=$workdir/data
export mg5dir=$workdir/MG5_aMC_v3_2_0/
```

To produce data use  .txt files with different configurations

```
cd $datadir
if [ -d '/cvmfs' ]; then source /cvmfs/sft.cern.ch/lcg/views/LCG_97a/x86_64-centos7-gcc8-opt/setup.sh; fi
export PYTHIA8DATA=$mg5dir/HEPTools/pythia8/share/Pythia8/xmldoc/
python $mg5dir/bin/mg5_aMC $workdir/c3-zhh-pythia8-delphes.txt
```

With above .txt file the event files will be in the directory $mg5dir/c3-zhh-pythia8-delphes/run_01/
However, by default Delphes does not use ILC detector configuation

For ILC detector configuration do the following to rerun with correct delphes configuration

```
cp $mg5dir/Delphes/cards/delphes_card_ILCgen.tcl $datadir/c3-zhh-pythia8-delphes/Cards/delphes_card.dat
cp -r $mg5dir/Delphes/cards/ILCgen $datadir/c3-zhh-pythia8-delphes/Cards/
cd $datadir/c3-zhh-pythia8-delphes/
python bin/madevent
```

From within madevent prompt do:

```
   delphes run_01
```

There should now be two root files in your directory $datadir/Events/run_01/
tag_1_delphes_events.root is generated using default configuration of Delphes
tag_2_delphes_events.root is generated using ILC configuration of Delphes
Use root to analyze the data, with the following as an example to follow:

https://twiki.cern.ch/twiki/bin/view/CMSPublic/MadgraphTutorial


