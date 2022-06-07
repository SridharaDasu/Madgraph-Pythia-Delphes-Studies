# Madgraph-Pythia-Delphes-Studies

This package contains instructions for
1) producing Madgraph5 based signals,
2) hadronize them using Pythia8,
3) simulate for ATLAS/CMS/C3/MuonCol (or other) detector's response using Delphes,
4) submit jobs to Condor in UW-Madison HEP group computing environment to produce large number of events
5) Sample analyzer using the resulting Delphes root file

The inputs are madgraph, pythia8 and Delphes control files for the signal root file producer are named *.txt.

## If you are working on login.hep.wisc.edu, you can use my existing madgraph:

### Code installation (On first login only):

```
export mg5dir=/nfs_scratch/dasu/CentOS7/MG5_aMC_v3_2_0/
export basedir=/nfs_scratch/$USER/`date +%Y-%m-%d`
mkdir -p $basedir
git clone git@github.com:SridharaDasu/Madgraph-Pythia-Delphes-Studies.git
source $basedir/Madgraph-Pythia-Delphes-Studies/setup.sh
```

### Setup of test data (On first login only):
```
mkdir $workdir/c3-zh-aa-bbtautau-analysis
cd $workdir/c3-zh-aa-bbtautau-analysis
ln -s /nfs_scratch/dasu/2022-03-07/Madgraph-Pythia-Delphes-Studies/c3-zh-aa-bbtautau-data/tag_1_delphes_events.root c3-zh-aa-bbtautau.root
```

### Analyzing data - digests delphes.root file and produces histos.root file:
```
source /nfs_scratch/$USER/<installation-date>/Madgraph-Pythia-Delphes-Studies/setup.sh
cd $workdir/c3-zh-aa-bbtautau-analysis
root -l -q $workdir/analyzer.C'("c3-zh-aa-bbtautau.root")'
```

### Generating data locally - result is a delphes.root file:
```
source /nfs_scratch/$USER/<installation-date>/Madgraph-Pythia-Delphes-Studies/setup.sh
cd $datadir
python $mg5dir/bin/mg5_aMC $workdir/c3-zh-aa-bbtautau.txt
```

## If you wish to use your own madgraph you need to install as below:

### First time only

If you are on the standard HEP machines you have access to /cvmfs and CentOS7 (login.hep.wisc.edu), you may use ROOT from there:

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh 
```

If you are on the newer HEP machines with /cvmfs and CentOS8 (mucol01.hep.wisc.edu), you may use ROOT from there:

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos8-gcc11-opt/setup.sh
```

Madgraph5 + Delphes installation needs to be done in a directory with plenty of space:

```
wget https://launchpad.net/mg5amcnlo/3.0/3.3.x/+download/MG5_aMC_v3.3.2.tar.gz
tar zxf MG5_aMC_v3.3.2.tar.gz 
export mg5dir=$PWD/MG5_aMC_v3_3_2/
python $mg5dir/bin/mg5_aMC
```

From within the MG5_aMC prompt execute the following. They take a long time 10 mins to finish.
Optionally, you can track the log files in the secondary login window, if you wish.

```
install lhapdf6
install pythia8
install Delphes
exit
```

### Once you have a Madgraph5 directory and tar.gz file already setup:

```
export mg5dir=<your MG5_aMC_v3_3_2 directory>
export mg5tar=<your MG5_aMC_v3_3_2 tar.gz file>
```

### Go to the dirctory where you wish to work and then install this code:

```
if [ -d '/nfs_scratch' ]; then export basedir=/nfs_scratch/$USER/`date +%Y-%m-%d`; else basedir=$PWD/`date +%Y-%m-%d`; fi
mkdir -p $basedir
git clone git@github.com:SridharaDasu/Madgraph-Pythia-Delphes-Studies.git
source $basedir/Madgraph-Pythia-Delphes-Studies/setup.sh
```

## On relogin use the base directory

For example the directory with the date of creation above, e.g., /nfs_scratch/dasu/2021-11-03/

```
source /nfs_scratch/dasu/2021-11-03/Madgraph-Pythia-Delphes-Studies/setup.sh
```

### To produce signal data (root files)

You can use appropriate *.txt files with different configurations; If you make your own signal process files, please share by making a pull request

```
cd $datadir
python $mg5dir/bin/mg5_aMC $workdir/cms-vbfh-pythia8-delphes.txt
```

### To produce data using Condor

Command to run Madgraph on the UW cluster:

```
runWiscJobs.py \
  --WorkFlow MG5Jobs \
  --Executable=runMG5JobOnWorker.sh \
  --Arguments=cms-vbfh-pythia8-delphes.txt \
  --nJobs=10 \
  --TransferInputFile=/nfs_scratch/dasu/CentOS7/MyMG5Dir.tar.gz,cms-vbfh-pythia8-delphes.txt \
  --OutputDir=/nfs_scratch/$USER \
  --HDFSProdDir None \
  --Experiment mucol \
  --MemoryRequirements 2048
```

## If you are on a Apple MacOS system

You need gcc, gfortran and root installed on your system, and you may need this additionally for Madgraph and Delphes work:

```
if [ `uname` == 'Darwin' ]; then echo export MACOSX_DEPLOYMENT_TARGET=10.15; fi
```

Once you have root, you can follow the instructions given for login machines above.

On your local system, you may want to skip Delphes first and check that Madgraph5 + Pythia are installed correctly. Here are just those steps:

```
wget https://launchpad.net/mg5amcnlo/3.0/3.3.x/+download/MG5_aMC_v3.3.2.tar.gz
tar zxf MG5_aMC_v3.3.2.tar.gz 
export mg5dir=$PWD/MG5_aMC_v3_3_2/
python $mg5dir/bin/mg5_aMC
```

follow up by:

```
install lhapdf6
install pythia8
exit
```

you can then test your setup using:

```
python $mg5dir/bin/mg5_aMC e+e-madgraph5-pythia-test.txt
```

You can then install Delphes:

```
python $mg5dir/bin/mg5_aMC
```

within mg5_aMC prompt:

```
install Delphes
```

You can then run the full chain locally on your Mac using the steps indicated above.
