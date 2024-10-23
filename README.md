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
export basedir=/nfs_scratch/$USER/2024-10
mkdir -p $basedir
cd $basedir
git clone https://github.com/SridharaDasu/Madgraph-Pythia-Delphes-Studies.git
```

### Creation of test data (On first login only):
```
source /nfs_scratch/$USER/2024-10/Madgraph-Pythia-Delphes-Studies/setup.sh
cd $datadir
python $mg5dir/bin/mg5_aMC $confdir/e+e-madgraph5-pythia-delphes-test.txt
ln -s $PWD/`find e+e-madgraph5-pythia-delphes-test -name '*.root' | head -1` e+e-madgraph5-pythia-delphes-test.root
```

### Analyzing data - digests delphes.root file and produces histos.root file:
```
source /nfs_scratch/$USER/2024-10/Madgraph-Pythia-Delphes-Studies/setup.sh
root -l -q $rootdir/analyzer.C'("$datadir/e+e-madgraph5-pythia-delphes-test.root")'
```

### Generating data locally - result is a delphes.root file - replace <mg5-input-file> with your input file name (without .txt):
```
source /nfs_scratch/$USER/2024-10/Madgraph-Pythia-Delphes-Studies/setup.sh
cd $datadir
python $mg5dir/bin/mg5_aMC $confdir/<mg5-input-file>.txt
ln -s $PWD/`find <mg5-input-file> -name '*.root' | head -1` <mg5-input-file>.root
```

### To produce data using Condor (to be updated)

Command to run Madgraph on the UW cluster:

```
source /nfs_scratch/$USER/2024-10/Madgraph-Pythia-Delphes-Studies/setup.sh
cd $conddir
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

You need gcc, gfortran and root installed on your system - I used miniconda3 to get those

On your local system, you may want to skip Delphes first and check that Madgraph5 + Pythia are installed correctly. Here are just those steps after you obtain the madgraph5 file from web (https://launchpad.net/mg5amcnlo/3.0/3.6.x/+download/MG5_aMC_v3.5.6.tar.gz):

```
tar zxf MG5_aMC_v3.5.6.tar.gz
export mg5dir=$PWD/MG5_aMC_v3_5_6/
python $mg5dir/bin/mg5_aMC
```

follow up by:

```
install lhapdf6
install pythia8
exit
```

you can then install this package:

```
export basedir=<a directory on your machine>
mkdir -p $basedir
git clone https://github.com/SridharaDasu/Madgraph-Pythia-Delphes-Studies.git
```

you can then test your setup using:

```
export basedir=<your base directory on your machine>
source $basedir/Madgraph-Pythia-Delphes-Studies/setup.sh
python $mg5dir/bin/mg5_aMC e+e-madgraph5-pythia-test.txt
```

You can then install Delphes - requires mucking with the Delphes as the first attempt will fail to install:

```
python $mg5dir/bin/mg5_aMC
```

within mg5_aMC prompt:

```
install Delphes
```

Above will likely fail to install fully presently on modern M3 based Mac computers.

You need to muck with the Delphes files - (Instructions to come to get that to work soon)

You can then run the full chain locally on your Mac using the same steps as on login machines:

### Generating data locally - result is a delphes.root file - replace <mg5-input-file> with your input file name (without .txt):
```
export basedir=<your base directory on your machine>
source $basedir/Madgraph-Pythia-Delphes-Studies/setup.sh
cd $datadir
python $mg5dir/bin/mg5_aMC $confdir/<mg5-input-file>.txt
ln -s $PWD/`find <mg5-input-file> -name '*.root' | head -1` <mg5-input-file>.root
```
