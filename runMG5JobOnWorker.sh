#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh
tar zxf MyMG5Dir.tar.gz
export mg5dir=$PWD/`tar ztf MyMG5Dir.tar.gz | head -1`
for file in `grep '/nfs_scratch/dasu/CentOS7/' $mg5dir -Ilr`; do
    cat $file | sed 's|/nfs_scratch/dasu/CentOS7|'$PWD'|g' > $file
done
if [[ $1 == cms-* ]]; then
    cp $mg5dir/Delphes/cards/delphes_card_CMS.tcl $mg5dir/Template/Common/Cards/delphes_card_default.dat
elif [[ $1 == mc-* ]]; then
    tar zxf CustomMuonCollider.tar.gz
    cp delphes_card_MuonColliderDet_MuonInJetShort.tcl $mg5dir/Template/Common/Cards/delphes_card_default.dat
    cp MuonCollider/* $mg5dir/Template/Common/Cards/MuonCollider/
    cp -r $mg5dir/Delphes/cards/MuonCollider $mg5dir/Template/Common/Cards/
elif [[ $1 == c3-* ]]; then
    cp $mg5dir/Delphes/cards/delphes_card_ILD.tcl $mg5dir/Template/Common/Cards/delphes_card_default.dat
elif [[ $1 == atlas-* ]]; then
    cp $mg5dir/Delphes/cards/delphes_card_ATLAS.tcl $mg5dir/Template/Common/Cards/delphes_card_default.dat
fi
export PATH=$mg5dir/bin:$mg5dir/HEPTools/bin:$PATH
export LD_LIBRARY_PATH=$mg5dir/HEPTools/lib/:$mg5dir/HEPTools/lhapdf6_py3/lib/:$mg5dir/HEPTools/lhapdf6_py3/lib/python3.9/site-packages/:$mg5dir/HEPTools/hepmc/lib/:$mg5dir/HEPTools/pythia8//lib:$mg5dir/HEPTools/zlib/lib/:$LD_LIBRARY_PATH
export PYTHIA8DATA=$mg5dir/HEPTools/pythia8/share/Pythia8/xmldoc
export LHAPDF_DATA_PATH=$mg5dir/HEPTools/lhapdf6_py3/share/LHAPDF
source $mg5dir/Delphes/DelphesEnv.sh
export workdir=$PWD
if [ "$MyRandomNumber" == "" ]; then export MyRandomNumber=`date +"%8N"`; fi
export datadir=$MyRandomNumber
echo "Using ROOTSYS=$ROOTSYS"
echo "Using mg5dir=$mg5dir"
echo "Using PATH=$PATH"
echo "Using LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "Using PYTHIA8DATA=$PYTHIA8DATA"
echo "Using LHAPDF_DATA_PATH=$LHAPDF_DATA_PATH"
echo "Using workdir=$workdir"
echo "Using datadir=$datadir"
echo "Current directory set to $datadir"
mkdir -p $datadir
cd $datadir
cat $workdir/$1 | sed 's/set iseed 0/set iseed '"$MyRandomNumber"'/g' > mg5_configuration_$MyRandomNumber.txt
python $mg5dir/bin/mg5_aMC mg5_configuration_$MyRandomNumber.txt
cd $workdir
ls -l $datadir
tar zcf $MyRandomNumber.tar.gz $datadir
