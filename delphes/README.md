This is a Delphes module to make Tau objects. Copy these files to your Delphes directories, and rebuild to get TauObjects in your Delphes root file.

```
cp DelphesClasses.cc DelphesClasses.h ClassesLinkDef.h $DELPHES_HOME/classes/
cp ModulesLinkDef.h TauReconstructor.* TreeWriter.cc TreeWriter.h $DELPHES_HOME/modules/
cp delphes_card_C3.tcl $DELPHES_HOME/Cards/
cd $DELPHES_HOME
./configure
make
cp delphes_card_default.dat $mg5dir/Template/Common/Cards/
```
