import ROOT

directory= str("/nfs_scratch/dasu/2022-03-07/Madgraph-Pythia-Delphes-Studies/data/c3-z-tautau/Events/run_01/tag_1_delphes_events.root")

Tree = ROOT.TFile(directory)
print('Tree = ', Tree)
Delphes = ROOT.gROOT.FindObject("Delphes;1")
print('Delphes = ', Delphes)

for i in Delphes:
    print('i = ', i)
    x = Delphes.Jet
    print('x = ', x)
    for z in x:
        print('z = ', z)

